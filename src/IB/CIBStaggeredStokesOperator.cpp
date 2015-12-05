// Filename: CIBStaggeredStokesOperator.cpp
// Created on 31 Oct 2014 by Amneet Bhalla
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <algorithm>
#include <map>

#include "ibamr/CIBStaggeredStokesOperator.h"
#include "ibamr/CIBStrategy.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/PETScMultiVec.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/HierarchyMathOps.h"
#include "CellVariable.h"
#include "PatchHierarchy.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "LocationIndexRobinBcCoefs.h"
#include "tbox/Database.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Types of refining and coarsening to perform prior to setting coarse-fine
// boundary and physical boundary ghost cell values.
static const std::string DATA_REFINE_TYPE = "NONE";
static const bool USE_CF_INTERPOLATION = true;
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Timers.
static Timer* t_apply;
static Timer* t_apply_vec;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBStaggeredStokesOperator::CIBStaggeredStokesOperator(const std::string& object_name,
                                                       Pointer<CIBStrategy> cib_strategy,
                                                       bool homogeneous_bc)
    : LinearOperator(object_name, homogeneous_bc), d_cib_strategy(cib_strategy),
      d_num_rigid_parts(d_cib_strategy->getNumberOfRigidStructures()),
      d_u_problem_coefs(d_object_name + "::u_problem_coefs"),
      d_default_u_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_u_bc_coef", Pointer<Database>(NULL))),
      d_u_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_u_bc_coef)),
      d_default_p_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_p_bc_coef", Pointer<Database>(NULL))),
      d_p_bc_coef(d_default_p_bc_coef), d_bc_helper(Pointer<StaggeredStokesPhysicalBoundaryHelper>(NULL)),
      d_u_fill_pattern(NULL), d_p_fill_pattern(NULL), d_transaction_comps(),
      d_hier_bdry_fill(Pointer<HierarchyGhostCellInterpolation>(NULL)),
      d_no_fill(Pointer<HierarchyGhostCellInterpolation>(NULL)), d_x(NULL), d_b(NULL), d_scale_interp(1.0),
      d_scale_spread(1.0), d_reg_mob_factor(1.0), d_normalize_spread_force(false)
{

    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions for the velocity and homogeneous Neumann
    // boundary conditions for the pressure.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        LocationIndexRobinBcCoefs<NDIM>* p_default_U_bc_coef =
            dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_u_bc_coef);
        p_default_U_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_U_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
        LocationIndexRobinBcCoefs<NDIM>* p_default_P_bc_coef =
            dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_p_bc_coef);
        p_default_P_bc_coef->setBoundarySlope(2 * d, 0.0);
        p_default_P_bc_coef->setBoundarySlope(2 * d + 1, 0.0);
    }

    // Initialize the boundary conditions objects.
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_u_bc_coef), d_default_p_bc_coef);

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_apply = TimerManager::getManager()->getTimer("IBAMR::CIBStaggeredStokesOperator::apply(SVR,SVR)");
        t_apply_vec = TimerManager::getManager()->getTimer("IBAMR::CIBStaggeredStokesOperator::apply(Vec,Vec)");
        t_initialize_operator_state =
            TimerManager::getManager()->getTimer("IBAMR::CIBStaggeredStokesOperator::initializeOperatorState()");
        t_deallocate_operator_state =
            TimerManager::getManager()->getTimer("IBAMR::CIBStaggeredStokesOperator::deallocateOperatorState()"););
    return;
} // CIBStaggeredStokesOperator

void CIBStaggeredStokesOperator::setInterpScaleFactor(const double beta)
{
    d_scale_interp = beta;
    return;
} // setInterpScaleFactor

void CIBStaggeredStokesOperator::setSpreadScaleFactor(const double gamma)
{
    d_scale_spread = gamma;
    return;
} // setSpreadScaleFactor

void CIBStaggeredStokesOperator::setRegularizeMobilityFactor(const double delta)
{
    d_reg_mob_factor = delta;
    return;
} // setRegularizeMobilityFactor

void CIBStaggeredStokesOperator::setNormalizeSpreadForce(const bool normalize_force)
{
    d_normalize_spread_force = normalize_force;
    return;
} // setNormalizeSpreadForce

CIBStaggeredStokesOperator::~CIBStaggeredStokesOperator()
{
    deallocateOperatorState();
    delete d_default_u_bc_coef;
    d_default_u_bc_coef = NULL;
    delete d_default_p_bc_coef;
    d_default_p_bc_coef = NULL;
    return;
} // ~CIBStaggeredStokesOperator

void CIBStaggeredStokesOperator::setVelocityPoissonSpecifications(const PoissonSpecifications& u_problem_coefs)
{
    d_u_problem_coefs = u_problem_coefs;
    return;
} // setVelocityPoissonSpecifications

void CIBStaggeredStokesOperator::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                                    RobinBcCoefStrategy<NDIM>* p_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(u_bc_coefs.size() == NDIM);
#endif
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (u_bc_coefs[d])
        {
            d_u_bc_coefs[d] = u_bc_coefs[d];
        }
        else
        {
            d_u_bc_coefs[d] = d_default_u_bc_coef;
        }
    }

    if (p_bc_coef)
    {
        d_p_bc_coef = p_bc_coef;
    }
    else
    {
        d_p_bc_coef = d_default_p_bc_coef;
    }
    return;
} // setPhysicalBcCoefs

void CIBStaggeredStokesOperator::setPhysicalBoundaryHelper(Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_helper);
#endif
    d_bc_helper = bc_helper;
    return;
} // setPhysicalBoundaryHelper

void CIBStaggeredStokesOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBAMR_TIMER_START(t_apply);
    const double half_time = 0.5 * (d_new_time + d_current_time);

    // Get the vector components.
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    const int A_u_idx = y.getComponentDescriptorIndex(0);
    const int A_p_idx = y.getComponentDescriptorIndex(1);
    const int u_scratch_idx = d_x->getComponentDescriptorIndex(0);

    Pointer<SideVariable<NDIM, double> > u_sc_var = x.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > p_cc_var = x.getComponentVariable(1);
    Pointer<SideVariable<NDIM, double> > A_u_sc_var = y.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > A_p_cc_var = y.getComponentVariable(1);

    // Simultaneously fill ghost cell values for all components.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(2);
    transaction_comps[0] = InterpolationTransactionComponent(u_scratch_idx, u_idx, DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY, d_u_bc_coefs, d_u_fill_pattern);
    transaction_comps[1] =
        InterpolationTransactionComponent(p_idx, DATA_REFINE_TYPE, USE_CF_INTERPOLATION, DATA_COARSEN_TYPE,
                                          BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_p_bc_coef, d_p_fill_pattern);
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_u_bc_coefs, d_p_bc_coef, u_scratch_idx, p_idx,
                                                              d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_u_bc_coefs, d_p_bc_coef);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the action of the operator:
    //
    // A*[u;p] := [A_u;A_p] = [(C*I+D*L)*u + Grad p; -Div u]
    d_hier_math_ops->grad(A_u_idx, A_u_sc_var, /*cf_bdry_synch*/ false, 1.0, p_idx, p_cc_var, d_no_fill, half_time);
    d_hier_math_ops->laplace(A_u_idx, A_u_sc_var, d_u_problem_coefs, u_scratch_idx, u_sc_var, d_no_fill, half_time, 1.0,
                             A_u_idx, A_u_sc_var);
    d_hier_math_ops->div(A_p_idx, A_p_cc_var, -1.0, u_scratch_idx, u_sc_var, d_no_fill, half_time,
                         /*cf_bdry_synch*/ true);
    d_bc_helper->copyDataAtDirichletBoundaries(A_u_idx, u_scratch_idx);

    IBAMR_TIMER_STOP(t_apply);
    return;
} // apply

void CIBStaggeredStokesOperator::apply(Vec x, Vec y)
{
#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    clock_t end_t = 0, start_med = 0;
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif

    IBAMR_TIMER_START(t_apply_vec);
    const double half_time = 0.5 * (d_new_time + d_current_time);
    Pointer<IBStrategy> ib_method_ops = d_cib_strategy;

    // Get some vectors and unpack them.
    Vec *vx, *vy;
    VecMultiVecGetSubVecs(x, &vx);
    VecMultiVecGetSubVecs(y, &vy);
    SAMRAIVectorReal<NDIM, double>& u_p = *PETScSAMRAIVectorReal::getSAMRAIVector(vx[0]);
    Vec L = vx[1];
    Vec U = vx[2];
    SAMRAIVectorReal<NDIM, double>& g_f = *PETScSAMRAIVectorReal::getSAMRAIVector(vy[0]);
    Vec V = vy[1];
    Vec F = vy[2];

    // VecView(U, PETSC_VIEWER_STDOUT_WORLD);
    // Temporary vectors.
    Vec Vrigid;
    VecDuplicate(V, &Vrigid);

    // Get the Eulerian vector components.
    const int u_idx = u_p.getComponentDescriptorIndex(0);
    const int p_idx = u_p.getComponentDescriptorIndex(1);
    const int A_u_idx = g_f.getComponentDescriptorIndex(0);
    const int A_p_idx = g_f.getComponentDescriptorIndex(1);
    const int u_scratch_idx = d_x->getComponentDescriptorIndex(0);

    Pointer<SideVariable<NDIM, double> > u_sc_var = u_p.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > p_cc_var = u_p.getComponentVariable(1);
    Pointer<SideVariable<NDIM, double> > A_u_sc_var = g_f.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > A_p_cc_var = g_f.getComponentVariable(1);
#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
        end_t = clock();
        pout << std::setprecision(4) << "      StokesOperator:Allocations CPU time taken for the time step is:"
             << double(end_t - start_med) / double(CLOCKS_PER_SEC) << std::endl;
        ;
    }
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif

    // Simultaneously fill ghost cell values for u and p.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(2);
    transaction_comps[0] = InterpolationTransactionComponent(u_scratch_idx, u_idx, DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY, d_u_bc_coefs, d_u_fill_pattern);
    transaction_comps[1] =
        InterpolationTransactionComponent(p_idx, DATA_REFINE_TYPE, USE_CF_INTERPOLATION, DATA_COARSEN_TYPE,
                                          BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_p_bc_coef, d_p_fill_pattern);
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_u_bc_coefs, d_p_bc_coef, u_scratch_idx, p_idx,
                                                              d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_u_bc_coefs, d_p_bc_coef);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
        end_t = clock();
        pout << std::setprecision(4)
             << "      StokesOperator:Interpolatiions + setups CPU time taken for the time step is:"
             << double(end_t - start_med) / double(CLOCKS_PER_SEC) << std::endl;
        ;
    }
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif

    // Compute the action of the operator:
    // A*[u;p;U;L] := [A_u;A_p;A_U;A_L] = [(C*I+D*L)*u + Grad P - gamma*S L; -Div u; T L;
    //                                     -beta*J u + beta*T^{*} U -beta*delta*Reg*L]

    // (a) Momentum equation.
    d_hier_math_ops->grad(A_u_idx, A_u_sc_var, /*cf_bdry_synch*/ false, 1.0, p_idx, p_cc_var, d_no_fill, half_time);
    d_hier_math_ops->laplace(A_u_idx, A_u_sc_var, d_u_problem_coefs, u_scratch_idx, u_sc_var, d_no_fill, half_time, 1.0,
                             A_u_idx, A_u_sc_var);
#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
        end_t = clock();
        pout << std::setprecision(4)
             << "      StokesOperator:grad + laplace operator, CPU time taken for the time step is:"
             << double(end_t - start_med) / double(CLOCKS_PER_SEC) << std::endl;
        ;
    }
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif

    d_cib_strategy->setConstraintForce(L, half_time, -1.0 * d_scale_spread);
    ib_method_ops->spreadForce(A_u_idx, NULL, std::vector<Pointer<RefineSchedule<NDIM> > >(), half_time);
    if (d_normalize_spread_force)
    {
        d_cib_strategy->subtractMeanConstraintForce(L, A_u_idx, -1 * d_scale_spread);
    }
#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
        end_t = clock();
        pout << std::setprecision(4)
             << "      StokesOperator:Force spreading + normalizations CPU time taken for the time step is:"
             << double(end_t - start_med) / double(CLOCKS_PER_SEC) << std::endl;
        ;
    }
#endif
    // (b) Divergence-free constraint.
    d_hier_math_ops->div(A_p_idx, A_p_cc_var, -1.0, u_scratch_idx, u_sc_var, d_no_fill, half_time,
                         /*cf_bdry_synch*/ true);
    d_bc_helper->copyDataAtDirichletBoundaries(A_u_idx, u_scratch_idx);

// (c) Rigid body velocity constraint.

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif

    std::vector<InterpolationTransactionComponent> u_transaction_comps(1);
    u_transaction_comps[0] =
        InterpolationTransactionComponent(u_idx, DATA_REFINE_TYPE, USE_CF_INTERPOLATION, DATA_COARSEN_TYPE,
                                          BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_u_bc_coefs, d_u_fill_pattern);
    Pointer<HierarchyGhostCellInterpolation> u_bdry_fill = new IBTK::HierarchyGhostCellInterpolation();
    u_bdry_fill->initializeOperatorState(u_transaction_comps, u_p.getPatchHierarchy());
    u_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    u_bdry_fill->fillData(d_solution_time);

    d_cib_strategy->setInterpolatedVelocityVector(V, half_time);
    ib_method_ops->interpolateVelocity(u_idx, std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                                       std::vector<Pointer<RefineSchedule<NDIM> > >(), half_time);
    d_cib_strategy->getInterpolatedVelocity(V, half_time, d_scale_interp);
    VecSet(Vrigid, 0.0);

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
        end_t = clock();
        pout << std::setprecision(4) << "      StokesOperator:Interpolate Velocity CPU time taken for the time step is:"
             << double(end_t - start_med) / double(CLOCKS_PER_SEC) << std::endl;
        ;
    }
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif

    d_cib_strategy->setRigidBodyVelocity(U, Vrigid, /*only_free_dofs*/ true,
                                         /*only_imposed_dofs*/ false);

    VecScale(Vrigid, d_scale_interp);
    VecAYPX(V, -1.0, Vrigid);

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
        end_t = clock();
        pout << std::setprecision(4) << "      StokesOperator:setRigidVelocity, CPU time taken for the time step is:"
             << double(end_t - start_med) / double(CLOCKS_PER_SEC) << std::endl;
        ;
    }
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif

    if (!MathUtilities<double>::equalEps(d_reg_mob_factor, 0.0))
    {
        d_cib_strategy->computeMobilityRegularization(Vrigid, L);
        VecAXPY(V, -1.0 * d_scale_interp * d_reg_mob_factor, Vrigid);
    }

#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
        end_t = clock();
        pout << std::setprecision(4)
             << "      StokesOperator: Mobility Regulization, CPU time taken for the time step is:"
             << double(end_t - start_med) / double(CLOCKS_PER_SEC) << std::endl;
        ;
    }
    if (SAMRAI_MPI::getRank() == 0) start_med = clock();
#endif

    // (d) Force and torque constraint.
    d_cib_strategy->computeNetRigidGeneralizedForce(L, F, /*only_free_dofs*/ true,
                                                    /*only_imposed_dofs*/ false);
#ifdef TIME_REPORT
    SAMRAI_MPI::barrier();
    if (SAMRAI_MPI::getRank() == 0)
    {
        end_t = clock();
        pout << std::setprecision(4) << "      StokesOperator: K^T*lambda, CPU time taken for the time step is:"
             << double(end_t - start_med) / double(CLOCKS_PER_SEC) << std::endl;
    }
#endif

    // Delete temporary vectors.
    VecDestroy(&Vrigid);

    IBAMR_TIMER_STOP(t_apply_vec);
    return;
} // apply

void CIBStaggeredStokesOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                         const SAMRAIVectorReal<NDIM, double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_x = in.cloneVector(in.getName());
    d_b = out.cloneVector(out.getName());
    d_x->allocateVectorData();

    // Setup the interpolation transaction information.
    d_u_fill_pattern = NULL;
    d_p_fill_pattern = NULL;
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    d_transaction_comps.resize(2);
    d_transaction_comps[0] = InterpolationTransactionComponent(
        d_x->getComponentDescriptorIndex(0), in.getComponentDescriptorIndex(0), DATA_REFINE_TYPE, USE_CF_INTERPOLATION,
        DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_u_bc_coefs, d_u_fill_pattern);
    d_transaction_comps[1] = InterpolationTransactionComponent(
        in.getComponentDescriptorIndex(1), DATA_REFINE_TYPE, USE_CF_INTERPOLATION, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE,
        CONSISTENT_TYPE_2_BDRY, d_p_bc_coef, d_p_fill_pattern);

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_x->getPatchHierarchy());

    // Initialize hierarchy math ops object.
    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops = new HierarchyMathOps(d_object_name + "::HierarchyMathOps", in.getPatchHierarchy(),
                                               in.getCoarsestLevelNumber(), in.getFinestLevelNumber());
    }
#if !defined(NDEBUG)
    else
    {
        TBOX_ASSERT(d_hier_math_ops);
    }
#endif

    // Indicate the operator is initialized.
    d_is_initialized = true;

    IBAMR_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void CIBStaggeredStokesOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    // Deallocate the interpolation operators.
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.setNull();
    d_transaction_comps.clear();
    d_u_fill_pattern.setNull();
    d_p_fill_pattern.setNull();

    // Delete the solution and rhs vectors.
    d_x->resetLevels(d_x->getCoarsestLevelNumber(),
                     std::min(d_x->getFinestLevelNumber(), d_x->getPatchHierarchy()->getFinestLevelNumber()));
    d_x->freeVectorComponents();

    d_b->resetLevels(d_b->getCoarsestLevelNumber(),
                     std::min(d_b->getFinestLevelNumber(), d_b->getPatchHierarchy()->getFinestLevelNumber()));
    d_b->freeVectorComponents();

    d_x.setNull();
    d_b.setNull();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
