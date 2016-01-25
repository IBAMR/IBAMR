// Filename: StaggeredStokesOperator.cpp
// Created on 29 Apr 2008 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#include <stddef.h>
#include <algorithm>
#include <ostream>
#include <string>
#include <vector>

#include "CellVariable.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "VariableFillPattern.h"
#include "ibamr/StaggeredStokesOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/SideNoCornersFillPattern.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;

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
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesOperator::StaggeredStokesOperator(const std::string& object_name, bool homogeneous_bc)
    : LinearOperator(object_name, homogeneous_bc),
      d_U_problem_coefs(d_object_name + "::U_problem_coefs"),
      d_default_U_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_U_bc_coef", Pointer<Database>(NULL))),
      d_U_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_U_bc_coef)),
      d_default_P_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_P_bc_coef", Pointer<Database>(NULL))),
      d_P_bc_coef(d_default_P_bc_coef),
      d_bc_helper(Pointer<StaggeredStokesPhysicalBoundaryHelper>(NULL)),
      d_U_fill_pattern(NULL),
      d_P_fill_pattern(NULL),
      d_transaction_comps(),
      d_hier_bdry_fill(Pointer<HierarchyGhostCellInterpolation>(NULL)),
      d_no_fill(Pointer<HierarchyGhostCellInterpolation>(NULL)),
      d_x(NULL),
      d_b(NULL)
{
    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions for the velocity and homogeneous Neumann
    // boundary conditions for the pressure.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        LocationIndexRobinBcCoefs<NDIM>* p_default_U_bc_coef =
            dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_U_bc_coef);
        p_default_U_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_U_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
        LocationIndexRobinBcCoefs<NDIM>* p_default_P_bc_coef =
            dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_P_bc_coef);
        p_default_P_bc_coef->setBoundarySlope(2 * d, 0.0);
        p_default_P_bc_coef->setBoundarySlope(2 * d + 1, 0.0);
    }

    // Initialize the boundary conditions objects.
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_U_bc_coef), d_default_P_bc_coef);

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBAMR::StaggeredStokesOperator::apply()");
                  t_initialize_operator_state =
                      TimerManager::getManager()->getTimer("IBAMR::StaggeredStokesOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::StaggeredStokesOperator::deallocateOperatorState()"););
    return;
} // StaggeredStokesOperator

StaggeredStokesOperator::~StaggeredStokesOperator()
{
    deallocateOperatorState();
    delete d_default_U_bc_coef;
    d_default_U_bc_coef = NULL;
    delete d_default_P_bc_coef;
    d_default_P_bc_coef = NULL;
    return;
} // ~StaggeredStokesOperator

void
StaggeredStokesOperator::setVelocityPoissonSpecifications(const PoissonSpecifications& U_problem_coefs)
{
    d_U_problem_coefs = U_problem_coefs;
    return;
} // setVelocityPoissonSpecifications

const PoissonSpecifications&
StaggeredStokesOperator::getVelocityPoissonSpecifications() const
{
    return d_U_problem_coefs;
} // getVelocityPoissonSpecifications

void
StaggeredStokesOperator::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                                            RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(U_bc_coefs.size() == NDIM);
#endif
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (U_bc_coefs[d])
        {
            d_U_bc_coefs[d] = U_bc_coefs[d];
        }
        else
        {
            d_U_bc_coefs[d] = d_default_U_bc_coef;
        }
    }

    if (P_bc_coef)
    {
        d_P_bc_coef = P_bc_coef;
    }
    else
    {
        d_P_bc_coef = d_default_P_bc_coef;
    }
    return;
} // setPhysicalBcCoefs

void
StaggeredStokesOperator::setPhysicalBoundaryHelper(Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(bc_helper);
#endif
    d_bc_helper = bc_helper;
    return;
} // setPhysicalBoundaryHelper

void
StaggeredStokesOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBAMR_TIMER_START(t_apply);

    // Allocate scratch data.
    d_x->allocateVectorData();

    // Get the vector components.
    const int U_idx = x.getComponentDescriptorIndex(0);
    const int P_idx = x.getComponentDescriptorIndex(1);
    const int A_U_idx = y.getComponentDescriptorIndex(0);
    const int A_P_idx = y.getComponentDescriptorIndex(1);
    const int U_scratch_idx = d_x->getComponentDescriptorIndex(0);

    Pointer<SideVariable<NDIM, double> > U_sc_var = x.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > P_cc_var = x.getComponentVariable(1);
    Pointer<SideVariable<NDIM, double> > A_U_sc_var = y.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > A_P_cc_var = y.getComponentVariable(1);

    // Simultaneously fill ghost cell values for all components.
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(2);
    transaction_comps[0] = InterpolationTransactionComponent(U_scratch_idx,
                                                             U_idx,
                                                             DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_U_bc_coefs,
                                                             d_U_fill_pattern);
    transaction_comps[1] = InterpolationTransactionComponent(P_idx,
                                                             DATA_REFINE_TYPE,
                                                             USE_CF_INTERPOLATION,
                                                             DATA_COARSEN_TYPE,
                                                             BDRY_EXTRAP_TYPE,
                                                             CONSISTENT_TYPE_2_BDRY,
                                                             d_P_bc_coef,
                                                             d_P_fill_pattern);
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        d_U_bc_coefs, d_P_bc_coef, U_scratch_idx, P_idx, d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the action of the operator:
    //
    // A*[U;P] := [A_U;A_P] = [(C*I+D*L)*U + Grad P; -Div U]
    d_hier_math_ops->grad(A_U_idx,
                          A_U_sc_var,
                          /*cf_bdry_synch*/ false,
                          1.0,
                          P_idx,
                          P_cc_var,
                          d_no_fill,
                          d_new_time);
    d_hier_math_ops->laplace(A_U_idx,
                             A_U_sc_var,
                             d_U_problem_coefs,
                             U_scratch_idx,
                             U_sc_var,
                             d_no_fill,
                             d_new_time,
                             1.0,
                             A_U_idx,
                             A_U_sc_var);
    d_hier_math_ops->div(A_P_idx,
                         A_P_cc_var,
                         -1.0,
                         U_scratch_idx,
                         U_sc_var,
                         d_no_fill,
                         d_new_time,
                         /*cf_bdry_synch*/ true);
    d_bc_helper->copyDataAtDirichletBoundaries(A_U_idx, U_scratch_idx);

    // Deallocate scratch data.
    d_x->deallocateVectorData();

    IBAMR_TIMER_STOP(t_apply);
    return;
} // apply

void
StaggeredStokesOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                 const SAMRAIVectorReal<NDIM, double>& out)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup solution and rhs vectors.
    d_x = in.cloneVector(in.getName());
    d_b = out.cloneVector(out.getName());

    // Setup the interpolation transaction information.
    d_U_fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);
    d_P_fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    d_transaction_comps.resize(2);
    d_transaction_comps[0] = InterpolationTransactionComponent(d_x->getComponentDescriptorIndex(0),
                                                               in.getComponentDescriptorIndex(0),
                                                               DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_U_bc_coefs,
                                                               d_U_fill_pattern);
    d_transaction_comps[1] = InterpolationTransactionComponent(in.getComponentDescriptorIndex(1),
                                                               DATA_REFINE_TYPE,
                                                               USE_CF_INTERPOLATION,
                                                               DATA_COARSEN_TYPE,
                                                               BDRY_EXTRAP_TYPE,
                                                               CONSISTENT_TYPE_2_BDRY,
                                                               d_P_bc_coef,
                                                               d_P_fill_pattern);

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_x->getPatchHierarchy());

    // Initialize hierarchy math ops object.
    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops = new HierarchyMathOps(d_object_name + "::HierarchyMathOps",
                                               in.getPatchHierarchy(),
                                               in.getCoarsestLevelNumber(),
                                               in.getFinestLevelNumber());
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

void
StaggeredStokesOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    // Deallocate the interpolation operators.
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.setNull();
    d_transaction_comps.clear();
    d_U_fill_pattern.setNull();
    d_P_fill_pattern.setNull();

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

void
StaggeredStokesOperator::modifyRhsForBcs(SAMRAIVectorReal<NDIM, double>& y)
{
    if (!d_homogeneous_bc)
    {
        // Set y := y - A*0, i.e., shift the right-hand-side vector to account for
        // inhomogeneous boundary conditions.
        Pointer<SAMRAIVectorReal<NDIM, double> > x = y.cloneVector("");
        Pointer<SAMRAIVectorReal<NDIM, double> > b = y.cloneVector("");
        x->allocateVectorData();
        b->allocateVectorData();
        x->setToScalar(0.0);
        const int U_idx = x->getComponentDescriptorIndex(0);
        const int P_idx = x->getComponentDescriptorIndex(1);
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_U_bc_coefs, d_P_bc_coef, U_idx, P_idx, d_homogeneous_bc);
        d_bc_helper->enforceNormalVelocityBoundaryConditions(U_idx, P_idx, d_U_bc_coefs, d_new_time, d_homogeneous_bc);
        StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
        apply(*x, *b);
        y.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&y, false), b);
        x->freeVectorComponents();
        b->freeVectorComponents();
    }
    const bool homogeneous_bc = true;
    const int U_idx = y.getComponentDescriptorIndex(0);
    const int P_idx = y.getComponentDescriptorIndex(1);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(d_U_bc_coefs, d_P_bc_coef, U_idx, P_idx, homogeneous_bc);
    d_bc_helper->enforceNormalVelocityBoundaryConditions(U_idx, P_idx, d_U_bc_coefs, d_new_time, homogeneous_bc);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
    return;
} // modifyRhsForBcs

void
StaggeredStokesOperator::imposeSolBcs(SAMRAIVectorReal<NDIM, double>& u)
{
    const int U_idx = u.getComponentDescriptorIndex(0);
    const int P_idx = u.getComponentDescriptorIndex(1);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        d_U_bc_coefs, d_P_bc_coef, U_idx, P_idx, d_homogeneous_bc);
    d_bc_helper->enforceNormalVelocityBoundaryConditions(U_idx, P_idx, d_U_bc_coefs, d_new_time, d_homogeneous_bc);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
    return;
} // imposeSolBcs

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
