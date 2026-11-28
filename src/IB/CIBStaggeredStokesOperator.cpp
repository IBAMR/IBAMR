// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/CIBStaggeredStokesOperator.h>
#include <ibamr/CIBStrategy.h>
#include <ibamr/IBStrategy.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/ibamr_utilities.h>

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/PETScSAMRAIVectorReal.h>

#include <tbox/MathUtilities.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

#include <petscvec.h>

#include <CoarsenSchedule.h>
#include <IntVector.h>
#include <MultiblockDataTranslator.h>
#include <PatchHierarchy.h>
#include <RefineSchedule.h>
#include <SAMRAIVectorReal.h>
#include <VariableFillPattern.h>

#include <utility>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_apply;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBStaggeredStokesOperator::CIBStaggeredStokesOperator(std::string object_name,
                                                       Pointer<CIBStrategy> cib_strategy,
                                                       bool homogeneous_bc,
                                                       Pointer<Database> input_db)
    : StaggeredStokesOperator(std::move(object_name), homogeneous_bc, input_db), d_cib_strategy(cib_strategy)
{
    // Setup Timers.
    IBAMR_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBAMR::CIBStaggeredStokesOperator::apply()"););
} // CIBStaggeredStokesOperator

void
CIBStaggeredStokesOperator::setInterpScaleFactor(const double beta)
{
    d_scale_interp = beta;
} // setInterpScaleFactor

void
CIBStaggeredStokesOperator::setSpreadScaleFactor(const double gamma)
{
    d_scale_spread = gamma;
} // setSpreadScaleFactor

void
CIBStaggeredStokesOperator::setRegularizeMobilityFactor(const double delta)
{
    d_reg_mob_factor = delta;
} // setRegularizeMobilityFactor

void
CIBStaggeredStokesOperator::setNormalizeSpreadForce(const bool normalize_force)
{
    d_normalize_spread_force = normalize_force;
} // setNormalizeSpreadForce

CIBStaggeredStokesOperator::~CIBStaggeredStokesOperator()
{
    deallocateOperatorState();
} // ~CIBStaggeredStokesOperator

void
CIBStaggeredStokesOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                    const SAMRAIVectorReal<NDIM, double>& out)
{
    StaggeredStokesOperator::initializeOperatorState(in, out);
    d_U_interp_fill = new HierarchyGhostCellInterpolation();
    d_U_interp_comp = getVelocityInterpComponent(in.getComponentDescriptorIndex(0));
    d_U_interp_fill->initializeOperatorState(d_U_interp_comp, in.getPatchHierarchy());
    d_U_rhs_interp_fill = new HierarchyGhostCellInterpolation();
    d_U_rhs_interp_comp = getVelocityRhsInterpComponent(in.getComponentDescriptorIndex(0));
    d_U_rhs_interp_fill->initializeOperatorState(d_U_rhs_interp_comp, in.getPatchHierarchy());
    return;
} // initializeOperatorState

void
CIBStaggeredStokesOperator::deallocateOperatorState()
{
    if (d_U_interp_fill)
    {
        d_U_interp_fill->deallocateOperatorState();
        d_U_interp_fill.setNull();
    }
    if (d_U_rhs_interp_fill)
    {
        d_U_rhs_interp_fill->deallocateOperatorState();
        d_U_rhs_interp_fill.setNull();
    }
    StaggeredStokesOperator::deallocateOperatorState();
    return;
} // deallocateOperatorState

void
CIBStaggeredStokesOperator::apply(Vec x, Vec y)
{
    IBAMR_TIMER_START(t_apply);

    const double half_time = 0.5 * (d_new_time + d_current_time);
    Pointer<IBStrategy> ib_method_ops = d_cib_strategy;

    // Get some vectors and unpack them.
    Vec *vx, *vy;
    VecNestGetSubVecs(x, nullptr, &vx);
    VecNestGetSubVecs(y, nullptr, &vy);
    Pointer<SAMRAIVectorReal<NDIM, double>> vx0, vy0;
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVectorRead(vx[0], &vx0);
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vy[0], &vy0);
    SAMRAIVectorReal<NDIM, double>& u_p = *vx0;
    Vec L = vx[1];
    Vec U = vx[2];
    SAMRAIVectorReal<NDIM, double>& g_f = *vy0;
    Vec V = vy[1];
    Vec F = vy[2];

    // Temporary vectors.
    Vec Vrigid;
    VecDuplicate(V, &Vrigid);

    // Get the Eulerian vector components.
    const int U_idx = u_p.getComponentDescriptorIndex(0);
    const int P_idx = u_p.getComponentDescriptorIndex(1);
    const int A_U_idx = g_f.getComponentDescriptorIndex(0);

    // Compute the action of the operator:
    // A*[u;p;U;L] :=
    //     [A_u;A_p;A_U;A_L] = [(C*I+D*L)*u + Grad P - gamma*S L; -Div u; T L;
    //                          -beta*J u + beta*T^{*} U
    //                          -beta*delta*Reg*L]

    // (a) Momentum equation and (b) divergence-free constraint, without the constraint force.
    StaggeredStokesOperator::apply(u_p, g_f);

    // Add the constraint force to the momentum equation.
    d_cib_strategy->setConstraintForce(L, half_time, -1.0 * d_scale_spread);
    ib_method_ops->spreadForce(A_U_idx, nullptr, std::vector<Pointer<RefineSchedule<NDIM>>>(), half_time);
    if (d_normalize_spread_force)
    {
        d_cib_strategy->subtractMeanConstraintForce(L, A_U_idx, -1 * d_scale_spread);
    }
    d_bc_helper->copyDataAtDirichletBoundaries(A_U_idx, U_idx);

    // Fill the same-level and physical-boundary ghost cells of u for interpolation.
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        d_U_bc_coefs, d_P_bc_coef, U_idx, P_idx, d_homogeneous_bc);
    fillVelocityGhostCells(U_idx);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);

    // (c) Rigid body velocity constraint.
    d_cib_strategy->setInterpolatedVelocityVector(V, half_time);
    d_cib_strategy->fillVelocityPhysBdryGhostValues(U_idx, d_homogeneous_bc, d_solution_time);
    ib_method_ops->interpolateVelocity(
        U_idx, std::vector<Pointer<CoarsenSchedule<NDIM>>>(), std::vector<Pointer<RefineSchedule<NDIM>>>(), half_time);

    d_cib_strategy->getInterpolatedVelocity(V, half_time, d_scale_interp);
    VecSet(Vrigid, 0.0);
    d_cib_strategy->setRigidBodyVelocity(U,
                                         Vrigid,
                                         /*only_free_dofs*/ true,
                                         /*only_imposed_dofs*/ false);

    VecScale(Vrigid, d_scale_interp);
    VecAYPX(V, -1.0, Vrigid);
    if (d_reg_mob_factor != 0.0)
    {
        d_cib_strategy->computeMobilityRegularization(Vrigid, L);
        VecAXPY(V, -1.0 * d_scale_interp * d_reg_mob_factor, Vrigid);
    }

    // (d) Force and torque constraint.
    d_cib_strategy->computeNetRigidGeneralizedForce(L,
                                                    F,
                                                    /*only_free_dofs*/ true,
                                                    /*only_imposed_dofs*/ false);
    // Delete temporary vectors.
    VecDestroy(&Vrigid);

    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVectorRead(vx[0], &vx0);
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVector(vy[0], &vy0);

    IBAMR_TIMER_STOP(t_apply);
} // apply

void
CIBStaggeredStokesOperator::modifyRhsForBcs(Vec y)
{
    const double half_time = 0.5 * (d_new_time + d_current_time);
    Pointer<IBStrategy> ib_method_ops = d_cib_strategy;

    // Get vectors corresponding to fluid and Lagrangian velocity.
    Vec* vy;
    VecNestGetSubVecs(y, nullptr, &vy);
    Pointer<SAMRAIVectorReal<NDIM, double>> vy0;
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vy[0], &vy0);
    SAMRAIVectorReal<NDIM, double>& b = *vy0;
    Vec W = vy[1];

    // Modify RHS for fluid Bcs.
    modifyRhsForBcs(b);

    // Modify RHS for the action of interpolation operator J on 0.
    if (!d_homogeneous_bc)
    {
        Vec V;
        VecDuplicate(W, &V);
        Pointer<SAMRAIVectorReal<NDIM, double>> x = b.cloneVector("");
        x->allocateVectorData();
        x->setToScalar(0.0);
        const int U_idx = x->getComponentDescriptorIndex(0);
        const int P_idx = x->getComponentDescriptorIndex(1);
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_U_bc_coefs, d_P_bc_coef, U_idx, P_idx, d_homogeneous_bc);
        d_bc_helper->enforceNormalVelocityBoundaryConditions(U_idx, P_idx, d_U_bc_coefs, d_new_time, d_homogeneous_bc);
        StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);

#if !defined(NDEBUG)
        TBOX_ASSERT(d_U_rhs_interp_fill);
#endif
        d_U_rhs_interp_fill->resetTransactionComponent(getVelocityRhsInterpComponent(U_idx));
        d_U_rhs_interp_fill->setHomogeneousBc(d_homogeneous_bc);
        d_U_rhs_interp_fill->fillData(d_solution_time);
        d_U_rhs_interp_fill->resetTransactionComponent(d_U_rhs_interp_comp);
        d_cib_strategy->fillVelocityPhysBdryGhostValues(U_idx, d_homogeneous_bc, d_solution_time);

        d_cib_strategy->setInterpolatedVelocityVector(V, half_time);
        ib_method_ops->interpolateVelocity(U_idx,
                                           std::vector<Pointer<CoarsenSchedule<NDIM>>>(),
                                           std::vector<Pointer<RefineSchedule<NDIM>>>(),
                                           half_time);

        d_cib_strategy->getInterpolatedVelocity(V, half_time, -1.0 * d_scale_interp);
        VecAXPY(W, -1.0, V);

        // Deallocate scratch data.
        free_vector_components(*x);
        VecDestroy(&V);
    }
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVector(vy[0], &vy0);
} // modifyRhsForBcs

void
CIBStaggeredStokesOperator::imposeSolBcs(Vec x)
{
    Vec* vx;
    VecNestGetSubVecs(x, nullptr, &vx);
    Pointer<SAMRAIVectorReal<NDIM, double>> vx0;
    IBTK::PETScSAMRAIVectorReal::getSAMRAIVector(vx[0], &vx0);
    SAMRAIVectorReal<NDIM, double>& u_p = *vx0;
    imposeSolBcs(u_p);
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVector(vx[0], &vx0);
} // imposeSolBcs

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

HierarchyGhostCellInterpolation::InterpolationTransactionComponent
CIBStaggeredStokesOperator::getVelocityInterpComponent(const int U_idx) const
{
    return HierarchyGhostCellInterpolation::InterpolationTransactionComponent(U_idx,
                                                                              "NONE",
                                                                              /*use_cf_bdry_interpolation*/ false,
                                                                              d_coarsen_type,
                                                                              d_bdry_extrap_type,
                                                                              /*consistent_type_2_bdry*/ false,
                                                                              d_U_bc_coefs,
                                                                              /*fill_pattern*/ nullptr,
                                                                              IBTK::enum_to_string(d_bdry_interp_type));
} // getVelocityInterpComponent

HierarchyGhostCellInterpolation::InterpolationTransactionComponent
CIBStaggeredStokesOperator::getVelocityRhsInterpComponent(const int U_idx) const
{
    // No refine operator and no Robin coefficients: the CIB strategy sets all the velocity ghost values outside the
    // physical boundaries.
    return HierarchyGhostCellInterpolation::InterpolationTransactionComponent(
        U_idx,
        "NONE",
        /*use_cf_bdry_interpolation*/ false,
        d_coarsen_type,
        d_bdry_extrap_type,
        /*consistent_type_2_bdry*/ false,
        static_cast<RobinBcCoefStrategy<NDIM>*>(nullptr),
        /*fill_pattern*/ nullptr,
        IBTK::enum_to_string(d_bdry_interp_type));
} // getVelocityRhsInterpComponent

void
CIBStaggeredStokesOperator::fillVelocityGhostCells(const int U_idx)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(d_U_interp_fill);
#endif
    d_U_interp_fill->resetTransactionComponent(getVelocityInterpComponent(U_idx));
    d_U_interp_fill->setHomogeneousBc(d_homogeneous_bc);
    d_U_interp_fill->fillData(d_solution_time);
    d_U_interp_fill->resetTransactionComponent(d_U_interp_comp);
    return;
} // fillVelocityGhostCells

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
