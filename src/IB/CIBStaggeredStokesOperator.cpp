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
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>

#include <tbox/MathUtilities.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>

#include <petscvec.h>

#include <CoarsenSchedule.h>
#include <IntVector.h>
#include <MultiblockDataTranslator.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
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
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

CIBStaggeredStokesOperator::CIBStaggeredStokesOperator(std::string object_name,
                                                       Pointer<CIBStrategy> cib_strategy,
                                                       bool homogeneous_bc,
                                                       Pointer<Database> input_db)
    : StaggeredStokesOperator(std::move(object_name), homogeneous_bc, input_db), d_cib_strategy(cib_strategy)
{
    // Setup Timers.
    IBAMR_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBAMR::CIBStaggeredStokesOperator::apply()");
                  t_initialize_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::CIBStaggeredStokesOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::CIBStaggeredStokesOperator::deallocateOperatorState()"););
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

    // Compute the action of the operator:
    // A*[u;p;U;L] :=
    //     [A_u;A_p;A_U;A_L] = [(C*I+D*L)*u + Grad P - gamma*S L; -Div u; T L;
    //                          -beta*J u + beta*T^{*} U
    //                          -beta*delta*Reg*L]

    // (a) Momentum equation and (b) divergence-free constraint.
    StaggeredStokesOperator::apply(u_p, g_f);
    const int U_idx = u_p.getComponentDescriptorIndex(0);
    const int A_U_idx = g_f.getComponentDescriptorIndex(0);

    d_cib_strategy->setConstraintForce(L, half_time, -1.0 * d_scale_spread);
    ib_method_ops->spreadForce(A_U_idx, nullptr, std::vector<Pointer<RefineSchedule<NDIM>>>(), half_time);
    if (d_normalize_spread_force)
    {
        d_cib_strategy->subtractMeanConstraintForce(L, A_U_idx, -1 * d_scale_spread);
    }

    d_bc_helper->copyDataAtDirichletBoundaries(A_U_idx, U_idx);

    // (c) Rigid body velocity constraint.
    interpolateVelocity(u_p, V, half_time, d_scale_interp);
    VecSet(Vrigid, 0.0);
    d_cib_strategy->setRigidBodyVelocity(U,
                                         Vrigid,
                                         /*only_free_dofs*/ true,
                                         /*only_imposed_dofs*/ false);

    VecScale(Vrigid, d_scale_interp);
    VecAYPX(V, -1.0, Vrigid);
    if (!IBTK::abs_equal_eps(d_reg_mob_factor, 0.0))
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
CIBStaggeredStokesOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                    const SAMRAIVectorReal<NDIM, double>& /*out*/)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup the interpolation transaction information.
    //
    // NOTE: This will ensure that all ghost cell values are filled in communication.  The default communication
    // patterns setup by StaggeredStokesOperator will omit corners and (in 3D) edges.
    d_U_fill_pattern = nullptr;
    d_P_fill_pattern = nullptr;
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_transaction_comps.resize(2);
    d_transaction_comps[0] = InterpolationTransactionComponent(in.getComponentDescriptorIndex(0),
                                                               d_refine_type,
                                                               d_use_cf_interpolation,
                                                               d_coarsen_type,
                                                               d_bdry_extrap_type,
                                                               d_consistent_type_2_bdry,
                                                               d_U_bc_coefs,
                                                               d_U_fill_pattern,
                                                               d_bdry_interp_type);
    d_transaction_comps[1] = InterpolationTransactionComponent(in.getComponentDescriptorIndex(1),
                                                               d_refine_type,
                                                               d_use_cf_interpolation,
                                                               d_coarsen_type,
                                                               d_bdry_extrap_type,
                                                               d_consistent_type_2_bdry,
                                                               d_P_bc_coef,
                                                               d_P_fill_pattern,
                                                               d_bdry_interp_type);

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, in.getPatchHierarchy());

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
} // initializeOperatorState

void
CIBStaggeredStokesOperator::deallocateOperatorState()
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

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    IBAMR_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

void
CIBStaggeredStokesOperator::modifyRhsForBcs(Vec y)
{
    const double half_time = 0.5 * (d_new_time + d_current_time);

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

        fillGhostCellValues(*x);
        interpolateVelocity(*x, V, half_time, -1.0 * d_scale_interp);
        VecAXPY(W, -1.0, V);

        // Deallocate scratch data.
        free_vector_components(*x);
        VecDestroy(&V);
    }
    IBTK::PETScSAMRAIVectorReal::restoreSAMRAIVector(vy[0], &vy0);
} // modifyRhsForBcs

void
CIBStaggeredStokesOperator::interpolateVelocity(SAMRAIVectorReal<NDIM, double>& u_p,
                                                Vec V,
                                                const double data_time,
                                                const double scale)
{
    // Spreading uses the IB velocity boundary conditions at physical boundaries,
    // so interpolate with the same ghost values to keep the two adjoint.
    const int U_idx = u_p.getComponentDescriptorIndex(0);
    RobinPhysBdryPatchStrategy* const u_phys_bdry_op = d_cib_strategy->getVelocityPhysBdryOp();
    if (u_phys_bdry_op)
    {
        u_phys_bdry_op->setPatchDataIndex(U_idx);
        u_phys_bdry_op->setHomogeneousBc(d_homogeneous_bc);
        Pointer<PatchHierarchy<NDIM>> hierarchy = u_p.getPatchHierarchy();
        for (int ln = u_p.getCoarsestLevelNumber(); ln <= u_p.getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                u_phys_bdry_op->setPhysicalBoundaryConditions(
                    *patch, d_solution_time, patch->getPatchData(U_idx)->getGhostCellWidth());
            }
        }
    }

    Pointer<IBStrategy> ib_method_ops = d_cib_strategy;
    d_cib_strategy->setInterpolatedVelocityVector(V, data_time);
    ib_method_ops->interpolateVelocity(
        U_idx, std::vector<Pointer<CoarsenSchedule<NDIM>>>(), std::vector<Pointer<RefineSchedule<NDIM>>>(), data_time);
    d_cib_strategy->getInterpolatedVelocity(V, data_time, scale);
    return;
} // interpolateVelocity

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

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
