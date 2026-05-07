// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
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
#include "ibamr/VCCollocatedStokesOperator.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/ProblemSpecification.h"

#include "CellVariable.h"
#include "LocationIndexRobinBcCoefs.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "VariableFillPattern.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

#include <string>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;

// Timers.
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

VCCollocatedStokesOperator::VCCollocatedStokesOperator(const std::string& object_name,
                                                       bool homogeneous_bc,
                                                       Pointer<Database> input_db)
    : LinearOperator(object_name, homogeneous_bc),
      d_default_U_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_U_bc_coef", Pointer<Database>(nullptr))),
      d_U_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_U_bc_coef)),
      d_default_P_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_P_bc_coef", Pointer<Database>(nullptr))),
      d_P_bc_coef(d_default_P_bc_coef)
{
    if (input_db)
    {
        d_refine_type = input_db->getStringWithDefault("refine_type", d_refine_type);
        d_coarsen_type = input_db->getStringWithDefault("coarsen_type", d_coarsen_type);
        d_bdry_extrap_type = input_db->getStringWithDefault("bdry_extrap_type", d_bdry_extrap_type);
        d_bdry_interp_type = input_db->getStringWithDefault("bdry_interp_type", d_bdry_interp_type);
        d_use_cf_interpolation = input_db->getBoolWithDefault("use_cf_interpolation", d_use_cf_interpolation);
        d_consistent_type_2_bdry = input_db->getBoolWithDefault("consistent_type_2_bdry", d_consistent_type_2_bdry);
        d_eps = input_db->getDoubleWithDefault("pressure_stabilization_coef", d_eps);
    }

    // Set a default interpolation type.
    d_D_interp_type = IBTK::VC_HARMONIC_INTERP;

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBAMR::VCCollocatedStokesOperator::apply()");
                  t_initialize_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::VCCollocatedStokesOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::VCCollocatedStokesOperator::deallocateOperatorState()"););

    return;
} // VCCollocatedStokesOperator

VCCollocatedStokesOperator::~VCCollocatedStokesOperator()
{
    deallocateOperatorState();
    return;
} // ~VCCollocatedStokesOperator

void
VCCollocatedStokesOperator::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                                               RobinBcCoefStrategy<NDIM>* P_bc_coef)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(U_bc_coefs.size() == NDIM);
#endif
    d_U_bc_coefs = U_bc_coefs;
    d_P_bc_coef = P_bc_coef;
    return;
} // setPhysicalBcCoefs

void
VCCollocatedStokesOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBAMR_TIMER_START(t_apply);

    auto& vc_op_spec = static_cast<const VCViscousDilatationalOpSpec&>(*d_problem_spec);

    // Get the vector components.
    const int U_idx = x.getComponentDescriptorIndex(0);
    const int P_idx = x.getComponentDescriptorIndex(1);
    const int A_U_idx = y.getComponentDescriptorIndex(0);
    const int A_P_idx = y.getComponentDescriptorIndex(1);

    Pointer<CellVariable<NDIM, double> > U_cc_var = x.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > P_cc_var = x.getComponentVariable(1);
    Pointer<CellVariable<NDIM, double> > A_U_cc_var = y.getComponentVariable(0);
    Pointer<CellVariable<NDIM, double> > A_P_cc_var = y.getComponentVariable(1);

    // Simultaneously fill ghost cell values for all components.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(2);
    transaction_comps[0] = InterpolationTransactionComponent(U_idx,
                                                             d_refine_type,
                                                             d_use_cf_interpolation,
                                                             d_coarsen_type,
                                                             d_bdry_extrap_type,
                                                             d_consistent_type_2_bdry,
                                                             d_U_bc_coefs,
                                                             d_U_fill_pattern,
                                                             d_bdry_interp_type);
    transaction_comps[1] = InterpolationTransactionComponent(P_idx,
                                                             d_refine_type,
                                                             d_use_cf_interpolation,
                                                             d_coarsen_type,
                                                             d_bdry_extrap_type,
                                                             d_consistent_type_2_bdry,
                                                             d_P_bc_coef,
                                                             d_P_fill_pattern,
                                                             d_bdry_interp_type);
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the action of the operator:
    //
    // A*[U;P] := [A_U;A_P] = [(C*I + L1(mu) + L2(lambda))*U + Grad P; -Div (U) eps*Lp(P)]
    d_hier_math_ops->grad(A_U_idx, A_U_cc_var, 1.0, P_idx, P_cc_var, d_no_fill, d_new_time);
    // A_U += (C*I + L1(mu))*U
    double alpha = 1.0;
    double beta = 1.0;
    if (vc_op_spec.d_C_is_const)
    {
        beta = vc_op_spec.d_C_const;
    }
    d_hier_math_ops->vc_laplace(A_U_idx,
                                A_U_cc_var,
                                alpha,
                                beta,
                                vc_op_spec.d_D_idx,
                                Pointer<CellVariable<NDIM, double> >(nullptr),
                                U_idx,
                                U_cc_var,
                                d_no_fill,
                                d_new_time,
                                d_D_interp_type,
                                vc_op_spec.d_C_idx,
                                Pointer<CellVariable<NDIM, double> >(nullptr),
                                1.0,
                                A_U_idx,
                                A_U_cc_var);

    // A_U += L2(lambda)*U
    if (vc_op_spec.d_L_idx > 0)
    {
        d_hier_math_ops->vc_dilatational(A_U_idx,
                                         A_U_cc_var,
                                         1.0,
                                         vc_op_spec.d_L_idx,
                                         Pointer<CellVariable<NDIM, double> >(nullptr),
                                         U_idx,
                                         U_cc_var,
                                         d_no_fill,
                                         d_new_time,
                                         1.0,
                                         A_U_idx,
                                         A_U_cc_var);
    }

    //-Div (U)
    d_hier_math_ops->div(A_P_idx, A_P_cc_var, -1.0, U_idx, U_cc_var, d_no_fill, d_new_time);

    // Add the pressure stabilization term (Brezzi–Pitkäranta–type)
    PoissonSpecifications pressure_stabilization_spec("pressure_stab_spec");
    pressure_stabilization_spec.setCZero();
    pressure_stabilization_spec.setDConstant(d_eps);
    d_hier_math_ops->laplace(A_P_idx,
                             A_P_cc_var,
                             pressure_stabilization_spec,
                             P_idx,
                             P_cc_var,
                             d_no_fill,
                             d_new_time,
                             1.0,
                             A_P_idx,
                             A_P_cc_var);

    IBAMR_TIMER_STOP(t_apply);
    return;
} // apply

void
VCCollocatedStokesOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                                    const SAMRAIVectorReal<NDIM, double>& /*out*/)
{
    IBAMR_TIMER_START(t_initialize_operator_state);

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup the interpolation transaction information.
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
    return;
} // initializeOperatorState

void
VCCollocatedStokesOperator::deallocateOperatorState()
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
VCCollocatedStokesOperator::modifyRhsForBcs(SAMRAIVectorReal<NDIM, double>& y)
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
        apply(*x, *b);
        y.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&y, false), b);
        free_vector_components(*x);
        free_vector_components(*b);
    }
    return;
} // modifyRhsForBcs

void
VCCollocatedStokesOperator::setPressureStabilizationCoeffient(double eps)
{
    d_eps = eps;
    return;
} // setPressureStabilizationCoeffient

void
VCCollocatedStokesOperator::setDPatchDataInterpolationType(const IBTK::VCInterpType D_interp_type)
{
    d_D_interp_type = D_interp_type;
    return;
} // setDPatchDataInterpolationType

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
