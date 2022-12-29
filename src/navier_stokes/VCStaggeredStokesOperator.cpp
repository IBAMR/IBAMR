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

#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/VCStaggeredStokesOperator.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"

#include "CellVariable.h"
#include "EdgeVariable.h"
#include "NodeVariable.h"
#include "PoissonSpecifications.h"
#include "RobinBcCoefStrategy.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
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
// Timers.
static Timer* t_apply;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

VCStaggeredStokesOperator::VCStaggeredStokesOperator(const std::string& object_name,
                                                     bool homogeneous_bc,
                                                     Pointer<Database> input_db)
    : StaggeredStokesOperator(object_name, homogeneous_bc, input_db)
{
    // Setup Timers.
    IBAMR_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBAMR::VCStaggeredStokesOperator::apply()"););

    // Set a default interpolation type.
    d_D_interp_type = IBTK::VC_HARMONIC_INTERP;

    return;
} // VCStaggeredStokesOperator

VCStaggeredStokesOperator::~VCStaggeredStokesOperator()
{
    deallocateOperatorState();
    delete d_default_U_bc_coef;
    d_default_U_bc_coef = nullptr;
    delete d_default_P_bc_coef;
    d_default_P_bc_coef = nullptr;
    return;
} // ~VCStaggeredStokesOperator

void
VCStaggeredStokesOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBAMR_TIMER_START(t_apply);

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
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps(2);
    transaction_comps[0] = InterpolationTransactionComponent(U_scratch_idx,
                                                             U_idx,
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
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        d_U_bc_coefs, d_P_bc_coef, U_scratch_idx, P_idx, d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);
    d_bc_helper->enforceDivergenceFreeConditionAtBoundary(U_scratch_idx);

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
    // A_U += (C*I*L(D))*U
    double alpha = 1.0;
    double beta = 1.0;
    if (d_U_problem_coefs.cIsZero() || d_U_problem_coefs.cIsConstant())
    {
        beta = d_U_problem_coefs.cIsZero() ? 0.0 : d_U_problem_coefs.getCConstant();
    }
    d_hier_math_ops->vc_laplace(A_U_idx,
                                A_U_sc_var,
                                alpha,
                                beta,
                                d_U_problem_coefs.getDPatchDataId(),
#if (NDIM == 2)
                                Pointer<NodeVariable<NDIM, double> >(nullptr),
#elif (NDIM == 3)
                                Pointer<EdgeVariable<NDIM, double> >(nullptr),
#endif
                                U_scratch_idx,
                                U_sc_var,
                                d_no_fill,
                                d_new_time,
                                d_D_interp_type,
                                d_U_problem_coefs.cIsVariable() ? d_U_problem_coefs.getCPatchDataId() : -1,
                                Pointer<SideVariable<NDIM, double> >(nullptr),
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

    IBAMR_TIMER_STOP(t_apply);
    return;
} // apply

void
VCStaggeredStokesOperator::setDPatchDataInterpolationType(const IBTK::VCInterpType D_interp_type)
{
    d_D_interp_type = D_interp_type;
    return;
} // setDPatchDataInterpolationType

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
