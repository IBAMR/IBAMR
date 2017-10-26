// Filename: VCStaggeredStokesOperator.cpp
// Created on 27 Sep 2017 by Nishant Nangia
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

#include <algorithm>
#include <ostream>
#include <stddef.h>
#include <string>
#include <vector>

#include "CellVariable.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "ibamr/VCStaggeredStokesOperator.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/SideNoCornersFillPattern.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

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

VCStaggeredStokesOperator::VCStaggeredStokesOperator(const std::string& object_name, bool homogeneous_bc)
    : StaggeredStokesOperator(object_name, homogeneous_bc),
      d_U_bdry_fill(Pointer<HierarchyGhostCellInterpolation>(NULL)),
      d_P_bdry_fill(Pointer<HierarchyGhostCellInterpolation>(NULL))
{
    // Setup Timers.
    IBAMR_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBAMR::VCStaggeredStokesOperator::apply()");
                  t_initialize_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::VCStaggeredStokesOperator::initializeOperatorState()");
                  t_deallocate_operator_state = TimerManager::getManager()->getTimer(
                      "IBAMR::VCStaggeredStokesOperator::deallocateOperatorState()"););
    return;
} // VCStaggeredStokesOperator

VCStaggeredStokesOperator::~VCStaggeredStokesOperator()
{
    deallocateOperatorState();
    delete d_default_U_bc_coef;
    d_default_U_bc_coef = NULL;
    delete d_default_P_bc_coef;
    d_default_P_bc_coef = NULL;
    return;
} // ~VCStaggeredStokesOperator

void
VCStaggeredStokesOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
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
    InterpolationTransactionComponent U_trans, P_trans;
    U_trans = InterpolationTransactionComponent(U_scratch_idx,
                                                U_idx,
                                                DATA_REFINE_TYPE,
                                                USE_CF_INTERPOLATION,
                                                DATA_COARSEN_TYPE,
                                                BDRY_EXTRAP_TYPE,
                                                CONSISTENT_TYPE_2_BDRY,
                                                d_U_bc_coefs,
                                                d_U_fill_pattern);
    P_trans = InterpolationTransactionComponent(P_idx,
                                                DATA_REFINE_TYPE,
                                                USE_CF_INTERPOLATION,
                                                DATA_COARSEN_TYPE,
                                                BDRY_EXTRAP_TYPE,
                                                CONSISTENT_TYPE_2_BDRY,
                                                d_P_bc_coef,
                                                d_P_fill_pattern);
    // d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    // d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    // StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
    //     d_U_bc_coefs, d_P_bc_coef, U_scratch_idx, P_idx, d_homogeneous_bc);
    // d_hier_bdry_fill->fillData(d_solution_time);
    // StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
    // d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);
    // d_bc_helper->enforceDivergenceFreeConditionAtBoundary(U_scratch_idx);

    // Fill velocities
    d_U_bdry_fill->resetTransactionComponents(std::vector<InterpolationTransactionComponent>(1, U_trans));
    d_U_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        d_U_bc_coefs, d_P_bc_coef, U_scratch_idx, P_idx, d_homogeneous_bc);
    d_U_bdry_fill->fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
    d_U_bdry_fill->resetTransactionComponents(std::vector<InterpolationTransactionComponent>(1, d_U_trans));

    // Enforce divergence free constraint for ghost cells at normal traction boundaries
    d_bc_helper->enforceDivergenceFreeConditionAtBoundary(U_scratch_idx);

    // Fill pressure
    d_P_bdry_fill->resetTransactionComponents(std::vector<InterpolationTransactionComponent>(1, P_trans));
    d_P_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        d_U_bc_coefs, d_P_bc_coef, U_scratch_idx, P_idx, d_homogeneous_bc);
    d_P_bdry_fill->fillData(d_solution_time);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
    d_P_bdry_fill->resetTransactionComponents(std::vector<InterpolationTransactionComponent>(1, d_P_trans));

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
                                Pointer<NodeVariable<NDIM, double> >(NULL),
#elif (NDIM == 3)
                                Pointer<EdgeVariable<NDIM, double> >(NULL),
#endif
                                U_scratch_idx,
                                U_sc_var,
                                d_no_fill,
                                d_new_time,
                                d_U_problem_coefs.cIsVariable() ? d_U_problem_coefs.getCPatchDataId() : -1,
                                Pointer<SideVariable<NDIM, double> >(NULL),
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
VCStaggeredStokesOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
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
    d_U_trans = InterpolationTransactionComponent(d_x->getComponentDescriptorIndex(0),
                                                  in.getComponentDescriptorIndex(0),
                                                  DATA_REFINE_TYPE,
                                                  USE_CF_INTERPOLATION,
                                                  DATA_COARSEN_TYPE,
                                                  BDRY_EXTRAP_TYPE,
                                                  CONSISTENT_TYPE_2_BDRY,
                                                  d_U_bc_coefs,
                                                  d_U_fill_pattern);
    d_P_trans = InterpolationTransactionComponent(in.getComponentDescriptorIndex(1),
                                                  DATA_REFINE_TYPE,
                                                  USE_CF_INTERPOLATION,
                                                  DATA_COARSEN_TYPE,
                                                  BDRY_EXTRAP_TYPE,
                                                  CONSISTENT_TYPE_2_BDRY,
                                                  d_P_bc_coef,
                                                  d_P_fill_pattern);

    // Initialize the interpolation operators.
    d_U_bdry_fill = new HierarchyGhostCellInterpolation();
    d_U_bdry_fill->initializeOperatorState(d_U_trans, d_x->getPatchHierarchy());

    d_P_bdry_fill = new HierarchyGhostCellInterpolation();
    d_P_bdry_fill->initializeOperatorState(d_P_trans, d_x->getPatchHierarchy());

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
VCStaggeredStokesOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBAMR_TIMER_START(t_deallocate_operator_state);

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    // Deallocate the interpolation operators.
    // d_hier_bdry_fill->deallocateOperatorState();
    // d_hier_bdry_fill.setNull();
    d_U_bdry_fill->deallocateOperatorState();
    d_U_bdry_fill.setNull();
    d_P_bdry_fill->deallocateOperatorState();
    d_P_bdry_fill.setNull();
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

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
