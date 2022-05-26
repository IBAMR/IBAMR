// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2022 by the IBAMR developers
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

#include "ibamr/StaggeredStokesOperator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_utilities.h"

#include "ibtk/CellNoCornersFillPattern.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LinearOperator.h"
#include "ibtk/SideNoCornersFillPattern.h"

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
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "tbox/Timer.h"
#include "tbox/TimerManager.h"

#include <algorithm>
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
static const int SIDEG = 1;

// Timers.
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

StaggeredStokesOperator::StaggeredStokesOperator(const std::string& object_name,
                                                 bool homogeneous_bc,
                                                 Pointer<Database> input_db)
    : LinearOperator(object_name, homogeneous_bc),
      d_U_problem_coefs(d_object_name + "::U_problem_coefs"),
      d_default_U_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_U_bc_coef", Pointer<Database>(nullptr))),
      d_U_bc_coefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_U_bc_coef)),
      d_default_P_bc_coef(
          new LocationIndexRobinBcCoefs<NDIM>(d_object_name + "::default_P_bc_coef", Pointer<Database>(nullptr))),
      d_P_bc_coef(d_default_P_bc_coef)
{
    // Setup a default boundary condition object that specifies homogeneous
    // Dirichlet boundary conditions for the velocity and homogeneous Neumann
    // boundary conditions for the pressure.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        auto p_default_U_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_U_bc_coef);
        p_default_U_bc_coef->setBoundaryValue(2 * d, 0.0);
        p_default_U_bc_coef->setBoundaryValue(2 * d + 1, 0.0);
        auto p_default_P_bc_coef = dynamic_cast<LocationIndexRobinBcCoefs<NDIM>*>(d_default_P_bc_coef);
        p_default_P_bc_coef->setBoundarySlope(2 * d, 0.0);
        p_default_P_bc_coef->setBoundarySlope(2 * d + 1, 0.0);
    }

    // Initialize the boundary conditions objects.
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, d_default_U_bc_coef), d_default_P_bc_coef);

    if (input_db)
    {
        d_refine_type = input_db->getStringWithDefault("refine_type", d_refine_type);
        d_coarsen_type = input_db->getStringWithDefault("coarsen_type", d_coarsen_type);
        d_bdry_extrap_type = input_db->getStringWithDefault("bdry_extrap_type", d_bdry_extrap_type);
        d_bdry_interp_type = input_db->getStringWithDefault("bdry_interp_type", d_bdry_interp_type);
        d_use_cf_interpolation = input_db->getBoolWithDefault("use_cf_interpolation", d_use_cf_interpolation);
        d_consistent_type_2_bdry = input_db->getBoolWithDefault("consistent_type_2_bdry", d_consistent_type_2_bdry);
    }

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
    d_default_U_bc_coef = nullptr;
    delete d_default_P_bc_coef;
    d_default_P_bc_coef = nullptr;
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
    if (d_bc_helper) d_bc_helper->copyDataAtDirichletBoundaries(A_U_idx, U_scratch_idx);

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

    // Allocate scratch data.
    d_x->allocateVectorData();

    // Setup the interpolation transaction information.
    d_U_fill_pattern = new SideNoCornersFillPattern(SIDEG, false, false, true);
    d_P_fill_pattern = new CellNoCornersFillPattern(CELLG, false, false, true);
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_transaction_comps.resize(2);
    d_transaction_comps[0] = InterpolationTransactionComponent(d_x->getComponentDescriptorIndex(0),
                                                               in.getComponentDescriptorIndex(0),
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

    // Deallocate scratch data.
    d_x->deallocateVectorData();

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
        if (d_bc_helper)
        {
            const int U_idx = x->getComponentDescriptorIndex(0);
            const int P_idx = x->getComponentDescriptorIndex(1);
            StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
                d_U_bc_coefs, d_P_bc_coef, U_idx, P_idx, d_homogeneous_bc);
            d_bc_helper->enforceNormalVelocityBoundaryConditions(
                U_idx, P_idx, d_U_bc_coefs, d_new_time, d_homogeneous_bc);
        }
        StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
        apply(*x, *b);
        y.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&y, false), b);
        x->freeVectorComponents();
        b->freeVectorComponents();
    }
    const bool homogeneous_bc = true;
    if (d_bc_helper)
    {
        const int U_idx = y.getComponentDescriptorIndex(0);
        const int P_idx = y.getComponentDescriptorIndex(1);
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_U_bc_coefs, d_P_bc_coef, U_idx, P_idx, homogeneous_bc);
        d_bc_helper->enforceNormalVelocityBoundaryConditions(U_idx, P_idx, d_U_bc_coefs, d_new_time, homogeneous_bc);
        StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
    }
    return;
} // modifyRhsForBcs

void
StaggeredStokesOperator::imposeSolBcs(SAMRAIVectorReal<NDIM, double>& u)
{
    if (d_bc_helper)
    {
        const int U_idx = u.getComponentDescriptorIndex(0);
        const int P_idx = u.getComponentDescriptorIndex(1);
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_U_bc_coefs, d_P_bc_coef, U_idx, P_idx, d_homogeneous_bc);
        d_bc_helper->enforceNormalVelocityBoundaryConditions(U_idx, P_idx, d_U_bc_coefs, d_new_time, d_homogeneous_bc);
        StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_U_bc_coefs, d_P_bc_coef);
    }
    return;
} // imposeSolBcs

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
