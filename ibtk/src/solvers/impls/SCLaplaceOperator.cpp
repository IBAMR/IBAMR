// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
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

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/LaplaceOperator.h"
#include "ibtk/SCLaplaceOperator.h"
#include "ibtk/SideNoCornersFillPattern.h"
#include "ibtk/StaggeredPhysicalBoundaryHelper.h"
#include "ibtk/ibtk_utilities.h"

#include "Box.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "VariableFillPattern.h"
#include "tbox/Timer.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int SIDEG = 1;

// Timers.
static Timer* t_apply;
static Timer* t_initialize_operator_state;
static Timer* t_deallocate_operator_state;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

SCLaplaceOperator::SCLaplaceOperator(std::string object_name, Pointer<Database> input_db, const bool homogeneous_bc)
    : LaplaceOperator(std::move(object_name), homogeneous_bc)
{
    if (input_db)
    {
        if (input_db->isString("data_refine_type")) d_data_refine_type = input_db->getString("data_refine_type");
        if (input_db->isBool("use_cf_interpolation"))
            d_use_cf_interpolation = input_db->getBool("use_cf_interpolation");
        if (input_db->isString("data_coarsen_type")) d_data_coarsen_type = input_db->getString("data_coarsen_type");
        if (input_db->isString("bdry_extrap_type")) d_bdry_extrap_type = input_db->getString("bdry_extrap_type");
        if (input_db->isBool("use_consistent_type_2_bdry"))
            d_use_consistent_type_2_bdry = input_db->getBool("use_consistent_type_2_bdry");
    }

    // Setup the operator to use default vector-valued boundary conditions.
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr));

    // Setup Timers.
    IBTK_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBTK::SCLaplaceOperator::apply()");
                 t_initialize_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::SCLaplaceOperator::initializeOperatorState()");
                 t_deallocate_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::SCLaplaceOperator::deallocateOperatorState()"););
    return;
} // SCLaplaceOperator()

SCLaplaceOperator::~SCLaplaceOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~SCLaplaceOperator()

void
SCLaplaceOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBTK_TIMER_START(t_apply);

#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
    TBOX_ASSERT(d_bc_coefs.size() == NDIM);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        Pointer<SideVariable<NDIM, double> > x_sc_var = x.getComponentVariable(comp);
        Pointer<SideVariable<NDIM, double> > y_sc_var = y.getComponentVariable(comp);
        if (!x_sc_var || !y_sc_var)
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  encountered non-side centered vector components" << std::endl);
        }
        Pointer<SideDataFactory<NDIM, double> > x_factory = x_sc_var->getPatchDataFactory();
        Pointer<SideDataFactory<NDIM, double> > y_factory = y_sc_var->getPatchDataFactory();
        TBOX_ASSERT(x_factory);
        TBOX_ASSERT(y_factory);
        const unsigned int x_depth = x_factory->getDefaultDepth();
        const unsigned int y_depth = y_factory->getDefaultDepth();
        TBOX_ASSERT(x_depth == y_depth);
        if (x_depth != 1 || y_depth != 1)
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  each vector component must have data depth == 1" << std::endl);
        }
    }
#endif

    // Simultaneously fill ghost cell values for all components.
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> transaction_comps;
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        InterpolationTransactionComponent x_component(x.getComponentDescriptorIndex(comp),
                                                      d_data_refine_type,
                                                      d_use_cf_interpolation,
                                                      d_data_coarsen_type,
                                                      d_bdry_extrap_type,
                                                      d_use_consistent_type_2_bdry,
                                                      d_bc_coefs,
                                                      d_fill_pattern);
        transaction_comps.push_back(x_component);
    }
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the action of the operator.
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        Pointer<SideVariable<NDIM, double> > x_sc_var = x.getComponentVariable(comp);
        Pointer<SideVariable<NDIM, double> > y_sc_var = y.getComponentVariable(comp);
        const int x_idx = x.getComponentDescriptorIndex(comp);
        const int y_idx = y.getComponentDescriptorIndex(comp);
        d_hier_math_ops->laplace(y_idx, y_sc_var, d_poisson_spec, x_idx, x_sc_var, d_no_fill, 0.0);
        d_bc_helpers[comp]->copyDataAtDirichletBoundaries(y_idx, x_idx);
    }

    IBTK_TIMER_STOP(t_apply);
    return;
} // apply

void
SCLaplaceOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                           const SAMRAIVectorReal<NDIM, double>& out)
{
    IBTK_TIMER_START(t_initialize_operator_state);

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Setup operator state.
    d_hierarchy = in.getPatchHierarchy();
    d_coarsest_ln = in.getCoarsestLevelNumber();
    d_finest_ln = in.getFinestLevelNumber();

    d_ncomp = in.getNumberOfComponents();

#if !defined(NDEBUG)
    TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
    TBOX_ASSERT(d_ncomp == out.getNumberOfComponents());
#endif

    if (!d_hier_math_ops_external)
    {
        d_hier_math_ops =
            new HierarchyMathOps(d_object_name + "::HierarchyMathOps", d_hierarchy, d_coarsest_ln, d_finest_ln);
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_hier_math_ops);
#endif
    }

    // Setup cached BC data.
    d_bc_helpers.resize(d_ncomp);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        d_bc_helpers[comp] = new StaggeredPhysicalBoundaryHelper();
        d_bc_helpers[comp]->cacheBcCoefData(d_bc_coefs, d_solution_time, d_hierarchy);
    }

    // Setup the interpolation transaction information.
    d_fill_pattern = nullptr;
    if (d_poisson_spec.dIsConstant())
    {
        d_fill_pattern = new SideNoCornersFillPattern(SIDEG, /*overwrite_interior*/ false);
    }
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    d_transaction_comps.clear();
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        InterpolationTransactionComponent component(in.getComponentDescriptorIndex(comp),
                                                    d_data_refine_type,
                                                    d_use_cf_interpolation,
                                                    d_data_coarsen_type,
                                                    d_bdry_extrap_type,
                                                    d_use_consistent_type_2_bdry,
                                                    d_bc_coefs,
                                                    d_fill_pattern);
        d_transaction_comps.push_back(component);
    }

    // Initialize the interpolation operators.
    d_hier_bdry_fill = new HierarchyGhostCellInterpolation();
    d_hier_bdry_fill->initializeOperatorState(d_transaction_comps, d_hierarchy, d_coarsest_ln, d_finest_ln);

    // Indicate the operator is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
SCLaplaceOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_operator_state);

    // Deallocate the interpolation operators.
    d_hier_bdry_fill->deallocateOperatorState();
    d_hier_bdry_fill.setNull();
    d_transaction_comps.clear();
    d_fill_pattern.setNull();

    // Deallocate hierarchy math operations object.
    if (!d_hier_math_ops_external) d_hier_math_ops.setNull();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
