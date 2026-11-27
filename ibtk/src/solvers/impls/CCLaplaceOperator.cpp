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

#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/ibtk_enums.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Timer.h>

#include <CellVariable.h>
#include <MultiblockDataTranslator.h>
#include <PatchHierarchy.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <VariableFillPattern.h>

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
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

CCLaplaceOperator::CCLaplaceOperator(std::string object_name, Pointer<Database> input_db, const bool homogeneous_bc)
    : LaplaceOperator(std::move(object_name), homogeneous_bc)
{
    if (input_db)
    {
        if (input_db->isString("data_refine_type")) d_data_refine_type = input_db->getString("data_refine_type");
        if (input_db->isBool("use_cf_interpolation"))
            d_use_cf_interpolation = input_db->getBool("use_cf_interpolation");
        if (input_db->isString("data_coarsen_type")) d_data_coarsen_type = input_db->getString("data_coarsen_type");
        if (input_db->isString("bdry_extrap_type")) d_bdry_extrap_type = input_db->getString("bdry_extrap_type");
        if (input_db->isString("bdry_interp_type"))
        {
            const std::string type_name = input_db->getString("bdry_interp_type");
            const BdryInterpType bdry_interp_type = string_to_enum<BdryInterpType>(type_name);
            if (bdry_interp_type == BdryInterpType::UNKNOWN_BDRY_INTERP_TYPE)
            {
                TBOX_ERROR(d_object_name << "::CCLaplaceOperator():\n"
                                         << "  unknown bdry_interp_type \"" << type_name << "\"\n"
                                         << "  valid values are \"LINEAR\" and \"QUADRATIC\"\n");
            }
            setBdryInterpType(bdry_interp_type);
        }
        if (input_db->getBoolWithDefault("use_consistent_type_2_bdry", false))
        {
            TBOX_ERROR(d_object_name << "::CCLaplaceOperator():\n"
                                     << "  the input key use_consistent_type_2_bdry was removed and has no effect\n"
                                     << "  set it to FALSE or remove it\n");
        }
    }

    // Setup the operator to use default scalar-valued boundary conditions.
    setPhysicalBcCoef(nullptr);

    // Setup Timers.
    IBTK_DO_ONCE(t_apply = TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator::apply()");
                 t_initialize_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator::initializeOperatorState()");
                 t_deallocate_operator_state =
                     TimerManager::getManager()->getTimer("IBTK::CCLaplaceOperator::deallocateOperatorState()"););
    return;
} // CCLaplaceOperator()

CCLaplaceOperator::~CCLaplaceOperator()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~CCLaplaceOperator()

bool
CCLaplaceOperator::supportsBdryInterpType(const BdryInterpType bdry_interp_type) const
{
    return bdry_interp_type == BdryInterpType::LINEAR || bdry_interp_type == BdryInterpType::QUADRATIC;
} // supportsBdryInterpType

void
CCLaplaceOperator::setBdryExtrapType(const std::string& bdry_extrap_type)
{
    d_bdry_extrap_type = bdry_extrap_type;
    return;
} // setBdryExtrapType

void
CCLaplaceOperator::apply(SAMRAIVectorReal<NDIM, double>& x, SAMRAIVectorReal<NDIM, double>& y)
{
    IBTK_TIMER_START(t_apply);

#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        Pointer<CellVariable<NDIM, double>> x_cc_var = x.getComponentVariable(comp);
        Pointer<CellVariable<NDIM, double>> y_cc_var = y.getComponentVariable(comp);
        if (!x_cc_var || !y_cc_var)
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  encountered non-cell centered vector components" << std::endl);
        }
        Pointer<CellDataFactory<NDIM, double>> x_factory = x_cc_var->getPatchDataFactory();
        Pointer<CellDataFactory<NDIM, double>> y_factory = y_cc_var->getPatchDataFactory();
        TBOX_ASSERT(x_factory);
        TBOX_ASSERT(y_factory);
        const unsigned int x_depth = x_factory->getDefaultDepth();
        const unsigned int y_depth = y_factory->getDefaultDepth();
        TBOX_ASSERT(x_depth == y_depth);
        if (x_depth != d_bc_coefs.size() || y_depth != d_bc_coefs.size())
        {
            TBOX_ERROR(d_object_name << "::apply()\n"
                                     << "  each vector component must have data depth == " << d_bc_coefs.size() << "\n"
                                     << "  since d_bc_coefs.size() == " << d_bc_coefs.size() << std::endl);
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
                                                      /*consistent_type_2_bdry*/ false,
                                                      d_bc_coefs,
                                                      d_fill_pattern,
                                                      enum_to_string(d_bdry_interp_type));
        transaction_comps.push_back(x_component);
    }
    d_hier_bdry_fill->resetTransactionComponents(transaction_comps);
    d_hier_bdry_fill->setHomogeneousBc(d_homogeneous_bc);
    d_hier_bdry_fill->fillData(d_solution_time);
    d_hier_bdry_fill->resetTransactionComponents(d_transaction_comps);

    // Compute the action of the operator.
    for (int comp = 0; comp < d_ncomp; ++comp)
    {
        Pointer<CellVariable<NDIM, double>> x_cc_var = x.getComponentVariable(comp);
        Pointer<CellVariable<NDIM, double>> y_cc_var = y.getComponentVariable(comp);
        const int x_idx = x.getComponentDescriptorIndex(comp);
        const int y_idx = y.getComponentDescriptorIndex(comp);
        for (unsigned int l = 0; l < d_bc_coefs.size(); ++l)
        {
            d_hier_math_ops->laplace(y_idx,
                                     y_cc_var,
                                     d_poisson_spec,
                                     x_idx,
                                     x_cc_var,
                                     d_no_fill,
                                     0.0,
                                     0.0,
                                     -1,
                                     Pointer<CellVariable<NDIM, double>>(nullptr),
                                     l,
                                     l);
        }
    }

    IBTK_TIMER_STOP(t_apply);
    return;
} // apply

void
CCLaplaceOperator::initializeOperatorState(const SAMRAIVectorReal<NDIM, double>& in,
                                           const SAMRAIVectorReal<NDIM, double>& out)
{
    IBTK_TIMER_START(t_initialize_operator_state);

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
#else
    NULL_USE(out);
#endif

    // Reuse the HierarchyMathOps object unless a level was replaced.
    if (!d_hier_math_ops_external)
    {
        const bool levels_changed = d_own_hier_math_ops_levels.changed(*d_hierarchy, d_coarsest_ln, d_finest_ln);
        if (!d_own_hier_math_ops || levels_changed)
        {
            d_own_hier_math_ops =
                new HierarchyMathOps(d_object_name + "::HierarchyMathOps", d_hierarchy, d_coarsest_ln, d_finest_ln);
        }
        d_hier_math_ops = d_own_hier_math_ops;
    }
    else
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_hier_math_ops);
#endif
    }

    // Setup the interpolation transaction information.
    d_fill_pattern = nullptr;
    if (d_poisson_spec.dIsConstant())
    {
        d_fill_pattern = new CellNoCornersFillPattern(CELLG, /*overwrite_interior*/ false);
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
                                                    /*consistent_type_2_bdry*/ false,
                                                    d_bc_coefs,
                                                    d_fill_pattern,
                                                    enum_to_string(d_bdry_interp_type));
        d_transaction_comps.push_back(component);
    }

    // Reuse the interpolation operators while the levels and components are unchanged.
    d_hier_bdry_fill_cache.clearIfIndicesDiffer(d_transaction_comps);
    d_hier_bdry_fill = Pointer<HierarchyGhostCellInterpolation>(
        &d_hier_bdry_fill_cache.get(d_transaction_comps, d_hierarchy, d_coarsest_ln, d_finest_ln, d_homogeneous_bc),
        /*managed*/ false);

    // Indicate the operator is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
CCLaplaceOperator::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_operator_state);

    // Deallocate the interpolation operators.
    d_hier_bdry_fill.setNull();
    d_hier_bdry_fill_cache.clear();
    d_transaction_comps.clear();
    d_fill_pattern.setNull();

    // Deallocate hierarchy math operations object.
    d_own_hier_math_ops.setNull();
    d_own_hier_math_ops_levels = detail::HierarchyLevelsIdentity();
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
