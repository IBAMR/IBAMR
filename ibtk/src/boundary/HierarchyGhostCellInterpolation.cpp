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

#include <ibtk/CartCellDoubleCubicCoarsen.h>
#include <ibtk/CartCellDoubleQuadraticCFInterpolation.h>
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/CartSideDoubleQuadraticCFInterpolation.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/CoarseFineBoundaryRefinePatchStrategy.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/RefinePatchStrategySet.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellDataFactory.h>
#include <CellVariable.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenOperator.h>
#include <CoarsenPatchStrategy.h>
#include <CoarsenSchedule.h>
#include <EdgeVariable.h>
#include <NodeVariable.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchGeometry.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefinePatchStrategy.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideDataFactory.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableDatabase.h>
#include <VariableFillPattern.h>

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "../coarsen_ops/CubicCoarsenLayer.h"

#include <ibtk/namespaces.h> // IWYU pragma: keep

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
// Timers.
static Timer* t_initialize_operator_state;
static Timer* t_reset_transaction_component;
static Timer* t_reset_transaction_components;
static Timer* t_reinitialize_operator_state;
static Timer* t_deallocate_operator_state;
static Timer* t_fill_data;
static Timer* t_fill_data_coarsen;
static Timer* t_fill_data_refine;
static Timer* t_fill_data_set_physical_bcs;
} // namespace

// Extends the temporary coarse data quadratically beyond physical boundaries before the coarse-fine interpolation reads
// it.
class HierarchyGhostCellInterpolation::CoarseDataExtension : public RefinePatchStrategy<NDIM>
{
public:
    explicit CoarseDataExtension(const HierarchyGhostCellInterpolation& ghost_fill);

    void setPhysicalBoundaryConditions(Patch<NDIM>& patch,
                                       double fill_time,
                                       const IntVector<NDIM>& ghost_width_to_fill) override;

    IntVector<NDIM> getRefineOpStencilWidth() const override;

    void preprocessRefine(Patch<NDIM>& fine,
                          const Patch<NDIM>& coarse,
                          const Box<NDIM>& fine_box,
                          const IntVector<NDIM>& ratio) override;

    void postprocessRefine(Patch<NDIM>& fine,
                           const Patch<NDIM>& coarse,
                           const Box<NDIM>& fine_box,
                           const IntVector<NDIM>& ratio) override;

private:
    const HierarchyGhostCellInterpolation& d_ghost_fill;
    CartExtrapPhysBdryOp d_extrap_op;
};

HierarchyGhostCellInterpolation::CoarseDataExtension::CoarseDataExtension(
    const HierarchyGhostCellInterpolation& ghost_fill)
    : d_ghost_fill(ghost_fill)
{
    d_extrap_op.setExtrapolationType("QUADRATIC");
    return;
} // CoarseDataExtension

void
HierarchyGhostCellInterpolation::CoarseDataExtension::setPhysicalBoundaryConditions(
    Patch<NDIM>& /*patch*/,
    const double /*fill_time*/,
    const IntVector<NDIM>& /*ghost_width_to_fill*/)
{
    // intentionally blank
    return;
} // setPhysicalBoundaryConditions

IntVector<NDIM>
HierarchyGhostCellInterpolation::CoarseDataExtension::getRefineOpStencilWidth() const
{
    return 0;
} // getRefineOpStencilWidth

void
HierarchyGhostCellInterpolation::CoarseDataExtension::preprocessRefine(Patch<NDIM>& /*fine*/,
                                                                       const Patch<NDIM>& /*coarse*/,
                                                                       const Box<NDIM>& /*fine_box*/,
                                                                       const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
HierarchyGhostCellInterpolation::CoarseDataExtension::postprocessRefine(Patch<NDIM>& fine,
                                                                        const Patch<NDIM>& coarse,
                                                                        const Box<NDIM>& /*fine_box*/,
                                                                        const IntVector<NDIM>& /*ratio*/)
{
    // Only the patches of the hierarchy interpolate from coarse data at coarse-fine interfaces.
    if (!fine.inHierarchy())
    {
        return;
    }

    // A coarse patch narrower than three cells normal to a physical boundary keeps the extension set by the
    // extrapolation type of the transaction component.
    Pointer<PatchGeometry<NDIM>> pgeom = coarse.getPatchGeometry();
    if (!pgeom->getTouchesRegularBoundary())
    {
        return;
    }
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if ((pgeom->getTouchesRegularBoundary(d, 0) || pgeom->getTouchesRegularBoundary(d, 1)) &&
            coarse.getBox().numberCells(d) < 3)
        {
            return;
        }
    }

    // SAMRAI passes the temporary coarse patch as const; its data are scratch data of the schedule.
    Patch<NDIM>& coarse_patch = const_cast<Patch<NDIM>&>(coarse);
    for (unsigned int comp_idx = 0; comp_idx < d_ghost_fill.d_transaction_comps.size(); ++comp_idx)
    {
        if (d_ghost_fill.d_cf_bdry_ops[comp_idx])
        {
            const int data_idx = d_ghost_fill.d_transaction_comps[comp_idx].d_dst_data_idx;
            d_extrap_op.setPatchDataIndex(data_idx);
            d_extrap_op.setPhysicalBoundaryConditions(
                coarse_patch, 0.0, coarse.getPatchData(data_idx)->getGhostCellWidth());
        }
    }
    return;
} // postprocessRefine

/////////////////////////////// PUBLIC ///////////////////////////////////////

HierarchyGhostCellInterpolation::HierarchyGhostCellInterpolation()
{
    // Setup Timers.
    IBTK_DO_ONCE(
        t_initialize_operator_state =
            TimerManager::getManager()->getTimer("IBTK::HierarchyGhostCellInterpolation::initializeOperatorState()");
        t_reset_transaction_component =
            TimerManager::getManager()->getTimer("IBTK::HierarchyGhostCellInterpolation::resetTransactionComponent()");
        t_reset_transaction_components =
            TimerManager::getManager()->getTimer("IBTK::HierarchyGhostCellInterpolation::resetTransactionComponents()");
        t_reinitialize_operator_state =
            TimerManager::getManager()->getTimer("IBTK::HierarchyGhostCellInterpolation::reinitializeOperatorState()");
        t_deallocate_operator_state =
            TimerManager::getManager()->getTimer("IBTK::HierarchyGhostCellInterpolation::deallocateOperatorState()");
        t_fill_data = TimerManager::getManager()->getTimer("IBTK::HierarchyGhostCellInterpolation::fillData()");
        t_fill_data_coarsen =
            TimerManager::getManager()->getTimer("IBTK::HierarchyGhostCellInterpolation::fillData()[coarsen]");
        t_fill_data_refine =
            TimerManager::getManager()->getTimer("IBTK::HierarchyGhostCellInterpolation::fillData()[refine]");
        t_fill_data_set_physical_bcs = TimerManager::getManager()->getTimer(
            "IBTK::HierarchyGhostCellInterpolation::fillData()[set_physical_bcs]"););
    return;
} // HierarchyGhostCellInterpolation

HierarchyGhostCellInterpolation::~HierarchyGhostCellInterpolation()
{
    if (d_is_initialized) deallocateOperatorState();
    return;
} // ~HierarchyGhostCellInterpolation

void
HierarchyGhostCellInterpolation::setHomogeneousBc(const bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
} // setHomogeneousBc

void
HierarchyGhostCellInterpolation::initializeOperatorState(const InterpolationTransactionComponent transaction_comp,
                                                         const Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                         const int coarsest_ln,
                                                         const int finest_ln)
{
    IBTK_TIMER_START(t_initialize_operator_state);

    initializeOperatorState(
        std::vector<InterpolationTransactionComponent>(1, transaction_comp), hierarchy, coarsest_ln, finest_ln);

    IBTK_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
HierarchyGhostCellInterpolation::initializeOperatorState(
    const std::vector<InterpolationTransactionComponent>& transaction_comps,
    const Pointer<PatchHierarchy<NDIM>> hierarchy,
    const int coarsest_ln,
    const int finest_ln)
{
    IBTK_TIMER_START(t_initialize_operator_state);

    // Deallocate the operator state if the operator is already initialized.
    if (d_is_initialized) deallocateOperatorState();

    // Reset the transaction components.
    d_transaction_comps = transaction_comps;

    // Cache hierarchy data.
    d_hierarchy = hierarchy;
    d_grid_geom = d_hierarchy->getGridGeometry();
    d_coarsest_ln = coarsest_ln == invalid_level_number ? 0 : coarsest_ln;
    d_finest_ln = finest_ln == invalid_level_number ? d_hierarchy->getFinestLevelNumber() : finest_ln;

    // Register the cubic coarsen operators with the grid geometry object. They serve the objects that do not mark
    // the invalid ghost values of their data.
    IBTK_DO_ONCE(d_grid_geom->addSpatialCoarsenOperator(new CartCellDoubleCubicCoarsen());
                 d_grid_geom->addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen()););

    // Create the cubic coarsen operators of this object.
    d_cc_cubic_coarsen_op = new CartCellDoubleCubicCoarsen(true);
    d_sc_cubic_coarsen_op = new CartSideDoubleCubicCoarsen(true);

    // Setup cached coarsen algorithms and schedules.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    bool registered_coarsen_op = false;
    d_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    for (const auto& transaction_comp : d_transaction_comps)
    {
        const std::string& coarsen_op_name = transaction_comp.d_coarsen_op_name;
        if (coarsen_op_name != "NONE")
        {
            const int src_data_idx = transaction_comp.d_src_data_idx;
            Pointer<Variable<NDIM>> var;
            var_db->mapIndexToVariable(src_data_idx, var);
#if !defined(NDEBUG)
            TBOX_ASSERT(var);
#endif
            Pointer<CoarsenOperator<NDIM>> coarsen_op = getCoarsenOperator(var, transaction_comp);
#if !defined(NDEBUG)
            TBOX_ASSERT(coarsen_op);
#endif
            d_coarsen_alg->registerCoarsen(src_data_idx, src_data_idx, coarsen_op);
            registered_coarsen_op = true;
        }
    }

    d_coarsen_strategy = nullptr;

    setUpCubicScratch();
    d_coarsen_scheds.resize(d_finest_ln + 1);
    d_scratch_fill_scheds.resize(d_finest_ln + 1);
    if (registered_coarsen_op)
    {
        for (int src_ln = std::max(1, d_coarsest_ln); src_ln <= d_finest_ln; ++src_ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(src_ln);
            Pointer<PatchLevel<NDIM>> coarser_level = d_hierarchy->getPatchLevel(src_ln - 1);
            if (usesCubicScratch(src_ln))
            {
                d_coarsen_scheds[src_ln] =
                    d_scratch_coarsen_alg->createSchedule(coarser_level, level, d_coarsen_strategy.get());
                d_scratch_fill_scheds[src_ln] = d_scratch_fill_alg->createSchedule(level);
            }
            else
            {
                d_coarsen_scheds[src_ln] =
                    d_coarsen_alg->createSchedule(coarser_level, level, d_coarsen_strategy.get());
            }
        }
    }

    // Setup cached refine algorithms and schedules.
    d_cf_bdry_ops.resize(d_transaction_comps.size());
    d_extrap_bc_ops.resize(d_transaction_comps.size());
    d_cc_robin_bc_ops.resize(d_transaction_comps.size());
    d_sc_robin_bc_ops.resize(d_transaction_comps.size());
    d_refine_alg = new RefineAlgorithm<NDIM>();
    std::vector<RefinePatchStrategy<NDIM>*> refine_patch_strategies;
    for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const int dst_data_idx = d_transaction_comps[comp_idx].d_dst_data_idx;
        const int src_data_idx = d_transaction_comps[comp_idx].d_src_data_idx;
        const std::string& phys_bdry_type = d_transaction_comps[comp_idx].d_phys_bdry_type;
        Pointer<Variable<NDIM>> var;
        var_db->mapIndexToVariable(src_data_idx, var);
        Pointer<CellVariable<NDIM, double>> cc_var = var;
        Pointer<NodeVariable<NDIM, double>> nc_var = var;
        Pointer<SideVariable<NDIM, double>> sc_var = var;
        Pointer<EdgeVariable<NDIM, double>> ec_var = var;
        Pointer<FaceVariable<NDIM, double>> fc_var = var;
        Pointer<RefineOperator<NDIM>> refine_op = nullptr;
        d_cf_bdry_ops[comp_idx] = nullptr;
        Pointer<VariableFillPattern<NDIM>> fill_pattern = d_transaction_comps[comp_idx].d_fill_pattern;
        if (cc_var)
        {
            if (d_transaction_comps[comp_idx].d_refine_op_name != "NONE")
            {
                refine_op = d_grid_geom->lookupRefineOperator(var, d_transaction_comps[comp_idx].d_refine_op_name);
            }
            if (d_transaction_comps[comp_idx].d_use_cf_bdry_interpolation)
            {
                d_cf_bdry_ops[comp_idx] = new CartCellDoubleQuadraticCFInterpolation();
                d_cf_bdry_ops[comp_idx]->setConsistentInterpolationScheme(
                    d_transaction_comps[comp_idx].d_consistent_type_2_bdry);
                d_cf_bdry_ops[comp_idx]->setPatchDataIndex(dst_data_idx);
                d_cf_bdry_ops[comp_idx]->setPatchHierarchy(d_hierarchy);
                refine_patch_strategies.push_back(d_cf_bdry_ops[comp_idx]);
            }
        }
        else if (nc_var)
        {
            if (d_transaction_comps[comp_idx].d_refine_op_name != "NONE")
            {
                refine_op = d_grid_geom->lookupRefineOperator(var, d_transaction_comps[comp_idx].d_refine_op_name);
            }
            if (d_transaction_comps[comp_idx].d_use_cf_bdry_interpolation)
            {
                TBOX_ERROR("not supported yet.\n");
            }
        }
        else if (sc_var)
        {
            if (d_transaction_comps[comp_idx].d_refine_op_name != "NONE")
            {
                refine_op = d_grid_geom->lookupRefineOperator(var, d_transaction_comps[comp_idx].d_refine_op_name);
            }
            if (d_transaction_comps[comp_idx].d_use_cf_bdry_interpolation)
            {
                d_cf_bdry_ops[comp_idx] = new CartSideDoubleQuadraticCFInterpolation();
                d_cf_bdry_ops[comp_idx]->setConsistentInterpolationScheme(
                    d_transaction_comps[comp_idx].d_consistent_type_2_bdry);
                d_cf_bdry_ops[comp_idx]->setPatchDataIndex(dst_data_idx);
                d_cf_bdry_ops[comp_idx]->setPatchHierarchy(d_hierarchy);
                refine_patch_strategies.push_back(d_cf_bdry_ops[comp_idx]);
            }
        }
        else if (ec_var)
        {
            if (d_transaction_comps[comp_idx].d_refine_op_name != "NONE")
            {
                refine_op = d_grid_geom->lookupRefineOperator(var, d_transaction_comps[comp_idx].d_refine_op_name);
            }
            if (d_transaction_comps[comp_idx].d_use_cf_bdry_interpolation)
            {
                TBOX_ERROR("not supported yet.\n");
            }
        }
        else if (fc_var)
        {
            if (d_transaction_comps[comp_idx].d_refine_op_name != "NONE")
            {
                refine_op = d_grid_geom->lookupRefineOperator(var, d_transaction_comps[comp_idx].d_refine_op_name);
            }
            if (d_transaction_comps[comp_idx].d_use_cf_bdry_interpolation)
            {
                TBOX_ERROR("not supported yet.\n");
            }
        }
        else
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::initializeOperatorState():\n"
                       << "  only double-precision cell-, node-, side-, or edge-centered data is "
                          "presently supported."
                       << std::endl);
        }

        d_refine_alg->registerRefine(dst_data_idx, src_data_idx, dst_data_idx, refine_op, fill_pattern);

        const std::string& phys_bdry_extrap_type = d_transaction_comps[comp_idx].d_phys_bdry_extrap_type;
        if (phys_bdry_extrap_type != "NONE")
        {
            d_extrap_bc_ops[comp_idx] = new CartExtrapPhysBdryOp(dst_data_idx, phys_bdry_extrap_type);
            refine_patch_strategies.push_back(d_extrap_bc_ops[comp_idx]);
        }

        const std::vector<RobinBcCoefStrategy<NDIM>*>& robin_bc_coefs = d_transaction_comps[comp_idx].d_robin_bc_coefs;
        bool null_bc_coefs = true;
        for (const auto& robin_bc_coef : robin_bc_coefs)
        {
            if (robin_bc_coef) null_bc_coefs = false;
        }
        if (!null_bc_coefs && cc_var)
        {
            d_cc_robin_bc_ops[comp_idx] =
                new CartCellRobinPhysBdryOp(dst_data_idx, robin_bc_coefs, d_homogeneous_bc, phys_bdry_type);
        }
        if (!null_bc_coefs && sc_var)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(robin_bc_coefs.size() % NDIM == 0);
#endif
            d_sc_robin_bc_ops[comp_idx] =
                new CartSideRobinPhysBdryOp(dst_data_idx, robin_bc_coefs, d_homogeneous_bc, phys_bdry_type);
        }
    }

    // The coarse data is extended before the coarse-fine interpolations read it.
    d_coarse_extension = nullptr;
    if (std::any_of(
            d_cf_bdry_ops.begin(), d_cf_bdry_ops.end(), [](const auto& cf_bdry_op) { return !cf_bdry_op.isNull(); }))
    {
        d_coarse_extension = std::make_unique<CoarseDataExtension>(*this);
        refine_patch_strategies.insert(refine_patch_strategies.begin(), d_coarse_extension.get());
    }

    d_refine_strategy =
        std::make_unique<RefinePatchStrategySet>(refine_patch_strategies.begin(), refine_patch_strategies.end(), false);

    d_refine_scheds.resize(d_finest_ln + 1);
    for (int dst_ln = d_coarsest_ln; dst_ln <= d_finest_ln; ++dst_ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(dst_ln);
        d_refine_scheds[dst_ln] = d_refine_alg->createSchedule(level, dst_ln - 1, d_hierarchy, d_refine_strategy.get());
    }

    // Setup physical BC type.
    setHomogeneousBc(d_homogeneous_bc);

    // Indicate the operator is initialized.
    d_is_initialized = true;

    IBTK_TIMER_STOP(t_initialize_operator_state);
    return;
} // initializeOperatorState

void
HierarchyGhostCellInterpolation::resetTransactionComponent(const InterpolationTransactionComponent& transaction_comp)
{
    IBTK_TIMER_START(t_reset_transaction_component);

#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif
    if (d_transaction_comps.size() != 1)
    {
        TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponent():"
                   << "  invalid reset operation.  attempting to change the number of registered "
                      "interpolation transaction components.\n");
    }
    resetTransactionComponents(std::vector<InterpolationTransactionComponent>(1, transaction_comp));

    IBTK_TIMER_STOP(t_reset_transaction_component);
    return;
} // resetTransactionComponent

void
HierarchyGhostCellInterpolation::resetTransactionComponents(
    const std::vector<InterpolationTransactionComponent>& transaction_comps)
{
    IBTK_TIMER_START(t_reset_transaction_components);

#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif
    if (d_transaction_comps.size() != transaction_comps.size())
    {
        TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():"
                   << "  invalid reset operation.  attempting to change the number of registered "
                      "interpolation transaction components.\n");
    }

    // Reset the transaction components.
    d_transaction_comps = transaction_comps;

    // Reset cached coarsen algorithms and schedules.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    bool registered_coarsen_op = false;
    d_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    for (const auto& transaction_comp : d_transaction_comps)
    {
        const std::string& coarsen_op_name = transaction_comp.d_coarsen_op_name;
        if (coarsen_op_name != "NONE")
        {
            const int src_data_idx = transaction_comp.d_src_data_idx;
            Pointer<Variable<NDIM>> var;
            var_db->mapIndexToVariable(src_data_idx, var);
#if !defined(NDEBUG)
            TBOX_ASSERT(var);
#endif
            Pointer<CoarsenOperator<NDIM>> coarsen_op = getCoarsenOperator(var, transaction_comp);
#if !defined(NDEBUG)
            TBOX_ASSERT(coarsen_op);
#endif
            d_coarsen_alg->registerCoarsen(src_data_idx, src_data_idx, coarsen_op);
            registered_coarsen_op = true;
        }
    }

    setUpCubicScratch();
    if (registered_coarsen_op)
    {
        for (int src_ln = std::max(1, d_coarsest_ln); src_ln <= d_finest_ln; ++src_ln)
        {
            if (usesCubicScratch(src_ln) != static_cast<bool>(d_scratch_fill_scheds[src_ln]))
            {
                TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():"
                           << "  invalid reset operation.  attempting to change the cubic coarsening of a level "
                              "with a refinement ratio less than 4.\n");
            }
            if (usesCubicScratch(src_ln))
            {
                d_scratch_coarsen_alg->resetSchedule(d_coarsen_scheds[src_ln]);
                d_scratch_fill_alg->resetSchedule(d_scratch_fill_scheds[src_ln]);
            }
            else
            {
                d_coarsen_alg->resetSchedule(d_coarsen_scheds[src_ln]);
            }
        }
    }

    // Reset cached refine algorithms and schedules.
    d_refine_alg = new RefineAlgorithm<NDIM>();
    for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const int dst_data_idx = d_transaction_comps[comp_idx].d_dst_data_idx;
        const int src_data_idx = d_transaction_comps[comp_idx].d_src_data_idx;
        Pointer<Variable<NDIM>> var;
        var_db->mapIndexToVariable(src_data_idx, var);
        Pointer<CellVariable<NDIM, double>> cc_var = var;
        Pointer<NodeVariable<NDIM, double>> nc_var = var;
        Pointer<SideVariable<NDIM, double>> sc_var = var;
        Pointer<EdgeVariable<NDIM, double>> ec_var = var;
        Pointer<RefineOperator<NDIM>> refine_op = nullptr;
        Pointer<VariableFillPattern<NDIM>> fill_pattern = d_transaction_comps[comp_idx].d_fill_pattern;
        if (d_cf_bdry_ops[comp_idx]) d_cf_bdry_ops[comp_idx]->setPatchDataIndex(dst_data_idx);
        if (cc_var)
        {
            if (d_transaction_comps[comp_idx].d_refine_op_name != "NONE")
            {
                refine_op = d_grid_geom->lookupRefineOperator(var, d_transaction_comps[comp_idx].d_refine_op_name);
            }
        }
        else if (nc_var)
        {
            if (d_transaction_comps[comp_idx].d_refine_op_name != "NONE")
            {
                refine_op = d_grid_geom->lookupRefineOperator(var, d_transaction_comps[comp_idx].d_refine_op_name);
            }
        }
        else if (sc_var)
        {
            if (d_transaction_comps[comp_idx].d_refine_op_name != "NONE")
            {
                refine_op = d_grid_geom->lookupRefineOperator(var, d_transaction_comps[comp_idx].d_refine_op_name);
            }
        }
        else if (ec_var)
        {
            if (d_transaction_comps[comp_idx].d_refine_op_name != "NONE")
            {
                refine_op = d_grid_geom->lookupRefineOperator(var, d_transaction_comps[comp_idx].d_refine_op_name);
            }
        }
        else
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                       << "  only double-precision cell-, node-, side-, or edge-centered data is "
                          "presently supported."
                       << std::endl);
        }

        d_refine_alg->registerRefine(dst_data_idx, src_data_idx, dst_data_idx, refine_op, fill_pattern);

        const std::string& phys_bdry_extrap_type = d_transaction_comps[comp_idx].d_phys_bdry_extrap_type;
        if (d_extrap_bc_ops[comp_idx])
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(phys_bdry_extrap_type != "NONE");
#endif
            d_extrap_bc_ops[comp_idx]->setPatchDataIndex(dst_data_idx);
            d_extrap_bc_ops[comp_idx]->setExtrapolationType(phys_bdry_extrap_type);
        }
#if !defined(NDEBUG)
        else
        {
            TBOX_ASSERT(phys_bdry_extrap_type == "NONE");
        }

#endif
        const std::vector<RobinBcCoefStrategy<NDIM>*>& robin_bc_coefs = d_transaction_comps[comp_idx].d_robin_bc_coefs;
#if !defined(NDEBUG)
        bool null_bc_coefs = true;
        for (const auto& robin_bc_coef : robin_bc_coefs)
        {
            if (robin_bc_coef) null_bc_coefs = false;
        }
#endif
        if (d_cc_robin_bc_ops[comp_idx])
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(!null_bc_coefs);
            TBOX_ASSERT(cc_var);
#endif
            d_cc_robin_bc_ops[comp_idx]->setPhysicalBcCoefs(robin_bc_coefs);
            d_cc_robin_bc_ops[comp_idx]->setPatchDataIndex(dst_data_idx);
        }
        if (d_sc_robin_bc_ops[comp_idx])
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(!null_bc_coefs);
            TBOX_ASSERT(sc_var);
            TBOX_ASSERT(robin_bc_coefs.size() % NDIM == 0);
#endif
            d_sc_robin_bc_ops[comp_idx]->setPhysicalBcCoefs(robin_bc_coefs);
            d_sc_robin_bc_ops[comp_idx]->setPatchDataIndex(dst_data_idx);
        }
    }

    for (int dst_ln = d_coarsest_ln; dst_ln <= d_finest_ln; ++dst_ln)
    {
        d_refine_alg->resetSchedule(d_refine_scheds[dst_ln]);
    }

    IBTK_TIMER_STOP(t_reset_transaction_components);
    return;
} // resetTransactionComponents

void
HierarchyGhostCellInterpolation::reinitializeOperatorState(Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_reinitialize_operator_state);

    initializeOperatorState(d_transaction_comps, hierarchy);

    IBTK_TIMER_STOP(t_reinitialize_operator_state);
    return;
} // reinitializeOperatorState

void
HierarchyGhostCellInterpolation::deallocateOperatorState()
{
    if (!d_is_initialized) return;

    IBTK_TIMER_START(t_deallocate_operator_state);

    // Clear cached refinement operators.
    d_cf_bdry_ops.clear();
    d_extrap_bc_ops.clear();
    d_cc_robin_bc_ops.clear();
    d_sc_robin_bc_ops.clear();

    // Clear cached communication schedules.
    d_coarsen_alg.setNull();
    d_coarsen_strategy = nullptr;
    d_coarsen_scheds.clear();
    d_scratch_coarsen_alg.setNull();
    d_scratch_fill_alg.setNull();
    d_scratch_fill_scheds.clear();
    d_cubic_scratch_idxs.clear();
    d_cc_cubic_coarsen_op.setNull();
    d_sc_cubic_coarsen_op.setNull();

    d_refine_alg.setNull();
    d_refine_strategy = nullptr;
    d_coarse_extension = nullptr;
    d_refine_scheds.clear();

    // Indicate that the operator is NOT initialized.
    d_is_initialized = false;

    IBTK_TIMER_STOP(t_deallocate_operator_state);
    return;
} // deallocateOperatorState

void
HierarchyGhostCellInterpolation::fillData(double fill_time)
{
    IBTK_TIMER_START(t_fill_data);

#if !defined(NDEBUG)
    TBOX_ASSERT(d_is_initialized);
#endif
    // Ensure the boundary condition objects are in the correct state.
    for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        if (d_cc_robin_bc_ops[comp_idx]) d_cc_robin_bc_ops[comp_idx]->setHomogeneousBc(d_homogeneous_bc);
        if (d_sc_robin_bc_ops[comp_idx]) d_sc_robin_bc_ops[comp_idx]->setHomogeneousBc(d_homogeneous_bc);
    }

    // Synchronize data on the patch hierarchy prior to filling ghost cell
    // values.
    IBTK_TIMER_START(t_fill_data_coarsen);
    for (int src_ln = d_finest_ln; src_ln >= std::max(1, d_coarsest_ln); --src_ln)
    {
        if (d_scratch_fill_scheds[src_ln])
        {
            fillCubicScratch(src_ln, fill_time);
        }
        if (d_coarsen_scheds[src_ln]) d_coarsen_scheds[src_ln]->coarsenData();
        if (d_scratch_fill_scheds[src_ln])
        {
            deallocateCubicScratch(src_ln);
        }
    }
    IBTK_TIMER_STOP(t_fill_data_coarsen);

    // Perform the initial data fill, using extrapolation to determine ghost
    // cell values at physical boundaries.
    IBTK_TIMER_START(t_fill_data_refine);
    for (int dst_ln = d_coarsest_ln; dst_ln <= d_finest_ln; ++dst_ln)
    {
        if (d_refine_scheds[dst_ln]) d_refine_scheds[dst_ln]->fillData(fill_time);
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(dst_ln);
        const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
            {
                if (d_cf_bdry_ops[comp_idx])
                {
                    const int dst_data_idx = d_transaction_comps[comp_idx].d_dst_data_idx;
                    const IntVector<NDIM>& ghost_width_to_fill = patch->getPatchData(dst_data_idx)->getGhostCellWidth();
                    d_cf_bdry_ops[comp_idx]->computeNormalExtension(*patch, ratio, ghost_width_to_fill);
                }
            }
        }
    }
    IBTK_TIMER_STOP(t_fill_data_refine);

    // Set Robin boundary conditions at physical boundaries.
    IBTK_TIMER_START(t_fill_data_set_physical_bcs);
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            if (patch->getPatchGeometry()->getTouchesRegularBoundary())
            {
                for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
                {
                    if (d_cc_robin_bc_ops[comp_idx])
                    {
                        const int dst_data_idx = d_transaction_comps[comp_idx].d_dst_data_idx;
                        const IntVector<NDIM>& ghost_width_to_fill =
                            patch->getPatchData(dst_data_idx)->getGhostCellWidth();
                        d_cc_robin_bc_ops[comp_idx]->setPhysicalBoundaryConditions(
                            *patch, fill_time, ghost_width_to_fill);
                    }
                    if (d_sc_robin_bc_ops[comp_idx])
                    {
                        const int dst_data_idx = d_transaction_comps[comp_idx].d_dst_data_idx;
                        const IntVector<NDIM>& ghost_width_to_fill =
                            patch->getPatchData(dst_data_idx)->getGhostCellWidth();
                        d_sc_robin_bc_ops[comp_idx]->setPhysicalBoundaryConditions(
                            *patch, fill_time, ghost_width_to_fill);
                    }
                }
            }
        }
    }
    IBTK_TIMER_STOP(t_fill_data_set_physical_bcs);

    IBTK_TIMER_STOP(t_fill_data);
    return;
} // fillData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

Pointer<CoarsenOperator<NDIM>>
HierarchyGhostCellInterpolation::getCoarsenOperator(const Pointer<Variable<NDIM>>& var,
                                                    const InterpolationTransactionComponent& comp) const
{
    if (needsCubicScratch(comp))
    {
        Pointer<CellVariable<NDIM, double>> cc_var = var;
        Pointer<SideVariable<NDIM, double>> sc_var = var;
        if (cc_var)
        {
            return d_cc_cubic_coarsen_op;
        }
        if (sc_var)
        {
            return d_sc_cubic_coarsen_op;
        }
    }
    return d_grid_geom->lookupCoarsenOperator(var, comp.d_coarsen_op_name);
} // getCoarsenOperator

bool
HierarchyGhostCellInterpolation::needsCubicScratch(const InterpolationTransactionComponent& comp) const
{
    // An in-place fill with a refine operator overwrites the fine values that the ratio 2 stencil reads, so such
    // components keep the registered operator.
    return comp.d_coarsen_op_name == "CUBIC_COARSEN" && comp.d_refine_op_name == "NONE";
} // needsCubicScratch

void
HierarchyGhostCellInterpolation::setUpCubicScratch()
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("IBTK::HierarchyGhostCellInterpolation::CUBIC_COARSEN");
    d_scratch_coarsen_alg.setNull();
    d_scratch_fill_alg.setNull();
    d_cubic_scratch_idxs.assign(d_transaction_comps.size(), invalid_index);

    // Register one scratch data index with a one-cell ghost layer for each source of cubic coarsening.
    bool registered_scratch = false;
    for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const int src_data_idx = d_transaction_comps[comp_idx].d_src_data_idx;
        if (!needsCubicScratch(d_transaction_comps[comp_idx]))
        {
            continue;
        }
        Pointer<Variable<NDIM>> var;
        var_db->mapIndexToVariable(src_data_idx, var);
        Pointer<CellVariable<NDIM, double>> cc_var = var;
        Pointer<SideVariable<NDIM, double>> sc_var = var;
        const std::string name = "IBTK::HierarchyGhostCellInterpolation::cubic_scratch_" + var->getName();
        Pointer<PatchDataFactory<NDIM>> factory = var->getPatchDataFactory();
        Pointer<CellDataFactory<NDIM, double>> cc_factory = factory;
        Pointer<SideDataFactory<NDIM, double>> sc_factory = factory;
        const int depth = cc_var ? cc_factory->getDefaultDepth() : sc_factory->getDefaultDepth();
        Pointer<Variable<NDIM>> scratch_var = var_db->checkVariableExists(name) ? var_db->getVariable(name) : nullptr;
        if (!scratch_var && cc_var)
        {
            scratch_var = new CellVariable<NDIM, double>(name, depth);
        }
        else if (!scratch_var && sc_var)
        {
            scratch_var = new SideVariable<NDIM, double>(name, depth);
        }
        if (!scratch_var)
        {
            continue;
        }
        d_cubic_scratch_idxs[comp_idx] = var_db->registerVariableAndContext(scratch_var, ctx, IntVector<NDIM>(1));
        registered_scratch = true;
    }
    if (!registered_scratch)
    {
        return;
    }

    // Coarsen the scratch data in place of the source data. Fill the scratch data from the source data of the
    // same level.
    d_scratch_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    d_scratch_fill_alg = new RefineAlgorithm<NDIM>();
    for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const std::string& coarsen_op_name = d_transaction_comps[comp_idx].d_coarsen_op_name;
        if (coarsen_op_name == "NONE")
        {
            continue;
        }
        const int src_data_idx = d_transaction_comps[comp_idx].d_src_data_idx;
        const int scratch_idx = d_cubic_scratch_idxs[comp_idx];
        Pointer<Variable<NDIM>> var;
        var_db->mapIndexToVariable(src_data_idx, var);
        Pointer<CoarsenOperator<NDIM>> coarsen_op = getCoarsenOperator(var, d_transaction_comps[comp_idx]);
        if (scratch_idx == invalid_index)
        {
            d_scratch_coarsen_alg->registerCoarsen(src_data_idx, src_data_idx, coarsen_op);
        }
        else
        {
            d_scratch_coarsen_alg->registerCoarsen(src_data_idx, scratch_idx, coarsen_op);
            d_scratch_fill_alg->registerRefine(scratch_idx, src_data_idx, scratch_idx, nullptr);
        }
    }
    return;
} // setUpCubicScratch

bool
HierarchyGhostCellInterpolation::usesCubicScratch(const int ln) const
{
    if (!d_scratch_coarsen_alg)
    {
        return false;
    }
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
    return level->getRatioToCoarserLevel().min() < 4;
} // usesCubicScratch

void
HierarchyGhostCellInterpolation::fillCubicScratch(const int ln, const double fill_time)
{
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
    const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();

    // Mark the ghost layer of the scratch data as invalid. The schedule fills the patch interiors and the ghost values
    // that other patches of the level cover.
    for (const int scratch_idx : d_cubic_scratch_idxs)
    {
        if (scratch_idx == invalid_index)
        {
            continue;
        }
        if (!level->checkAllocated(scratch_idx))
        {
            level->allocatePatchData(scratch_idx, fill_time);
        }
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double>> cc_scratch = patch->getPatchData(scratch_idx);
            Pointer<SideData<NDIM, double>> sc_scratch = patch->getPatchData(scratch_idx);
            if (cc_scratch)
            {
                mark_cubic_coarsen_layer(cc_scratch->getArrayData(), patch_box);
            }
            else
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    mark_cubic_coarsen_layer(sc_scratch->getArrayData(axis),
                                             SideGeometry<NDIM>::toSideBox(patch_box, axis));
                }
            }
        }
    }
    d_scratch_fill_scheds[ln]->fillData(fill_time, false);

    // Replace the invalid values of the layer by the extrapolation from the patch.
    for (const int scratch_idx : d_cubic_scratch_idxs)
    {
        if (scratch_idx == invalid_index)
        {
            continue;
        }
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<CellData<NDIM, double>> cc_scratch = patch->getPatchData(scratch_idx);
            Pointer<SideData<NDIM, double>> sc_scratch = patch->getPatchData(scratch_idx);
            if (cc_scratch)
            {
                IntVector<NDIM> width(0);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    width(d) = ratio(d) < 4 ? 1 : 0;
                }
                extrapolate_cubic_coarsen_layer(cc_scratch->getArrayData(), patch_box, width);
            }
            else
            {
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    IntVector<NDIM> width(0);
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        width(d) = (d != axis && ratio(d) < 4) ? 1 : 0;
                    }
                    extrapolate_cubic_coarsen_layer(
                        sc_scratch->getArrayData(axis), SideGeometry<NDIM>::toSideBox(patch_box, axis), width);
                }
            }
        }
    }
    return;
} // fillCubicScratch

void
HierarchyGhostCellInterpolation::deallocateCubicScratch(const int ln)
{
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
    for (const int scratch_idx : d_cubic_scratch_idxs)
    {
        if (scratch_idx != invalid_index && level->checkAllocated(scratch_idx))
        {
            level->deallocatePatchData(scratch_idx);
        }
    }
    return;
} // deallocateCubicScratch

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
