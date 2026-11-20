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

#include <ibtk/CFOwnerFill.h>
#include <ibtk/CartCellDoubleCubicCoarsen.h>
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/RefinePatchStrategySet.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Pointer.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CellDataFactory.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenOperator.h>
#include <CoarsenPatchStrategy.h>
#include <CoarsenSchedule.h>
#include <EdgeDataFactory.h>
#include <EdgeVariable.h>
#include <FaceDataFactory.h>
#include <NodeDataFactory.h>
#include <NodeVariable.h>
#include <Patch.h>
#include <PatchData.h>
#include <PatchDataFactory.h>
#include <PatchDescriptor.h>
#include <PatchGeometry.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefinePatchStrategy.h>
#include <RefineSchedule.h>
#include <SideDataFactory.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableDatabase.h>
#include <VariableFillPattern.h>

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

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

// Works around IBSAMRAI2#67: in an in-place refine fill, the side values of a fine patch on the coarse-fine interface
// are overwritten with refined coarse values. Remove once IBAMR requires an IBSAMRAI2 release that fixes #67.
bool
preserves_interface_faces(const HierarchyGhostCellInterpolation::InterpolationTransactionComponent& comp)
{
    return comp.d_dst_data_idx == comp.d_src_data_idx && comp.d_refine_op_name != "NONE";
}

// The centering (0 to 4 for cell, node, side, edge, and face data, -1 otherwise) and depth of the data with index idx.
std::pair<int, int>
get_data_kind(const int idx)
{
    Pointer<PatchDataFactory<NDIM>> factory =
        VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor()->getPatchDataFactory(idx);
    if (Pointer<CellDataFactory<NDIM, double>> f = factory) return { 0, f->getDefaultDepth() };
    if (Pointer<NodeDataFactory<NDIM, double>> f = factory) return { 1, f->getDefaultDepth() };
    if (Pointer<SideDataFactory<NDIM, double>> f = factory) return { 2, f->getDefaultDepth() };
    if (Pointer<EdgeDataFactory<NDIM, double>> f = factory) return { 3, f->getDefaultDepth() };
    if (Pointer<FaceDataFactory<NDIM, double>> f = factory) return { 4, f->getDefaultDepth() };
    return { -1, 0 };
}

std::string
data_kind_to_string(const std::pair<int, int>& kind)
{
    static const std::array<std::string, 5> names = { "cell", "node", "side", "edge", "face" };
    return (kind.first < 0 ? std::string("other") : names[kind.first]) + "-centered data of depth " +
           std::to_string(kind.second);
}

// In an in-place refine fill, IBSAMRAI2#67 also overwrites the patch-border values of node, edge, and face data, and
// only side data are restored. Remove this check once IBAMR requires an IBSAMRAI2 release that fixes #67.
void
check_in_place_fill_of_patch_border_data(const HierarchyGhostCellInterpolation::InterpolationTransactionComponent& comp,
                                         const Pointer<Variable<NDIM>>& var,
                                         const int finest_ln,
                                         const std::string& function_name)
{
    if (finest_ln < 1 || !preserves_interface_faces(comp)) return;
    Pointer<SideVariable<NDIM, double>> sc_var = var;
    if (sc_var || !var->getPatchDataFactory()->dataLivesOnPatchBorder()) return;
    TBOX_ERROR("HierarchyGhostCellInterpolation::"
               << function_name << "():\n"
               << "  variable " << var->getName() << " is filled in place with refine type " << comp.d_refine_op_name
               << " on a level with a coarser level, which overwrites its patch-border values\n"
               << "  fill it from a separate source patch data index\n");
}

// The refine strategies of the fill, which also restore the faces along the coarse-fine interface of a patch right
// after the refine operator fills it, before the schedule copies the values to other patches and sets the physical
// boundary values.
class RefineStrategies : public RefinePatchStrategySet
{
public:
    template <typename InputIterator>
    RefineStrategies(InputIterator first,
                     InputIterator last,
                     const std::vector<std::unique_ptr<InterfaceFaceGhostCopy>>& ghost_copies)
        : RefinePatchStrategySet(first, last, false), d_ghost_copies(ghost_copies)
    {
    }

    void postprocessRefineBoxes(Patch<NDIM>& fine,
                                const Patch<NDIM>& coarse,
                                const BoxList<NDIM>& fine_boxes,
                                const IntVector<NDIM>& ratio) override
    {
        // The schedule also uses this object for temporary coarser levels, whose patches are not in the hierarchy.
        if (fine.inHierarchy())
        {
            for (const auto& ghost_copy : d_ghost_copies)
            {
                if (ghost_copy) ghost_copy->restore(fine);
            }
        }
        RefinePatchStrategySet::postprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }

private:
    const std::vector<std::unique_ptr<InterfaceFaceGhostCopy>>& d_ghost_copies;
};
} // namespace

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

    // Register the cubic coarsen operators with the grid geometry object.
    IBTK_DO_ONCE(d_grid_geom->addSpatialCoarsenOperator(new CartCellDoubleCubicCoarsen());
                 d_grid_geom->addSpatialCoarsenOperator(new CartSideDoubleCubicCoarsen()););

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
            Pointer<CoarsenOperator<NDIM>> coarsen_op = d_grid_geom->lookupCoarsenOperator(var, coarsen_op_name);
#if !defined(NDEBUG)
            TBOX_ASSERT(coarsen_op);
#endif
            d_coarsen_alg->registerCoarsen(src_data_idx, src_data_idx, coarsen_op);
            registered_coarsen_op = true;
        }
    }

    d_coarsen_strategy = nullptr;

    d_coarsen_scheds.resize(d_finest_ln + 1);
    if (registered_coarsen_op)
    {
        for (int src_ln = std::max(1, d_coarsest_ln); src_ln <= d_finest_ln; ++src_ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(src_ln);
            Pointer<PatchLevel<NDIM>> coarser_level = d_hierarchy->getPatchLevel(src_ln - 1);
            d_coarsen_scheds[src_ln] = d_coarsen_alg->createSchedule(coarser_level, level, d_coarsen_strategy.get());
        }
    }

    // Setup cached refine algorithms and schedules.
    d_data_kinds.clear();
    d_data_kinds.resize(d_transaction_comps.size());
    d_cf_fills.clear();
    d_cf_fills.resize(d_transaction_comps.size());
    d_face_ghost_copies.clear();
    d_face_ghost_copies.resize(d_transaction_comps.size());
    d_extrap_bc_ops.resize(d_transaction_comps.size());
    d_cc_robin_bc_ops.resize(d_transaction_comps.size());
    d_sc_robin_bc_ops.resize(d_transaction_comps.size());
    d_refine_alg = new RefineAlgorithm<NDIM>();
    std::vector<RefinePatchStrategy<NDIM>*> refine_patch_strategies;
    for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const int dst_data_idx = d_transaction_comps[comp_idx].d_dst_data_idx;
        const int src_data_idx = d_transaction_comps[comp_idx].d_src_data_idx;
        d_data_kinds[comp_idx] = { get_data_kind(dst_data_idx), get_data_kind(src_data_idx) };
        const std::string& phys_bdry_type = d_transaction_comps[comp_idx].d_phys_bdry_type;
        Pointer<Variable<NDIM>> var;
        var_db->mapIndexToVariable(src_data_idx, var);
        Pointer<CellVariable<NDIM, double>> cc_var = var;
        Pointer<NodeVariable<NDIM, double>> nc_var = var;
        Pointer<SideVariable<NDIM, double>> sc_var = var;
        Pointer<EdgeVariable<NDIM, double>> ec_var = var;
        Pointer<FaceVariable<NDIM, double>> fc_var = var;
        Pointer<RefineOperator<NDIM>> refine_op = nullptr;
        Pointer<VariableFillPattern<NDIM>> fill_pattern = d_transaction_comps[comp_idx].d_fill_pattern;
        check_in_place_fill_of_patch_border_data(
            d_transaction_comps[comp_idx], var, d_finest_ln, "initializeOperatorState");
        if (cc_var)
        {
            if (d_transaction_comps[comp_idx].d_refine_op_name != "NONE")
            {
                refine_op = d_grid_geom->lookupRefineOperator(var, d_transaction_comps[comp_idx].d_refine_op_name);
            }
            if (d_transaction_comps[comp_idx].d_use_cf_bdry_interpolation)
            {
                d_cf_fills[comp_idx] = std::make_unique<CFOwnerFill>(
                    d_hierarchy, d_coarsest_ln, d_finest_ln, dst_data_idx, src_data_idx, fill_pattern);
                refine_patch_strategies.push_back(d_cf_fills[comp_idx].get());
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
            // Every side-centered component gets the faces of the fine patches along the coarse-fine interface in the
            // ghost cells of the other patches, whatever the refine operator.
            d_face_ghost_copies[comp_idx] =
                std::make_unique<InterfaceFaceGhostCopy>(d_hierarchy,
                                                         d_coarsest_ln,
                                                         d_finest_ln,
                                                         dst_data_idx,
                                                         src_data_idx,
                                                         fill_pattern,
                                                         preserves_interface_faces(d_transaction_comps[comp_idx]));
            if (d_transaction_comps[comp_idx].d_use_cf_bdry_interpolation)
            {
                d_cf_fills[comp_idx] = std::make_unique<CFOwnerFill>(
                    d_hierarchy, d_coarsest_ln, d_finest_ln, dst_data_idx, src_data_idx, fill_pattern);
                refine_patch_strategies.push_back(d_cf_fills[comp_idx].get());
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
        if (!null_bc_coefs && sc_var && !d_transaction_comps[comp_idx].d_sc_phys_bdry_op)
        {
#if !defined(NDEBUG)
            TBOX_ASSERT(robin_bc_coefs.size() % NDIM == 0);
#endif
            d_sc_robin_bc_ops[comp_idx] =
                new CartSideRobinPhysBdryOp(dst_data_idx, robin_bc_coefs, d_homogeneous_bc, phys_bdry_type);
        }
    }

    d_refine_strategy = std::make_unique<RefineStrategies>(
        refine_patch_strategies.begin(), refine_patch_strategies.end(), d_face_ghost_copies);

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

    // The data and the coarse-fine fills are set up once, in initializeOperatorState().
    for (unsigned int comp_idx = 0; comp_idx < transaction_comps.size(); ++comp_idx)
    {
        const std::array<std::pair<int, int>, 2> kinds = { get_data_kind(transaction_comps[comp_idx].d_dst_data_idx),
                                                           get_data_kind(transaction_comps[comp_idx].d_src_data_idx) };
        for (unsigned int k = 0; k < 2; ++k)
        {
            if (kinds[k] != d_data_kinds[comp_idx][k])
            {
                TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                           << "  new transaction component " << comp_idx << " has "
                           << (k == 0 ? "destination" : "source") << " " << data_kind_to_string(kinds[k])
                           << ", but initializeOperatorState() set up "
                           << data_kind_to_string(d_data_kinds[comp_idx][k]) << "\n"
                           << "  call initializeOperatorState() with the new transaction components instead\n");
            }
        }
        const bool use_cf = transaction_comps[comp_idx].d_use_cf_bdry_interpolation;
        if (use_cf != d_transaction_comps[comp_idx].d_use_cf_bdry_interpolation)
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                       << "  new transaction component " << comp_idx << " sets d_use_cf_bdry_interpolation to "
                       << use_cf << ", but initializeOperatorState() set it to " << !use_cf << "\n"
                       << "  call initializeOperatorState() with the new transaction components instead\n");
        }

        // A reset must not remove the operator: initializeOperatorState() sets up the Robin fill of a side-centered
        // component only if the component has no physical boundary operator.
        if (d_transaction_comps[comp_idx].d_sc_phys_bdry_op && !transaction_comps[comp_idx].d_sc_phys_bdry_op)
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                       << "  new transaction component " << comp_idx
                       << " has no physical boundary operator, unlike the current one\n"
                       << "  call initializeOperatorState() with the new transaction components instead\n");
        }
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
            Pointer<CoarsenOperator<NDIM>> coarsen_op = d_grid_geom->lookupCoarsenOperator(var, coarsen_op_name);
#if !defined(NDEBUG)
            TBOX_ASSERT(coarsen_op);
#endif
            d_coarsen_alg->registerCoarsen(src_data_idx, src_data_idx, coarsen_op);
            registered_coarsen_op = true;
        }
    }

    if (registered_coarsen_op)
    {
        for (int src_ln = std::max(1, d_coarsest_ln); src_ln <= d_finest_ln; ++src_ln)
        {
            d_coarsen_alg->resetSchedule(d_coarsen_scheds[src_ln]);
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
        check_in_place_fill_of_patch_border_data(
            d_transaction_comps[comp_idx], var, d_finest_ln, "resetTransactionComponents");
        if (sc_var)
        {
            const bool protect = preserves_interface_faces(d_transaction_comps[comp_idx]);
            if (protect != d_face_ghost_copies[comp_idx]->protectsFaces())
            {
                TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                           << "  new transaction component " << comp_idx
                           << " changes whether it fills in place with a refine operator\n"
                           << "  call initializeOperatorState() with the new transaction components instead\n");
            }
            else
            {
                d_face_ghost_copies[comp_idx]->resetPatchDataIndices(dst_data_idx, src_data_idx, fill_pattern);
            }
        }
        if (d_cf_fills[comp_idx])
        {
            d_cf_fills[comp_idx]->resetPatchDataIndices(
                dst_data_idx, src_data_idx, d_transaction_comps[comp_idx].d_fill_pattern);
        }
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
    d_data_kinds.clear();
    d_cf_fills.clear();
    d_face_ghost_copies.clear();
    d_extrap_bc_ops.clear();
    d_cc_robin_bc_ops.clear();
    d_sc_robin_bc_ops.clear();

    // Clear cached communication schedules.
    d_coarsen_alg.setNull();
    d_coarsen_strategy = nullptr;
    d_coarsen_scheds.clear();

    d_refine_alg.setNull();
    d_refine_strategy = nullptr;
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
        if (d_coarsen_scheds[src_ln]) d_coarsen_scheds[src_ln]->coarsenData();
    }
    IBTK_TIMER_STOP(t_fill_data_coarsen);

    // Perform the initial data fill, using extrapolation to determine ghost
    // cell values at physical boundaries.
    IBTK_TIMER_START(t_fill_data_refine);
    for (int dst_ln = d_coarsest_ln; dst_ln <= d_finest_ln; ++dst_ln)
    {
        // A refine operator overwrites the faces of the patch boundaries along the coarse-fine interface. Restore their
        // values, which the coarse-fine fill also uses.
        for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
        {
            if (d_face_ghost_copies[comp_idx]) d_face_ghost_copies[comp_idx]->save(dst_ln);
        }
        // The refine schedule calls the coarse-fine fills to compute the tangential values.
        if (d_refine_scheds[dst_ln]) d_refine_scheds[dst_ln]->fillData(fill_time);
        for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
        {
            if (d_face_ghost_copies[comp_idx]) d_face_ghost_copies[comp_idx]->copyToGhostCells(dst_ln, fill_time);
            if (d_cf_fills[comp_idx] && d_cf_fills[comp_idx]->hasValues(dst_ln))
            {
                d_cf_fills[comp_idx]->computeValues(dst_ln);
                d_cf_fills[comp_idx]->copyValues(dst_ln, fill_time);
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
                    RobinPhysBdryPatchStrategy* const sc_phys_bdry_op = d_transaction_comps[comp_idx].d_sc_phys_bdry_op;
                    if (sc_phys_bdry_op)
                    {
                        const int dst_data_idx = d_transaction_comps[comp_idx].d_dst_data_idx;
                        const IntVector<NDIM>& ghost_width_to_fill =
                            patch->getPatchData(dst_data_idx)->getGhostCellWidth();
                        sc_phys_bdry_op->setPatchDataIndex(dst_data_idx);
                        sc_phys_bdry_op->setHomogeneousBc(d_homogeneous_bc);
                        sc_phys_bdry_op->setPhysicalBoundaryConditions(*patch, fill_time, ghost_width_to_fill);
                    }
                    else if (d_sc_robin_bc_ops[comp_idx])
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

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
