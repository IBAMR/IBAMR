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
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/CartSideDoubleCubicCoarsen.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/RefinePatchStrategySet.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/private/CoarseFineOwnerFill.h>
#include <ibtk/private/CubicCoarsenStencils.h>
#include <ibtk/private/InterfaceFaceGhostCopy.h>
#include <ibtk/private/data_kind.h>

#include <tbox/Pointer.h>
#include <tbox/Schedule.h>
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
#include <PatchDataFactory.h>
#include <PatchGeometry.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefinePatchStrategy.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideDataFactory.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableDatabase.h>
#include <VariableFillPattern.h>

#include <algorithm>
#include <limits>
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
} // preserves_interface_faces

// In an in-place refine fill, IBSAMRAI2#67 also overwrites the patch-border values of node, edge, and face data, and
// only side data are restored. Remove this check once IBAMR requires an IBSAMRAI2 release that fixes #67.
void
check_in_place_fill_of_patch_border_data(const HierarchyGhostCellInterpolation::InterpolationTransactionComponent& comp,
                                         const Pointer<Variable<NDIM>>& var,
                                         const int finest_ln)
{
    if (finest_ln < 1 || !preserves_interface_faces(comp))
    {
        return;
    }
    Pointer<SideVariable<NDIM, double>> sc_var = var;
    if (sc_var || !var->getPatchDataFactory()->dataLivesOnPatchBorder())
    {
        return;
    }
    TBOX_ERROR("HierarchyGhostCellInterpolation::check_in_place_fill_of_patch_border_data():\n"
               << "  variable " << var->getName() << " is filled in place with refine type " << comp.d_refine_op_name
               << " on a level with a coarser level, which overwrites its patch-border values\n"
               << "  fill it from a separate source patch data index\n");
} // check_in_place_fill_of_patch_border_data

// The refine strategies of the fill, which also restore the faces along the coarse-fine interface of a patch right
// after the refine operator fills it, before the schedule copies the values to other patches and sets the physical
// boundary values.
class RefineStrategies : public RefinePatchStrategySet
{
public:
    template <typename InputIterator>
    RefineStrategies(InputIterator first,
                     InputIterator last,
                     const std::vector<std::unique_ptr<detail::InterfaceFaceGhostCopy>>& ghost_copies)
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
                if (ghost_copy)
                {
                    ghost_copy->restore(fine);
                }
            }
        }
        RefinePatchStrategySet::postprocessRefineBoxes(fine, coarse, fine_boxes, ratio);
    }

private:
    const std::vector<std::unique_ptr<detail::InterfaceFaceGhostCopy>>& d_ghost_copies;
};

// Whether the component has Robin coefficients that its coarse data can be extended with.
bool
uses_robin_coefs(const HierarchyGhostCellInterpolation::InterpolationTransactionComponent& comp)
{
    return !comp.d_sc_phys_bdry_op &&
           std::any_of(comp.d_robin_bc_coefs.begin(),
                       comp.d_robin_bc_coefs.end(),
                       [](const RobinBcCoefStrategy<NDIM>* coef) { return coef != nullptr; });
} // uses_robin_coefs

// Whether the coarse data of a component are extended with the Robin ghost value.
bool
has_robin_extension(const HierarchyGhostCellInterpolation::InterpolationTransactionComponent& comp)
{
    return uses_robin_coefs(comp) && !comp.d_extend_coarse_data_linearly;
} // has_robin_extension

// Whether the coarse data of a component are extended linearly in place of the Robin ghost value.
bool
has_linear_extension(const HierarchyGhostCellInterpolation::InterpolationTransactionComponent& comp)
{
    return uses_robin_coefs(comp) && comp.d_extend_coarse_data_linearly;
} // has_linear_extension

// Determine whether the coarse data of the components are extended by quadratic extrapolation, and the stencil width
// that getRefineOpStencilWidth() returns.
void
determine_extension(const std::vector<HierarchyGhostCellInterpolation::InterpolationTransactionComponent>& comps,
                    bool& extrapolate_quadratically,
                    int& stencil_width)
{
    // Quadratic extrapolation uses three layers next to a boundary and requires stencil width 2, which the schedule
    // supports only if all destination data have that many ghost cells.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    bool has_extrapolated_component = false;
    bool wide_data = true;
    for (const auto& comp : comps)
    {
        has_extrapolated_component =
            has_extrapolated_component || (comp.d_use_cf_bdry_interpolation && !uses_robin_coefs(comp));
        wide_data =
            wide_data &&
            var_db->getPatchDescriptor()->getPatchDataFactory(comp.d_dst_data_idx)->getGhostCellWidth().min() >= 2;
    }
    extrapolate_quadratically = has_extrapolated_component && wide_data;
    stencil_width = extrapolate_quadratically ? 2 : 1;
    return;
} // determine_extension

// Return true if the cubic coarsening of a level with the refinement ratio uses a ghost layer of its source: every
// ratio is even, and some ratio is 2.
bool
ratio_uses_cubic_scratch(const IntVector<NDIM>& ratio)
{
    return detail::ratios_are_even(ratio) && ratio.min() < 4;
} // ratio_uses_cubic_scratch

// Gives the fine level to the cubic coarsen operators, which SAMRAI's coarsen() does not provide, and releases it when
// the object goes out of scope.
class CubicCoarsenFineLevel
{
public:
    CubicCoarsenFineLevel(CartCellDoubleCubicCoarsen& cc_op,
                          CartSideDoubleCubicCoarsen& sc_op,
                          const Pointer<PatchLevel<NDIM>>& fine_level)
        : d_cc_op(cc_op), d_sc_op(sc_op)
    {
        d_cc_op.setFineLevel(fine_level);
        d_sc_op.setFineLevel(fine_level);
    }

    ~CubicCoarsenFineLevel()
    {
        d_cc_op.setFineLevel(nullptr);
        d_sc_op.setFineLevel(nullptr);
    }

    CubicCoarsenFineLevel(const CubicCoarsenFineLevel&) = delete;
    CubicCoarsenFineLevel& operator=(const CubicCoarsenFineLevel&) = delete;

private:
    CartCellDoubleCubicCoarsen& d_cc_op;
    CartSideDoubleCubicCoarsen& d_sc_op;
};
} // namespace

// Extends the temporary coarse data of the coarse-fine interpolation beyond physical boundaries with the Robin ghost
// value, or with linear or quadratic extrapolation, as each component requires.
class HierarchyGhostCellInterpolation::CoarseDataExtension : public RefinePatchStrategy<NDIM>
{
public:
    explicit CoarseDataExtension(const HierarchyGhostCellInterpolation& ghost_fill);

    void setFillTime(const double fill_time)
    {
        d_fill_time = fill_time;
    }

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

    // Stop if the extension that the components determine differs from the one made at set-up.
    void checkReset(const std::vector<InterpolationTransactionComponent>& comps) const;

private:
    const HierarchyGhostCellInterpolation& d_ghost_fill;
    CartExtrapPhysBdryOp d_extrap_op;
    int d_stencil_width = 1;
    bool d_extrapolate_quadratically = false;
    double d_fill_time = 0.0;
};

HierarchyGhostCellInterpolation::CoarseDataExtension::CoarseDataExtension(
    const HierarchyGhostCellInterpolation& ghost_fill)
    : d_ghost_fill(ghost_fill)
{
    d_extrap_op.setExtrapolationType("QUADRATIC");
    determine_extension(d_ghost_fill.d_transaction_comps, d_extrapolate_quadratically, d_stencil_width);
    return;
} // CoarseDataExtension

void
HierarchyGhostCellInterpolation::CoarseDataExtension::checkReset(
    const std::vector<InterpolationTransactionComponent>& comps) const
{
    bool extrapolate_quadratically = false;
    int stencil_width = 1;
    determine_extension(comps, extrapolate_quadratically, stencil_width);
    if (extrapolate_quadratically != d_extrapolate_quadratically || stencil_width != d_stencil_width)
    {
        TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                   << "  the new transaction components require a different extension of the coarse data beyond "
                      "physical boundaries\n"
                   << "  call initializeOperatorState() with the new transaction components instead\n");
    }
    return;
} // checkReset

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
    return d_stencil_width;
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
    // Interpolation from coarse data occurs only on patches of the hierarchy.
    if (!fine.inHierarchy())
    {
        return;
    }

    Pointer<PatchGeometry<NDIM>> pgeom = coarse.getPatchGeometry();
    if (!pgeom->getTouchesRegularBoundary())
    {
        return;
    }

    // The data of the temporary coarse patch, which SAMRAI passes as const, are scratch data of the schedule.
    Patch<NDIM>& coarse_patch = const_cast<Patch<NDIM>&>(coarse);
    for (unsigned int comp_idx = 0; comp_idx < d_ghost_fill.d_transaction_comps.size(); ++comp_idx)
    {
        if (!d_ghost_fill.d_coarse_fine_fills[comp_idx])
        {
            continue;
        }
        const InterpolationTransactionComponent& comp = d_ghost_fill.d_transaction_comps[comp_idx];
        const int dst_data_idx = comp.d_dst_data_idx;
        Pointer<PatchData<NDIM>> coarse_data = coarse.getPatchData(dst_data_idx);
        if (has_robin_extension(comp))
        {
            // Computing the Robin value uses layers 0 and 1, which exist on every temporary patch next to a boundary.
            IntVector<NDIM> width = coarse_data->getGhostCellWidth();
            width.min(IntVector<NDIM>(1));
            if (Pointer<CellData<NDIM, double>>(coarse_data))
            {
                CartCellRobinPhysBdryOp op(
                    dst_data_idx, comp.d_robin_bc_coefs, d_ghost_fill.d_homogeneous_bc, "QUADRATIC");
                op.setPhysicalBoundaryConditions(coarse_patch, d_fill_time, width);
            }
            else
            {
                CartSideRobinPhysBdryOp op(
                    dst_data_idx, comp.d_robin_bc_coefs, d_ghost_fill.d_homogeneous_bc, "QUADRATIC");
                // Set the data at Dirichlet boundaries here: if the source differs from the destination, the pre-step
                // of fillData() sets these values only in the destination data.
                op.setDataAtDirichletBoundaries(coarse_patch, d_fill_time);
                op.setPhysicalBoundaryConditions(coarse_patch, d_fill_time, width);
            }
        }
        else if (has_linear_extension(comp))
        {
            // Linear extrapolation, like the computation of the Robin value, uses layers 0 and 1.
            IntVector<NDIM> width = coarse_data->getGhostCellWidth();
            width.min(IntVector<NDIM>(1));
            CartExtrapPhysBdryOp op(dst_data_idx, "LINEAR");
            op.setPhysicalBoundaryConditions(coarse_patch, d_fill_time, width);
        }
        else if (d_extrapolate_quadratically)
        {
            d_extrap_op.setPatchDataIndex(dst_data_idx);
            d_extrap_op.setPhysicalBoundaryConditions(coarse_patch, 0.0, coarse_data->getGhostCellWidth());
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

    // Register the IBTK transfer operators with the grid geometry object. The cubic coarsen operators among them serve
    // the objects that do not mark the invalid ghost values of their data.
    register_ibtk_transfer_operators(*d_grid_geom);

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
            Pointer<CoarsenOperator<NDIM>> coarsen_op = d_grid_geom->lookupCoarsenOperator(var, coarsen_op_name);
#if !defined(NDEBUG)
            TBOX_ASSERT(coarsen_op);
#endif
            d_coarsen_alg->registerCoarsen(src_data_idx, src_data_idx, coarsen_op);
            registered_coarsen_op = true;
        }
    }

    d_coarsen_strategy = nullptr;

    setUpCubicScratch(/*initialize*/ true);
    d_coarsen_scheds.resize(d_finest_ln + 1);
    d_cubic_scratch_fill_scheds.resize(d_finest_ln + 1);
    if (registered_coarsen_op)
    {
        for (int src_ln = std::max(1, d_coarsest_ln); src_ln <= d_finest_ln; ++src_ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(src_ln);
            Pointer<PatchLevel<NDIM>> coarser_level = d_hierarchy->getPatchLevel(src_ln - 1);
            if (levelUsesCubicScratch(src_ln))
            {
                d_coarsen_scheds[src_ln] =
                    d_cubic_scratch_coarsen_alg->createSchedule(coarser_level, level, d_coarsen_strategy.get());
                d_cubic_scratch_fill_scheds[src_ln] = d_cubic_scratch_fill_alg->createSchedule(level);
            }
            else
            {
                d_coarsen_scheds[src_ln] =
                    d_coarsen_alg->createSchedule(coarser_level, level, d_coarsen_strategy.get());
            }
        }
    }
    allocateCubicScratch();

    // Setup cached refine algorithms and schedules.
    d_data_kinds.clear();
    d_data_kinds.resize(d_transaction_comps.size());
    d_coarse_fine_fills.clear();
    d_coarse_fine_fills.resize(d_transaction_comps.size());
    d_face_ghost_copies.clear();
    d_face_ghost_copies.resize(d_transaction_comps.size());
    d_extrap_bc_ops.resize(d_transaction_comps.size());
    d_cc_robin_bc_ops.resize(d_transaction_comps.size());
    d_sc_robin_bc_ops.resize(d_transaction_comps.size());
    d_refine_scheds_built_in_place.resize(d_transaction_comps.size());
    d_refine_alg = new RefineAlgorithm<NDIM>();
    std::vector<RefinePatchStrategy<NDIM>*> refine_patch_strategies;
    for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
    {
        const int dst_data_idx = d_transaction_comps[comp_idx].d_dst_data_idx;
        const int src_data_idx = d_transaction_comps[comp_idx].d_src_data_idx;
        if (d_transaction_comps[comp_idx].d_consistent_type_2_bdry)
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::initializeOperatorState():\n"
                       << "  transaction component " << comp_idx << " sets consistent_type_2_bdry to true, but the "
                       << "option was removed and has no effect\n"
                       << "  set consistent_type_2_bdry to false\n");
        }
        d_data_kinds[comp_idx] = { detail::get_data_kind(dst_data_idx), detail::get_data_kind(src_data_idx) };
        d_refine_scheds_built_in_place[comp_idx] = src_data_idx == dst_data_idx;
        const std::string& phys_bdry_type = d_transaction_comps[comp_idx].d_phys_bdry_type;
        Pointer<Variable<NDIM>> var;
        var_db->mapIndexToVariable(src_data_idx, var);
        Pointer<CellVariable<NDIM, double>> cc_var = var;
        Pointer<NodeVariable<NDIM, double>> nc_var = var;
        Pointer<SideVariable<NDIM, double>> sc_var = var;
        Pointer<EdgeVariable<NDIM, double>> ec_var = var;
        Pointer<FaceVariable<NDIM, double>> fc_var = var;
        if (d_transaction_comps[comp_idx].d_extend_coarse_data_linearly && !cc_var)
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::initializeOperatorState():\n"
                       << "  d_extend_coarse_data_linearly is true for transaction component " << comp_idx
                       << ", whose data are not cell-centered\n"
                       << "  set d_extend_coarse_data_linearly only for cell-centered data\n");
        }
        Pointer<RefineOperator<NDIM>> refine_op = nullptr;
        Pointer<VariableFillPattern<NDIM>> fill_pattern = d_transaction_comps[comp_idx].d_fill_pattern;
        check_in_place_fill_of_patch_border_data(d_transaction_comps[comp_idx], var, d_finest_ln);
        if (cc_var)
        {
            if (d_transaction_comps[comp_idx].d_refine_op_name != "NONE")
            {
                refine_op = d_grid_geom->lookupRefineOperator(var, d_transaction_comps[comp_idx].d_refine_op_name);
            }
            if (d_transaction_comps[comp_idx].d_use_cf_bdry_interpolation)
            {
                d_coarse_fine_fills[comp_idx] = std::make_unique<detail::CoarseFineOwnerFill>(
                    d_hierarchy, d_coarsest_ln, d_finest_ln, dst_data_idx, src_data_idx, fill_pattern);
                refine_patch_strategies.push_back(d_coarse_fine_fills[comp_idx].get());
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
            d_face_ghost_copies[comp_idx] = std::make_unique<detail::InterfaceFaceGhostCopy>(
                d_hierarchy,
                d_coarsest_ln,
                d_finest_ln,
                dst_data_idx,
                src_data_idx,
                fill_pattern,
                preserves_interface_faces(d_transaction_comps[comp_idx]));
            if (d_transaction_comps[comp_idx].d_use_cf_bdry_interpolation)
            {
                d_coarse_fine_fills[comp_idx] = std::make_unique<detail::CoarseFineOwnerFill>(
                    d_hierarchy, d_coarsest_ln, d_finest_ln, dst_data_idx, src_data_idx, fill_pattern);
                refine_patch_strategies.push_back(d_coarse_fine_fills[comp_idx].get());
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

    // Extend the coarse data before the coarse-fine interpolations use them.
    d_coarse_extension = nullptr;
    if (std::any_of(d_coarse_fine_fills.begin(),
                    d_coarse_fine_fills.end(),
                    [](const auto& coarse_fine_fill) { return coarse_fine_fill != nullptr; }))
    {
        d_coarse_extension = std::make_unique<CoarseDataExtension>(*this);
        refine_patch_strategies.insert(refine_patch_strategies.begin(), d_coarse_extension.get());
    }

    d_refine_strategy = std::make_unique<RefineStrategies>(
        refine_patch_strategies.begin(), refine_patch_strategies.end(), d_face_ghost_copies);

    d_refine_scheds.resize(d_finest_ln + 1);
    for (int dst_ln = d_coarsest_ln; dst_ln <= d_finest_ln; ++dst_ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(dst_ln);
        d_refine_scheds[dst_ln] = d_refine_alg->createSchedule(level, dst_ln - 1, d_hierarchy, d_refine_strategy.get());
    }

    // One schedule for each level sends the faces along the coarse-fine interface and the coarse-fine values of all
    // components.
    d_interface_scheds.clear();
    d_interface_scheds.resize(d_finest_ln + 1);
    const bool has_interface_data =
        std::any_of(
            d_face_ghost_copies.begin(), d_face_ghost_copies.end(), [](const auto& c) { return c != nullptr; }) ||
        std::any_of(d_coarse_fine_fills.begin(), d_coarse_fine_fills.end(), [](const auto& f) { return f != nullptr; });
    for (int dst_ln = std::max(1, d_coarsest_ln); has_interface_data && dst_ln <= d_finest_ln; ++dst_ln)
    {
        d_interface_scheds[dst_ln] = std::make_unique<Schedule>();
        for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
        {
            if (d_face_ghost_copies[comp_idx])
            {
                d_face_ghost_copies[comp_idx]->appendTransactions(dst_ln, *d_interface_scheds[dst_ln]);
            }
            if (d_coarse_fine_fills[comp_idx])
            {
                d_coarse_fine_fills[comp_idx]->appendTransactions(dst_ln, *d_interface_scheds[dst_ln]);
            }
        }
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
        if (transaction_comps[comp_idx].d_consistent_type_2_bdry)
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                       << "  new transaction component " << comp_idx << " sets consistent_type_2_bdry to true, but the "
                       << "option was removed and has no effect\n"
                       << "  set consistent_type_2_bdry to false\n");
        }
        const std::array<detail::DataKind, 2> kinds = {
            detail::get_data_kind(transaction_comps[comp_idx].d_dst_data_idx),
            detail::get_data_kind(transaction_comps[comp_idx].d_src_data_idx)
        };
        for (unsigned int k = 0; k < 2; ++k)
        {
            if (kinds[k] != d_data_kinds[comp_idx][k])
            {
                TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                           << "  new transaction component " << comp_idx << " has "
                           << (k == 0 ? "destination" : "source") << " " << detail::data_kind_to_string(kinds[k])
                           << ", but initializeOperatorState() set up "
                           << detail::data_kind_to_string(d_data_kinds[comp_idx][k]) << "\n"
                           << "  call initializeOperatorState() with the new transaction components instead\n");
            }
        }
        const bool use_coarse_fine_interpolation = transaction_comps[comp_idx].d_use_cf_bdry_interpolation;
        if (use_coarse_fine_interpolation != d_transaction_comps[comp_idx].d_use_cf_bdry_interpolation)
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                       << "  new transaction component " << comp_idx << " sets d_use_cf_bdry_interpolation to "
                       << use_coarse_fine_interpolation << ", but initializeOperatorState() set it to "
                       << !use_coarse_fine_interpolation << "\n"
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

        // A reset cannot add or remove the extrapolation or Robin fill set up by initializeOperatorState().
        const std::string& phys_bdry_extrap_type = transaction_comps[comp_idx].d_phys_bdry_extrap_type;
        if ((phys_bdry_extrap_type != "NONE") != static_cast<bool>(d_extrap_bc_ops[comp_idx]))
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                       << "  new transaction component " << comp_idx << " sets d_phys_bdry_extrap_type to "
                       << phys_bdry_extrap_type << ", but initializeOperatorState() set up "
                       << (d_extrap_bc_ops[comp_idx] ? "an" : "no") << " extrapolation at physical boundaries\n"
                       << "  call initializeOperatorState() with the new transaction components instead\n");
        }
        const std::vector<RobinBcCoefStrategy<NDIM>*>& robin_bc_coefs = transaction_comps[comp_idx].d_robin_bc_coefs;
        const bool has_robin_bc_coefs =
            std::any_of(robin_bc_coefs.begin(),
                        robin_bc_coefs.end(),
                        [](const RobinBcCoefStrategy<NDIM>* const bc_coef) { return bc_coef != nullptr; });
        const bool has_sc_phys_bdry_op = transaction_comps[comp_idx].d_sc_phys_bdry_op != nullptr;
        const bool needs_robin_fill =
            has_robin_bc_coefs && (kinds[1].centering == detail::DataCentering::CELL ||
                                   (kinds[1].centering == detail::DataCentering::SIDE && !has_sc_phys_bdry_op));
        const bool has_robin_fill =
            d_cc_robin_bc_ops[comp_idx] || (d_sc_robin_bc_ops[comp_idx] && !has_sc_phys_bdry_op);
        if (needs_robin_fill != has_robin_fill)
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                       << "  new transaction component " << comp_idx << " has " << (has_robin_bc_coefs ? "" : "no ")
                       << "Robin coefficients, but initializeOperatorState() set up " << (has_robin_fill ? "a" : "no")
                       << " Robin fill at physical boundaries\n"
                       << "  call initializeOperatorState() with the new transaction components instead\n");
        }
    }

    // SAMRAI leaves out the copy within each patch from a refine schedule built with the source as the destination.
    for (unsigned int comp_idx = 0; comp_idx < transaction_comps.size(); ++comp_idx)
    {
        const int dst_data_idx = transaction_comps[comp_idx].d_dst_data_idx;
        const int src_data_idx = transaction_comps[comp_idx].d_src_data_idx;
        if (d_refine_scheds_built_in_place[comp_idx] && src_data_idx != dst_data_idx)
        {
            // Name each patch data by its variable and context, which do not depend on the data indices.
            VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
            const auto name_data = [var_db](const int data_idx)
            {
                Pointer<Variable<NDIM>> var;
                Pointer<VariableContext> ctx;
                var_db->mapIndexToVariableAndContext(data_idx, var, ctx);
                if (!var) return std::string("with no variable");
                return ctx ? var->getName() + " in context " + ctx->getName() : var->getName() + " with no context";
            };
            TBOX_ERROR(
                "HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                << "  new transaction component " << comp_idx << " has source data " << name_data(src_data_idx)
                << " and destination data " << name_data(dst_data_idx) << "\n"
                << "  the schedules of this component were built with equal source and destination data indices\n"
                << "  call initializeOperatorState() with the new transaction components instead\n");
        }
    }

    // The extension of the coarse data is set up once, in initializeOperatorState().
    if (d_coarse_extension)
    {
        d_coarse_extension->checkReset(transaction_comps);
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

    setUpCubicScratch(/*initialize*/ false);
    for (int src_ln = std::max(1, d_coarsest_ln); src_ln <= d_finest_ln; ++src_ln)
    {
        if (levelUsesCubicScratch(src_ln) != static_cast<bool>(d_cubic_scratch_fill_scheds[src_ln]))
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                       << "  the new transaction components require a different cubic coarsening of a level whose "
                          "refinement ratios are all even and include a ratio of 2\n"
                       << "  call initializeOperatorState() with the new transaction components instead\n");
        }
    }
    if (registered_coarsen_op)
    {
        for (int src_ln = std::max(1, d_coarsest_ln); src_ln <= d_finest_ln; ++src_ln)
        {
            if (levelUsesCubicScratch(src_ln))
            {
                d_cubic_scratch_coarsen_alg->resetSchedule(d_coarsen_scheds[src_ln]);
                d_cubic_scratch_fill_alg->resetSchedule(d_cubic_scratch_fill_scheds[src_ln]);
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
        if (d_transaction_comps[comp_idx].d_extend_coarse_data_linearly && !cc_var)
        {
            TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                       << "  d_extend_coarse_data_linearly is true for transaction component " << comp_idx
                       << ", whose data are not cell-centered\n"
                       << "  set d_extend_coarse_data_linearly only for cell-centered data\n");
        }
        Pointer<RefineOperator<NDIM>> refine_op = nullptr;
        Pointer<VariableFillPattern<NDIM>> fill_pattern = d_transaction_comps[comp_idx].d_fill_pattern;
        check_in_place_fill_of_patch_border_data(d_transaction_comps[comp_idx], var, d_finest_ln);
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
        if (d_coarse_fine_fills[comp_idx])
        {
            d_coarse_fine_fills[comp_idx]->resetPatchDataIndices(
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
            d_extrap_bc_ops[comp_idx]->setPatchDataIndex(dst_data_idx);
            d_extrap_bc_ops[comp_idx]->setExtrapolationType(phys_bdry_extrap_type);
        }
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

    // Clear cached refinement operators. The schedules refer to the coarse-fine fills.
    d_interface_scheds.clear();
    d_data_kinds.clear();
    d_coarse_fine_fills.clear();
    d_face_ghost_copies.clear();
    d_extrap_bc_ops.clear();
    d_cc_robin_bc_ops.clear();
    d_sc_robin_bc_ops.clear();
    d_refine_scheds_built_in_place.clear();

    // Clear cached communication schedules.
    d_coarsen_alg.setNull();
    d_coarsen_strategy = nullptr;
    d_coarsen_scheds.clear();
    deallocateCubicScratch();
    d_cubic_scratch_coarsen_alg.setNull();
    d_cubic_scratch_fill_alg.setNull();
    d_cubic_scratch_fill_scheds.clear();
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
    if (d_coarse_extension)
    {
        d_coarse_extension->setFillTime(fill_time);
    }
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

    // In an in-place fill, set the side-centered data at Dirichlet boundaries before the coarsening and the
    // coarse-fine interpolation use them. A fill of level 0 alone uses none of them.
    if (d_finest_ln > 0)
    {
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                if (!patch->getPatchGeometry()->getTouchesRegularBoundary())
                {
                    continue;
                }
                for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
                {
                    if (d_sc_robin_bc_ops[comp_idx] && !d_transaction_comps[comp_idx].d_sc_phys_bdry_op)
                    {
                        d_sc_robin_bc_ops[comp_idx]->setDataAtDirichletBoundaries(*patch, fill_time);
                    }
                }
            }
        }
    }

    // Synchronize data on the patch hierarchy prior to filling ghost cell
    // values.
    IBTK_TIMER_START(t_fill_data_coarsen);
    for (int src_ln = d_finest_ln; src_ln >= std::max(1, d_coarsest_ln); --src_ln)
    {
        if (d_cubic_scratch_fill_scheds[src_ln])
        {
            // Copy the source data of the level into the scratch data. The cubic coarsen operators use the boxes of
            // the level.
            d_cubic_scratch_fill_scheds[src_ln]->fillData(fill_time, false);
            const CubicCoarsenFineLevel fine_level_scope(
                *d_cc_cubic_coarsen_op, *d_sc_cubic_coarsen_op, d_hierarchy->getPatchLevel(src_ln));
            d_coarsen_scheds[src_ln]->coarsenData();
        }
        else if (d_coarsen_scheds[src_ln])
        {
            d_coarsen_scheds[src_ln]->coarsenData();
        }
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
            if (d_face_ghost_copies[comp_idx])
            {
                d_face_ghost_copies[comp_idx]->save(dst_ln);
            }
        }
        // The refine schedule calls the coarse-fine fills to compute the tangential values.
        if (d_refine_scheds[dst_ln])
        {
            d_refine_scheds[dst_ln]->fillData(fill_time);
        }
        if (d_interface_scheds[dst_ln])
        {
            for (const auto& coarse_fine_fill : d_coarse_fine_fills)
            {
                if (coarse_fine_fill)
                {
                    coarse_fine_fill->computeValues(dst_ln);
                }
            }
            d_interface_scheds[dst_ln]->communicate();
            for (const auto& coarse_fine_fill : d_coarse_fine_fills)
            {
                if (coarse_fine_fill)
                {
                    coarse_fine_fill->copyValues(dst_ln);
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

Pointer<CoarsenOperator<NDIM>>
HierarchyGhostCellInterpolation::getCoarsenOperator(const Pointer<Variable<NDIM>>& var,
                                                    const InterpolationTransactionComponent& comp) const
{
    if (componentUsesCubicScratch(comp))
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
HierarchyGhostCellInterpolation::componentUsesCubicScratch(const InterpolationTransactionComponent& comp) const
{
    // An in-place fill with a refine operator overwrites the fine values that the ratio 2 stencil involves. Such
    // components are therefore coarsened with the registered operator.
    return comp.d_coarsen_op_name == "CUBIC_COARSEN" && comp.d_refine_op_name == "NONE";
} // componentUsesCubicScratch

void
HierarchyGhostCellInterpolation::setUpCubicScratch(const bool initialize)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_cubic_scratch_coarsen_alg.setNull();
    d_cubic_scratch_fill_alg.setNull();

    // A reset keeps the scratch data of initializeOperatorState(), which depend only on the components that use them
    // and on the patch layout.
    if (!initialize)
    {
        for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
        {
            if (componentUsesCubicScratch(d_transaction_comps[comp_idx]) !=
                (d_cubic_scratch_idxs[comp_idx] != invalid_index))
            {
                TBOX_ERROR("HierarchyGhostCellInterpolation::resetTransactionComponents():\n"
                           << "  new transaction component " << comp_idx
                           << " changes whether it is coarsened cubically with refine type NONE\n"
                           << "  call initializeOperatorState() with the new transaction components instead\n");
            }
        }
    }
    else
    {
        // Give each source of cubic coarsening a scratch data index with one ghost layer of its own, cloned from an
        // index shared by all objects, and removed by deallocateOperatorState().
        d_cubic_scratch_idxs.assign(d_transaction_comps.size(), invalid_index);
        bool registered_scratch = false;
        for (unsigned int comp_idx = 0; comp_idx < d_transaction_comps.size(); ++comp_idx)
        {
            const int src_data_idx = d_transaction_comps[comp_idx].d_src_data_idx;
            if (!componentUsesCubicScratch(d_transaction_comps[comp_idx]))
            {
                continue;
            }
            Pointer<Variable<NDIM>> var;
            var_db->mapIndexToVariable(src_data_idx, var);
            Pointer<CellVariable<NDIM, double>> cc_var = var;
            Pointer<SideVariable<NDIM, double>> sc_var = var;
            Pointer<PatchDataFactory<NDIM>> factory = var->getPatchDataFactory();
            Pointer<CellDataFactory<NDIM, double>> cc_factory = factory;
            Pointer<SideDataFactory<NDIM, double>> sc_factory = factory;
            const int depth = cc_var ? cc_factory->getDefaultDepth() : sc_factory->getDefaultDepth();
            const std::string name = std::string("IBTK::HierarchyGhostCellInterpolation::cubic_scratch_") +
                                     (cc_var ? "cc" : "sc") + std::to_string(depth);
            Pointer<Variable<NDIM>> scratch_var =
                var_db->checkVariableExists(name) ? var_db->getVariable(name) : nullptr;
            if (!scratch_var && cc_var)
            {
                scratch_var = new CellVariable<NDIM, double>(name, depth);
            }
            else if (!scratch_var && sc_var)
            {
                scratch_var = new SideVariable<NDIM, double>(name, depth);
            }
            Pointer<VariableContext> ctx = var_db->getContext("IBTK::HierarchyGhostCellInterpolation::CUBIC_COARSEN");
            const int shared_idx = var_db->registerVariableAndContext(scratch_var, ctx, IntVector<NDIM>(1));
            d_cubic_scratch_idxs[comp_idx] = var_db->registerClonedPatchDataIndex(scratch_var, shared_idx);
            registered_scratch = true;
        }
        if (!registered_scratch)
        {
            return;
        }

        // The stencil of a coarse cell must lie in the ghost layer of the fine patch over it.
        for (int ln = std::max(1, d_coarsest_ln); ln <= d_finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            const IntVector<NDIM>& ratio = level->getRatioToCoarserLevel();
            if (!ratio_uses_cubic_scratch(ratio))
            {
                continue;
            }
            const BoxArray<NDIM>& boxes = level->getBoxes();
            for (int k = 0; k < boxes.size(); ++k)
            {
                const Box<NDIM> aligned = Box<NDIM>::refine(Box<NDIM>::coarsen(boxes[k], ratio), ratio);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    if (ratio(d) < 4 &&
                        (aligned.lower(d) != boxes[k].lower(d) || aligned.upper(d) != boxes[k].upper(d)))
                    {
                        TBOX_ERROR(
                            "HierarchyGhostCellInterpolation::setUpCubicScratch():\n"
                            << "  patch " << boxes[k] << " of level " << ln
                            << " is not a whole number of coarse cells along axis " << d << "\n"
                            << "  cubic coarsening requires patches that start and end on coarse cell boundaries\n");
                    }
                }
            }
        }
    }
    if (std::none_of(
            d_cubic_scratch_idxs.begin(), d_cubic_scratch_idxs.end(), [](int idx) { return idx != invalid_index; }))
    {
        return;
    }

    // Coarsen the scratch data in place of the source data. Fill the scratch data from the source data of the
    // same level.
    d_cubic_scratch_coarsen_alg = new CoarsenAlgorithm<NDIM>();
    d_cubic_scratch_fill_alg = new RefineAlgorithm<NDIM>();
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
            d_cubic_scratch_coarsen_alg->registerCoarsen(src_data_idx, src_data_idx, coarsen_op);
        }
        else
        {
            d_cubic_scratch_coarsen_alg->registerCoarsen(src_data_idx, scratch_idx, coarsen_op);
            d_cubic_scratch_fill_alg->registerRefine(scratch_idx, src_data_idx, scratch_idx, nullptr);
        }
    }
    return;
} // setUpCubicScratch

bool
HierarchyGhostCellInterpolation::levelUsesCubicScratch(const int ln) const
{
    if (!d_cubic_scratch_coarsen_alg)
    {
        return false;
    }
    Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
    return ratio_uses_cubic_scratch(level->getRatioToCoarserLevel());
} // levelUsesCubicScratch

void
HierarchyGhostCellInterpolation::allocateCubicScratch()
{
    // The ghost values that no patch of the level covers are NaN. The fills set the other values.
    for (int ln = std::max(1, d_coarsest_ln); ln <= d_finest_ln; ++ln)
    {
        if (!d_cubic_scratch_fill_scheds[ln])
        {
            continue;
        }
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (const int scratch_idx : d_cubic_scratch_idxs)
        {
            if (scratch_idx == invalid_index)
            {
                continue;
            }
            if (!level->checkAllocated(scratch_idx))
            {
                level->allocatePatchData(scratch_idx);
            }
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> cc_scratch = patch->getPatchData(scratch_idx);
                Pointer<SideData<NDIM, double>> sc_scratch = patch->getPatchData(scratch_idx);
                if (cc_scratch)
                {
                    cc_scratch->fillAll(std::numeric_limits<double>::quiet_NaN());
                }
                else
                {
                    sc_scratch->fillAll(std::numeric_limits<double>::quiet_NaN());
                }
            }
        }
    }
    return;
} // allocateCubicScratch

void
HierarchyGhostCellInterpolation::deallocateCubicScratch()
{
    for (int ln = std::max(1, d_coarsest_ln); ln <= std::min(d_finest_ln, d_hierarchy->getFinestLevelNumber()); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
        for (const int scratch_idx : d_cubic_scratch_idxs)
        {
            if (scratch_idx != invalid_index && level->checkAllocated(scratch_idx))
            {
                level->deallocatePatchData(scratch_idx);
            }
        }
    }
    for (const int scratch_idx : d_cubic_scratch_idxs)
    {
        if (scratch_idx != invalid_index)
        {
            VariableDatabase<NDIM>::getDatabase()->removePatchDataIndex(scratch_idx);
        }
    }
    return;
} // deallocateCubicScratch

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
