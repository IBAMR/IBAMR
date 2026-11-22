// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2026 by the IBAMR developers
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
#include <ibamr/INSVCStaggeredConservativeMassMomentumSSPRKIntegrator.h>
#include <ibamr/ibamr_utilities.h>

#include <ibtk/CartGridFunction.h>
#include <ibtk/ExtendedRobinBcCoefStrategy.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

#include <ArrayData.h>
#include <BasePatchHierarchy.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <CartesianPatchGeometry.h>
#include <CoarseFineBoundary.h>
#include <FaceData.h>
#include <HierarchyDataOpsManager.h>
#include <HierarchySideDataOpsReal.h>
#include <Index.h>
#include <IntVector.h>
#include <MultiblockDataTranslator.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RobinBcCoefStrategy.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <Variable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

namespace SAMRAI
{
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Timer* t_apply_convective_operator;
static Timer* t_integrate;
static Timer* t_initialize_integrator;
static Timer* t_deallocate_integrator;

bool
rho_bcs_preserve_constant(const Pointer<PatchHierarchy<NDIM>>& hierarchy,
                          const int coarsest_ln,
                          const int finest_ln,
                          const int rho_scratch_idx,
                          const Pointer<Variable<NDIM>>& rho_var,
                          const std::vector<RobinBcCoefStrategy<NDIM>*>& rho_bc_coefs,
                          const std::string& extrap_type,
                          const double rho_const,
                          const std::array<double, 3>& eval_times,
                          const double rel_tol)
{
    if (!rho_var)
    {
        return false;
    }
    bool local_ok = true;
    const bool valid_extrap = (extrap_type == "CONSTANT" || extrap_type == "LINEAR" || extrap_type == "QUADRATIC");

    for (int ln = coarsest_ln; ln <= finest_ln && local_ok; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p && local_ok; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Array<BoundaryBox<NDIM>> physical_codim1_boxes =
                PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
            const int n_physical_codim1_boxes = physical_codim1_boxes.size();
            if (n_physical_codim1_boxes == 0)
            {
                continue;
            }

            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM> gcw_to_fill = patch->getPatchData(rho_scratch_idx)->getGhostCellWidth();
            const double* const dx = pgeom->getDx();
            const double* const patch_x_lower = pgeom->getXLower();
            const double* const patch_x_upper = pgeom->getXUpper();
            const IntVector<NDIM>& ratio_to_level_zero = pgeom->getRatio();
            Array<Array<bool>> touches_regular_bdry(NDIM), touches_periodic_bdry(NDIM);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                touches_regular_bdry[d].resizeArray(2);
                touches_periodic_bdry[d].resizeArray(2);
                for (int ul = 0; ul < 2; ++ul)
                {
                    touches_regular_bdry[d][ul] = pgeom->getTouchesRegularBoundary(d, ul);
                    touches_periodic_bdry[d][ul] = pgeom->getTouchesPeriodicBoundary(d, ul);
                }
            }

            for (unsigned int axis = 0; axis < NDIM && local_ok; ++axis)
            {
                RobinBcCoefStrategy<NDIM>* const bc_coef = (axis < rho_bc_coefs.size()) ? rho_bc_coefs[axis] : nullptr;
                if (!bc_coef)
                {
                    if (!valid_extrap)
                    {
                        local_ok = false;
                    }
                    continue;
                }
                if (dynamic_cast<ExtendedRobinBcCoefStrategy*>(bc_coef) != nullptr)
                {
                    local_ok = false;
                    break;
                }

                for (int n = 0; n < n_physical_codim1_boxes && local_ok; ++n)
                {
                    const BoundaryBox<NDIM>& bdry_box = physical_codim1_boxes[n];
                    const unsigned int location_index = bdry_box.getLocationIndex();
                    const unsigned int bdry_normal_axis = location_index / 2;
                    const double h = dx[bdry_normal_axis];
                    const int gcw_n = std::max(1, gcw_to_fill(bdry_normal_axis));
                    const Box<NDIM> bc_fill_box = pgeom->getBoundaryFillBox(bdry_box, patch_box, gcw_to_fill);
                    const BoundaryBox<NDIM> trimmed_bdry_box(
                        bdry_box.getBox() * bc_fill_box, bdry_box.getBoundaryType(), location_index);
                    if (trimmed_bdry_box.getBox().empty())
                    {
                        continue;
                    }
                    if (!std::isfinite(h) || h <= 0.0)
                    {
                        local_ok = false;
                        break;
                    }
                    Box<NDIM> bc_coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(trimmed_bdry_box);
                    const bool shifted = (bdry_normal_axis != axis);
                    if (shifted)
                    {
                        bc_coef_box.upper()(axis) += 1;
                        std::array<double, NDIM> shifted_lower, shifted_upper;
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            shifted_lower[d] = patch_x_lower[d];
                            shifted_upper[d] = patch_x_upper[d];
                        }
                        shifted_lower[axis] -= 0.5 * dx[axis];
                        shifted_upper[axis] -= 0.5 * dx[axis];
                        patch->setPatchGeometry(new CartesianPatchGeometry<NDIM>(ratio_to_level_zero,
                                                                                 touches_regular_bdry,
                                                                                 touches_periodic_bdry,
                                                                                 dx,
                                                                                 shifted_lower.data(),
                                                                                 shifted_upper.data()));
                    }

                    Pointer<ArrayData<NDIM, double>> acoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double>> bcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);
                    Pointer<ArrayData<NDIM, double>> gcoef_data = new ArrayData<NDIM, double>(bc_coef_box, 1);

                    for (double t_eval : eval_times)
                    {
                        bc_coef->setBcCoefs(
                            acoef_data, bcoef_data, gcoef_data, rho_var, *patch, trimmed_bdry_box, t_eval);
                        for (Box<NDIM>::Iterator bc(bc_coef_box); bc; bc++)
                        {
                            const hier::Index<NDIM>& idx = bc();
                            const double a = (*acoef_data)(idx, 0);
                            const double b = (*bcoef_data)(idx, 0);
                            const double g = (*gcoef_data)(idx, 0);
                            if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(g) || (a == 0.0 && b == 0.0))
                            {
                                local_ok = false;
                                break;
                            }
                            if (!shifted)
                            {
                                // Match cartphysbdryop3d.f.m4 side-normal branch (|b| < 1e-12 -> Dirichlet u_b = g /
                                // a).
                                if (std::abs(b) < 1.0e-12)
                                {
                                    if (a == 0.0)
                                    {
                                        local_ok = false;
                                        break;
                                    }
                                    const double u_b = g / a;
                                    const double u_scale = std::max(std::abs(rho_const), std::abs(u_b));
                                    if (!std::isfinite(u_b) || std::abs(u_b - rho_const) > rel_tol * u_scale)
                                    {
                                        local_ok = false;
                                        break;
                                    }
                                }
                                else
                                {
                                    const double res = std::abs(a * rho_const - g);
                                    const double rhs_scale = std::max(std::abs(a * rho_const), std::abs(g));
                                    const double ghost_err = (2.0 * static_cast<double>(gcw_n) * h / std::abs(b)) * res;
                                    const double rho_scale = std::max(1.0, std::abs(rho_const));
                                    if (res > rel_tol * rhs_scale || !std::isfinite(ghost_err) ||
                                        ghost_err > rel_tol * rho_scale)
                                    {
                                        local_ok = false;
                                        break;
                                    }
                                }
                            }
                            else
                            {
                                // Match cartphysbdryop3d.f.m4 tangential (shifted cell-centered) Robin branch.
                                const double res = std::abs(a * rho_const - g);
                                const double rhs_scale = std::max(std::abs(a * rho_const), std::abs(g));
                                if (res > rel_tol * rhs_scale)
                                {
                                    local_ok = false;
                                    break;
                                }
                                const double rho_scale = std::max(1.0, std::abs(rho_const));
                                for (int i = 0; i < gcw_n; ++i)
                                {
                                    const double n_off = 1.0 + 2.0 * static_cast<double>(i);
                                    const double denom = a * n_off * h + 2.0 * b;
                                    if (!std::isfinite(denom) || denom == 0.0)
                                    {
                                        local_ok = false;
                                        break;
                                    }
                                    const double ghost_err = std::abs((2.0 * n_off * h / denom) * (g - a * rho_const));
                                    if (!std::isfinite(ghost_err) || ghost_err > rel_tol * rho_scale)
                                    {
                                        local_ok = false;
                                        break;
                                    }
                                }
                                if (!local_ok)
                                {
                                    break;
                                }
                            }
                        }
                        if (!local_ok)
                        {
                            break;
                        }
                    }

                    if (shifted)
                    {
                        patch->setPatchGeometry(pgeom);
                    }
                }
            }
        }
    }

    return IBTK_MPI::minReduction(local_ok ? 1 : 0) != 0;
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::INSVCStaggeredConservativeMassMomentumSSPRKIntegrator(
    std::string object_name,
    Pointer<Database> input_db)
    : INSVCStaggeredConservativeMassMomentumRKIntegrator(object_name, input_db)
{
    switch (d_density_time_stepping_type)
    {
    case SSPRK2:
        d_num_steps = 2;
        break;
    case SSPRK3:
        d_num_steps = 3;
        break;
    default:
        TBOX_ERROR(
            "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::"
            "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator():\n"
            << "  unsupported density time stepping type: "
            << IBAMR::enum_to_string<TimeSteppingType>(d_density_time_stepping_type) << " \n"
            << "  valid choices are: SSPRK2 and SSPRK3\n");
    }

    if (input_db && input_db->keyExists("skip_constant_density_ssprk_stages"))
    {
        d_skip_const_density_stages = input_db->getBool("skip_constant_density_ssprk_stages");
    }
    if (input_db && input_db->keyExists("constant_density_max_stage_div_dt"))
    {
        d_const_density_max_stage_div_dt = input_db->getDouble("constant_density_max_stage_div_dt");
        if (!std::isfinite(d_const_density_max_stage_div_dt) || d_const_density_max_stage_div_dt < 0.0)
        {
            TBOX_ERROR(
                "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::"
                "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator():\n"
                << "  constant_density_max_stage_div_dt must be finite and >= 0.0, got "
                << d_const_density_max_stage_div_dt << "\n");
        }
    }

    // Setup Timers.
    IBAMR_DO_ONCE(t_apply_convective_operator = TimerManager::getManager()->getTimer(
                      "IBAMR::INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::"
                      "applyConvectiveOperator()");
                  t_integrate = TimerManager::getManager()->getTimer(
                      "IBAMR::INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate("
                      ")");
                  t_initialize_integrator = TimerManager::getManager()->getTimer(
                      "IBAMR::INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::"
                      "initializeSTSIntegrator()");
                  t_deallocate_integrator = TimerManager::getManager()->getTimer(
                      "IBAMR::INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::"
                      "deallocateSTSIntegrator()"););
    return;
} // INSVCStaggeredConservativeMassMomentumSSPRKIntegrator

INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::~INSVCStaggeredConservativeMassMomentumSSPRKIntegrator()
{
    if (d_skip_const_density_stages && (d_const_density_fastpath_hits > 0 || d_const_density_fallback_count > 0) &&
        IBTK_MPI::getRank() == 0)
    {
        std::ostringstream oss;
        oss << d_object_name << ": constant-density fastpath summary (hits=" << d_const_density_fastpath_hits
            << ", fallbacks=" << d_const_density_fallback_count << ", max_stage_div_dt=" << d_max_observed_stage_div_dt
            << ")\n";
        plog << oss.str();
        pout << oss.str();
    }
    return;
} // ~INSVCStaggeredConservativeMassMomentumSSPRKIntegrator

void
INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate(double dt)
{
    // Get hierarchy operation object
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_sc_data_ops =
        hier_ops_manager->getOperationsDouble(new SideVariable<NDIM, double>("sc_var"), d_hierarchy, true);

    IBAMR_TIMER_START(t_integrate)
#if !defined(NDEBUG)
    if (!d_is_initialized)
    {
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate():\n"
                   << "  time integrator must be initialized prior to call to "
                      "integrate()\n");
    }

    TBOX_ASSERT(d_rho_current_idx >= 0);
    TBOX_ASSERT(d_V_old_idx >= 0);
    TBOX_ASSERT(d_V_current_idx >= 0);
    TBOX_ASSERT(d_V_new_idx >= 0);
#endif

#if !defined(NDEBUG)
    TBOX_ASSERT(IBTK::rel_equal_eps(dt, getTimeStepSize()));
#endif

    if (d_V_old_idx == d_V_current_idx)
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(d_dt_prev <= 0.0);
#endif
        // Ensure that previous time step is set for initial times
        d_dt_prev = dt;
    }
#if !defined(NDEBUG)
    if (!(dt > 0.0))
    {
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate():\n"
                   << " invalid time step size dt = " << dt << "\n");
    }
#endif

// Assertions for velocity interpolation and extrapolation
#if !defined(NDEBUG)
    if (d_cycle_num < 0)
    {
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate():\n"
                   << "  invalid cycle number = " << d_cycle_num << "\n");
    }
    if (d_dt_prev <= 0.0)
    {
        TBOX_ERROR("INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate():\n"
                   << "  invalid previous time step size = " << d_dt_prev << "\n");
    }
#endif

    // Fill ghost cell values
    static const bool homogeneous_bc = false;
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;

    if (d_skip_const_density_stages)
    {
        int local_nonfinite = 0;
        long long local_count = 0;
        double rho_min = std::numeric_limits<double>::max();
        double rho_max = -std::numeric_limits<double>::max();
        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<SideData<NDIM, double>> R_cur_data = patch->getPatchData(d_rho_current_idx);
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
                    {
                        const SideIndex<NDIM> si(b(), axis, SideIndex<NDIM>::Lower);
                        const double val = (*R_cur_data)(si, 0);
                        if (!std::isfinite(val))
                        {
                            local_nonfinite = 1;
                            continue;
                        }
                        ++local_count;
                        rho_min = std::min(rho_min, val);
                        rho_max = std::max(rho_max, val);
                    }
                }
            }
        }
        const int global_nonfinite = IBTK_MPI::maxReduction(local_nonfinite);
        if (global_nonfinite != 0)
        {
            TBOX_ERROR(
                "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate(): "
                "non-finite density encountered with skip_constant_density_ssprk_stages enabled.\n");
        }
        const int global_has_samples = IBTK_MPI::maxReduction(local_count > 0 ? 1 : 0);
        rho_min = IBTK_MPI::minReduction(local_count > 0 ? rho_min : std::numeric_limits<double>::max());
        rho_max = IBTK_MPI::maxReduction(local_count > 0 ? rho_max : -std::numeric_limits<double>::max());
        const double rho_const = 0.5 * (rho_min + rho_max);
        const double scale = std::max(1.0, std::abs(rho_const));
        const double tol = 1.0e-12 * scale;
        const std::array<double, 3> bc_eval_times = {
            d_current_time,
            d_current_time + 0.5 * dt,
            d_current_time + dt,
        };

        bool can_fast_path = (!d_S_fcn) && (d_density_time_stepping_type == SSPRK3) && (global_has_samples != 0) &&
                             std::isfinite(rho_min) && std::isfinite(rho_max) && ((rho_max - rho_min) <= tol);
        Pointer<Variable<NDIM>> rho_var;
        if (can_fast_path)
        {
            VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(d_rho_scratch_idx, rho_var);
            can_fast_path = static_cast<bool>(rho_var) && rho_bcs_preserve_constant(d_hierarchy,
                                                                                    d_coarsest_ln,
                                                                                    d_finest_ln,
                                                                                    d_rho_scratch_idx,
                                                                                    rho_var,
                                                                                    d_rho_bc_coefs,
                                                                                    d_density_bdry_extrap_type,
                                                                                    rho_const,
                                                                                    bc_eval_times,
                                                                                    1.0e-12);
        }

        const double omega = dt / d_dt_prev;
        const double sum_dt = dt + d_dt_prev;
        struct StageVelSpec
        {
            double eval_time;
            double w0;
            double w1;
            double w2;
        };
        const std::array<StageVelSpec, 3> stage_specs = { {
            { d_current_time, 0.0, 1.0, 0.0 },
            { d_current_time + dt,
              (d_cycle_num > 0) ? 0.0 : (-omega),
              (d_cycle_num > 0) ? 0.0 : (1.0 + omega),
              (d_cycle_num > 0) ? 1.0 : 0.0 },
            { d_current_time + 0.5 * dt,
              (d_cycle_num > 0) ? (-0.25 * dt * dt / (d_dt_prev * sum_dt)) : (-0.5 * omega),
              (d_cycle_num > 0) ? (0.25 * (2.0 + omega)) : (1.0 + 0.5 * omega),
              (d_cycle_num > 0) ? (0.25 * (dt + 2.0 * d_dt_prev) / sum_dt) : 0.0 },
        } };

        for (int stage = 0; stage < 3 && can_fast_path; ++stage)
        {
            const StageVelSpec& spec = stage_specs[stage];
            if (stage == 0)
            {
                d_hier_sc_data_ops->copyData(d_V_composite_idx, d_V_current_idx, /*interior_only*/ true);
            }
            else
            {
                d_hier_sc_data_ops->linearSum(
                    d_V_composite_idx, spec.w0, d_V_old_idx, spec.w1, d_V_current_idx, /*interior_only*/ true);
                d_hier_sc_data_ops->axpy(
                    d_V_composite_idx, spec.w2, d_V_new_idx, d_V_composite_idx, /*interior_only*/ true);
            }
            std::vector<InterpolationTransactionComponent> v_stage_transaction_comps({ { d_V_scratch_idx,
                                                                                         d_V_composite_idx,
                                                                                         "CONSERVATIVE_LINEAR_REFINE",
                                                                                         false,
                                                                                         "CONSERVATIVE_COARSEN",
                                                                                         d_velocity_bdry_extrap_type,
                                                                                         false,
                                                                                         d_u_bc_coefs } });
            d_hier_v_bdry_fill->resetTransactionComponents(v_stage_transaction_comps);
            StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
                d_u_bc_coefs, nullptr, d_V_scratch_idx, -1, homogeneous_bc);
            d_hier_v_bdry_fill->setHomogeneousBc(homogeneous_bc);
            d_hier_v_bdry_fill->fillData(spec.eval_time);
            d_bc_helper->enforceDivergenceFreeConditionAtBoundary(
                d_V_scratch_idx, d_coarsest_ln, d_finest_ln, StaggeredStokesPhysicalBoundaryHelper::ALL_BDRY);
            enforceDivergenceFreeConditionAtCoarseFineInterface(d_V_scratch_idx);
            StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_u_bc_coefs, nullptr);
            d_hier_v_bdry_fill->resetTransactionComponents(d_v_transaction_comps);

            int local_nonfinite_div = 0;
            double local_max_div_dt = 0.0;
            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
                    const double* const dx = patch_geom->getDx();
                    const Box<NDIM>& patch_box = patch->getBox();
                    const IntVector<NDIM>& patch_lower = patch_box.lower();
                    const IntVector<NDIM>& patch_upper = patch_box.upper();

                    Pointer<SideData<NDIM, double>> V_data = patch->getPatchData(d_V_scratch_idx);
                    Pointer<SideData<NDIM, double>> E_data = patch->getPatchData(d_E_scratch_idx);

                    const IntVector<NDIM> ghosts = IntVector<NDIM>(1);
                    std::array<Box<NDIM>, NDIM> side_boxes;
                    std::array<Pointer<FaceData<NDIM, double>>, NDIM> V_adv_data;
                    std::array<Pointer<FaceData<NDIM, double>>, NDIM> R_one_data;
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                        V_adv_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                        R_one_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                        R_one_data[axis]->fillAll(1.0);
                    }
                    computeAdvectionVelocity(V_adv_data, V_data, patch_lower, patch_upper, side_boxes);
                    Pointer<SideData<NDIM, double>> R_zero_data =
                        new SideData<NDIM, double>(patch_box, 1, IntVector<NDIM>(0));
                    R_zero_data->fillAll(0.0);
                    computeErrorOfMassConservationEquation(
                        E_data, R_zero_data, R_zero_data, V_adv_data, R_one_data, side_boxes, dt, dx);
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (Box<NDIM>::Iterator it(side_boxes[axis]); it; it++)
                        {
                            const SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);
                            const double e_val = (*E_data)(si);
                            if (!std::isfinite(e_val))
                            {
                                local_nonfinite_div = 1;
                                continue;
                            }
                            local_max_div_dt = std::max(local_max_div_dt, dt * std::abs(e_val));
                        }
                    }
                }
            }
            const int global_nonfinite_div = IBTK_MPI::maxReduction(local_nonfinite_div);
            const double stage_max_div_dt = IBTK_MPI::maxReduction(local_max_div_dt);
            if (std::isfinite(stage_max_div_dt))
            {
                d_max_observed_stage_div_dt = std::max(d_max_observed_stage_div_dt, stage_max_div_dt);
            }
            if (global_nonfinite_div != 0 || !std::isfinite(stage_max_div_dt) ||
                stage_max_div_dt > d_const_density_max_stage_div_dt)
            {
                can_fast_path = false;
                break;
            }
        }

        if (!can_fast_path)
        {
            ++d_const_density_fallback_count;
        }
        else
        {
            ++d_const_density_fastpath_hits;

            // Stage 2 (midpoint) velocity is already filled in d_V_scratch_idx.
            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
                    const double* const dx = patch_geom->getDx();
                    const Box<NDIM>& patch_box = patch->getBox();
                    const IntVector<NDIM>& patch_lower = patch_box.lower();
                    const IntVector<NDIM>& patch_upper = patch_box.upper();

                    Pointer<SideData<NDIM, double>> N_data = patch->getPatchData(d_N_idx);
                    Pointer<SideData<NDIM, double>> V_data = patch->getPatchData(d_V_scratch_idx);
                    Pointer<SideData<NDIM, double>> V_cur_data = patch->getPatchData(d_V_current_idx);
                    Pointer<SideData<NDIM, double>> E_data = patch->getPatchData(d_E_scratch_idx);

                    const IntVector<NDIM> ghosts = IntVector<NDIM>(1);
                    std::array<Box<NDIM>, NDIM> side_boxes;
                    std::array<Pointer<FaceData<NDIM, double>>, NDIM> V_adv_data;
                    std::array<Pointer<FaceData<NDIM, double>>, NDIM> V_half_data;
                    std::array<Pointer<FaceData<NDIM, double>>, NDIM> R_half_data;
                    std::array<Pointer<FaceData<NDIM, double>>, NDIM> P_half_data;
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                        V_adv_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                        V_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                        R_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                        P_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                        R_half_data[axis]->fillAll(rho_const);
                    }

                    computeAdvectionVelocity(V_adv_data, V_data, patch_lower, patch_upper, side_boxes);
                    interpolateSideQuantity(V_half_data,
                                            V_adv_data,
                                            V_data,
                                            patch_lower,
                                            patch_upper,
                                            side_boxes,
                                            d_velocity_convective_limiter);

                    IBAMR_TIMER_START(t_apply_convective_operator);
                    computeConvectiveDerivative(
                        N_data, P_half_data, V_adv_data, R_half_data, V_half_data, side_boxes, dx, patch);
                    IBAMR_TIMER_STOP(t_apply_convective_operator);

                    // Cancel discrete mass-conservation error consistency term with uniform rho_new = rho_cur =
                    // rho_const.
                    Pointer<SideData<NDIM, double>> R_const_data =
                        new SideData<NDIM, double>(patch_box, 1, IntVector<NDIM>(0));
                    R_const_data->fillAll(rho_const);
                    computeErrorOfMassConservationEquation(
                        E_data, R_const_data, R_const_data, V_adv_data, R_half_data, side_boxes, dt, dx);
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (Box<NDIM>::Iterator it(side_boxes[axis]); it; it++)
                        {
                            SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);
                            (*N_data)(si) -= (*V_cur_data)(si) * (*E_data)(si);
                        }
                    }
                }
            }

            d_hier_sc_data_ops->setToScalar(d_rho_new_idx, rho_const, /*interior_only*/ false);

            d_N_idx = -1;
            d_rho_current_idx = -1;
            d_V_old_idx = -1;
            d_V_current_idx = -1;
            d_V_new_idx = -1;
            d_cycle_num = -1;
            d_dt_prev = -1.0;

            IBAMR_TIMER_STOP(t_integrate);
            return;
        }
    }

    // Fill ghost cells for current density
    std::vector<InterpolationTransactionComponent> rho_transaction_comps({ { d_rho_scratch_idx,
                                                                             d_rho_current_idx,
                                                                             "CONSERVATIVE_LINEAR_REFINE",
                                                                             false,
                                                                             "CONSERVATIVE_COARSEN",
                                                                             d_density_bdry_extrap_type,
                                                                             false,
                                                                             d_rho_bc_coefs } });
    d_hier_rho_bdry_fill->resetTransactionComponents(rho_transaction_comps);
    d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_rho_bdry_fill->fillData(d_current_time);
    d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

    // Fill ghost cells for the velocity used to compute the density update
    // Note, enforce divergence free condition on all physical boundaries to
    // ensure boundedness of density update
    d_hier_sc_data_ops->copyData(d_V_composite_idx,
                                 d_V_current_idx,
                                 /*interior_only*/ true);
    std::vector<InterpolationTransactionComponent> v_transaction_comps({ { d_V_scratch_idx,
                                                                           d_V_composite_idx,
                                                                           "CONSERVATIVE_LINEAR_REFINE",
                                                                           false,
                                                                           "CONSERVATIVE_COARSEN",
                                                                           d_velocity_bdry_extrap_type,
                                                                           false,
                                                                           d_u_bc_coefs } });
    d_hier_v_bdry_fill->resetTransactionComponents(v_transaction_comps);
    StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
        d_u_bc_coefs, nullptr, d_V_scratch_idx, -1, homogeneous_bc);
    d_hier_v_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_v_bdry_fill->fillData(d_current_time);
    d_bc_helper->enforceDivergenceFreeConditionAtBoundary(
        d_V_scratch_idx, d_coarsest_ln, d_finest_ln, StaggeredStokesPhysicalBoundaryHelper::ALL_BDRY);
    enforceDivergenceFreeConditionAtCoarseFineInterface(d_V_scratch_idx);
    StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_u_bc_coefs, nullptr);
    d_hier_v_bdry_fill->resetTransactionComponents(d_v_transaction_comps);

    // Compute the old mass
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    const double old_mass = d_hier_sc_data_ops->integral(d_rho_current_idx, wgt_sc_idx);
    if (d_enable_logging)
    {
        plog << "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate(): "
                "old mass in the domain = "
             << old_mass << "\n";
    }

    // Compute the convective derivative.
    for (int step = 0; step < d_num_steps; ++step)
    {
        double eval_time = std::numeric_limits<double>::quiet_NaN();
        double w0 = std::numeric_limits<double>::quiet_NaN();
        double w1 = std::numeric_limits<double>::quiet_NaN();
        double w2 = std::numeric_limits<double>::quiet_NaN();
        const double omega = dt / d_dt_prev;
        const double sum_dt = dt + d_dt_prev;

        switch (step)
        {
        case 0:
            eval_time = d_current_time;
            break;
        case 1:
            eval_time = d_current_time + dt;
            if (d_cycle_num > 0)
            {
                w0 = 0.0, w1 = 0.0, w2 = 1.0;
            }
            else
            {
                w0 = -1.0 * omega, w1 = 1.0 + omega, w2 = 0.0;
            }
            break;
        case 2:
            eval_time = d_current_time + dt / 2.0;
            if (d_cycle_num > 0)
            {
                w0 = -0.25 * dt * dt / (d_dt_prev * sum_dt);
                w1 = 0.25 * (2.0 + omega);
                w2 = 0.25 * (dt + 2.0 * d_dt_prev) / sum_dt;
            }
            else
            {
                w0 = -0.5 * omega, w1 = 1.0 + 0.5 * omega, w2 = 0.0;
            }
            break;
        default:
            TBOX_ERROR("This statement should not be reached");
        }
        // Fill ghost cells for new density and velocity, if needed
        if (step > 0)
        {
            std::vector<InterpolationTransactionComponent> update_transaction_comps({ { d_rho_scratch_idx,
                                                                                        d_rho_new_idx,
                                                                                        "CONSERVATIVE_LINEAR_REFINE",
                                                                                        false,
                                                                                        "CONSERVATIVE_COARSEN",
                                                                                        d_density_bdry_extrap_type,
                                                                                        false,
                                                                                        d_rho_bc_coefs } });
            d_hier_rho_bdry_fill->resetTransactionComponents(update_transaction_comps);
            d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
            d_hier_rho_bdry_fill->fillData(eval_time);
            d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

            // Compute an approximation to velocity at eval_time Note, enforce
            // divergence free condition on all physical boundaries to ensure
            // boundedness of density update
            d_hier_sc_data_ops->linearSum(
                d_V_composite_idx, w0, d_V_old_idx, w1, d_V_current_idx, /*interior_only*/ true);
            d_hier_sc_data_ops->axpy(d_V_composite_idx, w2, d_V_new_idx, d_V_composite_idx, /*interior_only*/ true);
            std::vector<InterpolationTransactionComponent> v_update_transaction_comps({ { d_V_scratch_idx,
                                                                                          d_V_composite_idx,
                                                                                          "CONSERVATIVE_LINEAR_REFINE",
                                                                                          false,
                                                                                          "CONSERVATIVE_COARSEN",
                                                                                          d_velocity_bdry_extrap_type,
                                                                                          false,
                                                                                          d_u_bc_coefs } });
            d_hier_v_bdry_fill->resetTransactionComponents(v_update_transaction_comps);
            StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
                d_u_bc_coefs, nullptr, d_V_scratch_idx, -1, homogeneous_bc);
            d_hier_v_bdry_fill->setHomogeneousBc(homogeneous_bc);
            d_hier_v_bdry_fill->fillData(eval_time);
            d_bc_helper->enforceDivergenceFreeConditionAtBoundary(
                d_V_scratch_idx, d_coarsest_ln, d_finest_ln, StaggeredStokesPhysicalBoundaryHelper::ALL_BDRY);
            enforceDivergenceFreeConditionAtCoarseFineInterface(d_V_scratch_idx);
            StaggeredStokesPhysicalBoundaryHelper::resetBcCoefObjects(d_u_bc_coefs, nullptr);
            d_hier_v_bdry_fill->resetTransactionComponents(d_v_transaction_comps);
        }

        // Compute the source term
        if (d_S_fcn)
        {
            d_S_fcn->setDataOnPatchHierarchy(d_S_scratch_idx, d_S_var, d_hierarchy, eval_time);
        }
        else
        {
            d_hier_sc_data_ops->setToScalar(d_S_scratch_idx, 0.0);
        }

        for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = d_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());

                const Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
                const double* const dx = patch_geom->getDx();

                const Box<NDIM>& patch_box = patch->getBox();
                const IntVector<NDIM>& patch_lower = patch_box.lower();
                const IntVector<NDIM>& patch_upper = patch_box.upper();

                Pointer<SideData<NDIM, double>> N_data = patch->getPatchData(d_N_idx);
                Pointer<SideData<NDIM, double>> V_data = patch->getPatchData(d_V_scratch_idx);
                Pointer<SideData<NDIM, double>> R_cur_data = patch->getPatchData(d_rho_current_idx);
                Pointer<SideData<NDIM, double>> R_pre_data = patch->getPatchData(d_rho_scratch_idx);
                Pointer<SideData<NDIM, double>> R_new_data = patch->getPatchData(d_rho_new_idx);
                Pointer<SideData<NDIM, double>> R_src_data = patch->getPatchData(d_S_scratch_idx);
                Pointer<SideData<NDIM, double>> E_data = patch->getPatchData(d_E_scratch_idx);

                // Define variables that live on the "faces" of control
                // volumes centered about side-centered staggered velocity
                // components
                const IntVector<NDIM> ghosts = IntVector<NDIM>(1);
                std::array<Box<NDIM>, NDIM> side_boxes;
                std::array<Pointer<FaceData<NDIM, double>>, NDIM> V_adv_data;
                std::array<Pointer<FaceData<NDIM, double>>, NDIM> V_half_data;
                std::array<Pointer<FaceData<NDIM, double>>, NDIM> R_half_data;
                std::array<Pointer<FaceData<NDIM, double>>, NDIM> P_half_data;
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    side_boxes[axis] = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                    V_adv_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                    V_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                    R_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                    P_half_data[axis] = new FaceData<NDIM, double>(side_boxes[axis], 1, ghosts);
                }
                // Interpolate velocity components onto "faces" using simple averages.
                computeAdvectionVelocity(V_adv_data, V_data, patch_lower, patch_upper, side_boxes);

                // Upwind side-centered densities onto faces.
                interpolateSideQuantity(R_half_data,
                                        V_adv_data,
                                        R_pre_data,
                                        patch_lower,
                                        patch_upper,
                                        side_boxes,
                                        d_density_convective_limiter);

                // Compute the convective derivative with the penultimate density and
                // velocity, if necessary
                if ((d_density_time_stepping_type == SSPRK2 && step == 1) ||
                    (d_density_time_stepping_type == SSPRK3 && step == 2))
                {
                    interpolateSideQuantity(V_half_data,
                                            V_adv_data,
                                            V_data,
                                            patch_lower,
                                            patch_upper,
                                            side_boxes,
                                            d_velocity_convective_limiter);

                    IBAMR_TIMER_START(t_apply_convective_operator);

                    computeConvectiveDerivative(
                        N_data, P_half_data, V_adv_data, R_half_data, V_half_data, side_boxes, dx, patch);

                    IBAMR_TIMER_STOP(t_apply_convective_operator);
                }

                // Compute the updated density
                double a0, a1, a2;
                switch (step)
                {
                case 0:
                    a0 = 0.5;
                    a1 = 0.5;
                    a2 = 1.0;
                    break;
                case 1:
                    if (d_density_time_stepping_type == SSPRK2)
                    {
                        a0 = 0.5;
                        a1 = 0.5;
                        a2 = 0.5;
                        break;
                    }
                    else if (d_density_time_stepping_type == SSPRK3)
                    {
                        a0 = 0.75;
                        a1 = 0.25;
                        a2 = 0.25;
                        break;
                    }
                    else
                    {
                        TBOX_ERROR("This statement should not be reached");
                        break;
                    }
                case 2:
                    a0 = 1.0 / 3.0;
                    a1 = 2.0 / 3.0;
                    a2 = 2.0 / 3.0;
                    break;
                default:
                    TBOX_ERROR("This statement should not be reached");
                }
                computeDensityUpdate(R_new_data,
                                     a0,
                                     R_cur_data,
                                     a1,
                                     R_pre_data,
                                     a2,
                                     V_adv_data,
                                     R_half_data,
                                     R_src_data,
                                     side_boxes,
                                     dt,
                                     dx);

                Pointer<SideData<NDIM, double>> V_cur_data = patch->getPatchData(d_V_current_idx);
                if ((d_density_time_stepping_type == SSPRK2 && step == 1) ||
                    (d_density_time_stepping_type == SSPRK3 && step == 2))
                {
                    computeErrorOfMassConservationEquation(
                        E_data, R_new_data, R_cur_data, V_adv_data, R_half_data, side_boxes, dt, dx);

                    // subtract Error*U from the convective operator.
                    for (unsigned int axis = 0; axis < NDIM; ++axis)
                    {
                        for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                        {
                            SideIndex<NDIM> si(it(), axis, SideIndex<NDIM>::Lower);

                            (*N_data)(si) -= (*V_cur_data)(si) * (*E_data)(si);
                        }
                    }
                }
            }
        }
    }

    // Refill boundary values of newest density
    const double new_time = d_current_time + dt;
    std::vector<InterpolationTransactionComponent> new_transaction_comps({ { d_rho_scratch_idx,
                                                                             d_rho_new_idx,
                                                                             "CONSERVATIVE_LINEAR_REFINE",
                                                                             false,
                                                                             "CONSERVATIVE_COARSEN",
                                                                             d_density_bdry_extrap_type,
                                                                             false,
                                                                             d_rho_bc_coefs } });
    d_hier_rho_bdry_fill->resetTransactionComponents(new_transaction_comps);
    d_hier_rho_bdry_fill->setHomogeneousBc(homogeneous_bc);
    d_hier_rho_bdry_fill->fillData(new_time);
    d_hier_rho_bdry_fill->resetTransactionComponents(d_rho_transaction_comps);

    d_hier_sc_data_ops->copyData(d_rho_new_idx,
                                 d_rho_scratch_idx,
                                 /*interior_only*/ true);

    // Compute the new mass
    const double new_mass = d_hier_sc_data_ops->integral(d_rho_new_idx, wgt_sc_idx);
    if (d_enable_logging)
    {
        plog << "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate(): "
                "new mass in the domain = "
             << new_mass << "\n";
        plog << "INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::integrate(): "
                "change in mass = "
             << new_mass - old_mass << "\n";
    }

    // Reset select options
    d_N_idx = -1;
    d_rho_current_idx = -1;
    d_V_old_idx = -1;
    d_V_current_idx = -1;
    d_V_new_idx = -1;
    d_cycle_num = -1;
    d_dt_prev = -1.0;

    IBAMR_TIMER_STOP(t_integrate);
    return;
} // integrate

long
INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::getConstantDensityFastPathHits() const
{
    return d_const_density_fastpath_hits;
} // getConstantDensityFastPathHits

long
INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::getConstantDensityFallbackCount() const
{
    return d_const_density_fallback_count;
} // getConstantDensityFallbackCount

double
INSVCStaggeredConservativeMassMomentumSSPRKIntegrator::getMaxObservedStageDivDt() const
{
    return d_max_observed_stage_div_dt;
} // getMaxObservedStageDivDt

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
