// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2026 by the IBAMR developers
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

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/LevelSetUtilities.h>

#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/RestartManager.h>

#include <algorithm>
#include <cmath>
#include <iterator>

#include <ibamr/app_namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
std::vector<double>
compute_heaviside_integrals(Pointer<HierarchyMathOps> hier_math_ops,
                            const int phi_idx,
                            const double ncells,
                            const double shift = 0.0)
{
    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
    Pointer<PatchHierarchy<NDIM>> patch_hier = hier_math_ops->getPatchHierarchy();

    const int hier_finest_ln = patch_hier->getFinestLevelNumber();
    double vol_phase1 = 0.0;
    double vol_phase2 = 0.0;
    double integral_delta = 0.0;
    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double>> phi_data = patch->getPatchData(phi_idx);
            Pointer<CellData<NDIM, double>> wgt_data = patch->getPatchData(wgt_cc_idx);

            // Get grid spacing information
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            double cell_size = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                cell_size *= patch_dx[d];
            }
            cell_size = std::pow(cell_size, 1.0 / static_cast<double>(NDIM));
            const double alpha = ncells * cell_size;

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double phi = (*phi_data)(ci) + shift;
                const double dv = (*wgt_data)(ci);

                // smoothed delta and Heaviside functions
                const double h_phi = IBTK::smooth_heaviside(phi, alpha);
                const double h_prime = IBTK::smooth_delta(phi, alpha);

                vol_phase1 += IBTK::smooth_heaviside(-phi, alpha) * dv;
                vol_phase2 += h_phi * dv;
                integral_delta += h_prime * dv;
            }
        }
    }

    std::vector<double> integrals{ vol_phase1, vol_phase2, integral_delta };
    IBTK_MPI::sumReduction(&integrals[0], integrals.size());

    return integrals;
} // compute_heaviside_integrals

std::vector<double>
compute_heaviside_integrals(Pointer<HierarchyMathOps> hier_math_ops,
                            const int phi_idx,
                            const int psi_idx,
                            const double ncells,
                            const double shift = 0.0)
{
    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
    Pointer<PatchHierarchy<NDIM>> patch_hier = hier_math_ops->getPatchHierarchy();

    const int hier_finest_ln = patch_hier->getFinestLevelNumber();
    double vol_phase1 = 0.0;
    double vol_phase2 = 0.0;
    double vol_phase3 = 0.0;
    double integral_delta = 0.0;
    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double>> phi_data = patch->getPatchData(phi_idx);
            Pointer<CellData<NDIM, double>> psi_data = patch->getPatchData(psi_idx);
            Pointer<CellData<NDIM, double>> wgt_data = patch->getPatchData(wgt_cc_idx);

            // Get grid spacing information
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            double cell_size = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                cell_size *= patch_dx[d];
            }
            cell_size = std::pow(cell_size, 1.0 / static_cast<double>(NDIM));
            const double alpha = ncells * cell_size;

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double phi = (*phi_data)(ci) + shift;
                const double psi = (*psi_data)(ci);
                const double dv = (*wgt_data)(ci);

                // smoothed delta and Heaviside functions
                const double h_phi = IBTK::smooth_heaviside(phi, alpha);
                const double h_phi_prime = IBTK::smooth_delta(phi, alpha);
                const double h_psi = IBTK::smooth_heaviside(psi, alpha);

                vol_phase1 += IBTK::smooth_heaviside(-phi, alpha) * h_psi * dv;
                vol_phase2 += h_phi * h_psi * dv;
                vol_phase3 += IBTK::smooth_heaviside(-psi, alpha) * dv;
                integral_delta += h_phi_prime * h_psi * dv;
            }
        }
    }

    std::vector<double> integrals{ vol_phase1, vol_phase2, vol_phase3, integral_delta };
    IBTK_MPI::sumReduction(&integrals[0], integrals.size());

    return integrals;
} // compute_heaviside_integrals

} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

namespace LevelSetUtilities
{
void
tagLSCells(Pointer<BasePatchHierarchy<NDIM>> hierarchy,
           const int level_number,
           const double /*error_data_time*/,
           const int tag_index,
           const bool initial_time,
           const bool /*uses_richardson_extrapolation_too*/,
           void* ctx)
{
    if (initial_time || level_number == hierarchy->getFinestLevelNumber()) return;

    TagLSRefinementCells* ls_tagger = static_cast<TagLSRefinementCells*>(ctx);

#if !defined(NDEBUG)
    TBOX_ASSERT(ls_tagger);
    TBOX_ASSERT(hierarchy);
    TBOX_ASSERT(level_number >= 0);
    TBOX_ASSERT(hierarchy->getPatchLevel(level_number));
#endif

    const LevelSetContainer& ls_container = ls_tagger->getLevelSetContainer();

    // Get the level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_idx = var_db->mapVariableAndContextToIndex(
        ls_container.getLevelSetVariable(), ls_container.getAdvDiffHierarchyIntegrator()->getCurrentContext());

    // Get the tagging criterion
    const double& tag_min_val = ls_tagger->getTagMinValue();
    const double& tag_max_val = ls_tagger->getTagMaxValue();

    // Tag cells based on the value of the level set variable
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, int>> tags_data = patch->getPatchData(tag_index);
        Pointer<CellData<NDIM, double>> ls_data = patch->getPatchData(ls_idx);

        for (CellIterator<NDIM> ic(patch_box); ic; ic++)
        {
            const hier::Index<NDIM>& i = ic();
            const double dist = (*ls_data)(i);

            if (dist >= tag_min_val && dist <= tag_max_val)
            {
                (*tags_data)(i) = 1;
            }
        }
    }

    return;
} // tagLSCells

LevelSetMassLossFixer::LevelSetMassLossFixer(std::string object_name,
                                             Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                                             std::vector<Pointer<CellVariable<NDIM, double>>> ls_vars,
                                             Pointer<Database> input_db,
                                             bool register_for_restart)
    : d_object_name(std::move(object_name)),
      d_ls_container(adv_diff_integrator, ls_vars),
      d_registered_for_restart(register_for_restart)
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    if (input_db)
    {
        getFromInput(input_db);
    }

    if (RestartManager::getManager()->isFromRestart())
    {
        getFromRestart();
    }
    const double ncells = d_ls_container.getInterfaceHalfWidth();
    if (d_interval <= 0 || d_max_its <= 0 || !std::isfinite(d_rel_tol) || d_rel_tol < 0.0 || d_rel_tol >= 1.0 ||
        (d_abs_tol && (!std::isfinite(*d_abs_tol) || *d_abs_tol < 0.0)) || !std::isfinite(ncells) || ncells <= 0.0)
    {
        TBOX_ERROR(d_object_name
                   << "::LevelSetMassLossFixer(): invalid correction controls\n"
                   << "  correction_interval and max_its must be positive; half_width must be finite and positive;\n"
                   << "  rel_tol must be finite in [0,1), and abs_tol must be finite and nonnegative.\n");
    }
    return;
} // LevelSetMassLossFixer

LevelSetMassLossFixer::~LevelSetMassLossFixer()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~LevelSetMassLossFixer

void
LevelSetMassLossFixer::putToDatabase(Pointer<Database> db)
{
    const LevelSetContainer& ls_container = getLevelSetContainer();
    db->putDouble("vol_init", d_vol_init);
    db->putDouble("ncells", ls_container.getInterfaceHalfWidth());
    db->putDouble("vol_target", d_vol_target);
} // putToDatabase

void
LevelSetMassLossFixer::setInitialVolume(double v0)
{
    if (RestartManager::getManager()->isFromRestart()) return;

    d_vol_init = v0;

    // Set the default target volume as the initial
    // phase volume.
    d_vol_target = v0;
    return;
} // setInitialVolume

void
LevelSetMassLossFixer::correctVolume(const double new_time, const bool three_phase)
{
    Pointer<AdvDiffHierarchyIntegrator> integrator = d_ls_container.getAdvDiffHierarchyIntegrator();
    const int step = integrator->getIntegratorStep();
    if (IBTK_MPI::minReduction(d_interval) != IBTK_MPI::maxReduction(d_interval))
    {
        TBOX_ERROR(d_object_name << "::correctVolume(): inconsistent correction intervals across ranks\n");
    }
    if (step % d_interval != 0)
    {
        return;
    }

    double capacity = std::numeric_limits<double>::quiet_NaN();
    double residual = std::numeric_limits<double>::quiet_NaN();
    int trials = 0;
    const auto fail = [&](const char* reason)
    {
        TBOX_ERROR(d_object_name << "::correctVolume(): " << reason << '\n'
                                 << "  phase = " << (three_phase ? "three-phase liquid" : "two-phase gas")
                                 << ", target = " << d_vol_target << ", capacity = " << capacity
                                 << ", residual = " << residual << ", trials = " << trials << '\n');
    };
    const double ncells = d_ls_container.getInterfaceHalfWidth();
    if (!std::isfinite(ncells) || ncells <= 0.0)
    {
        fail("half_width must be finite and positive");
    }
    if (!std::isfinite(d_vol_target))
    {
        fail("target volume must be finite");
    }
    double controls_min[] = { ncells,
                              d_vol_target,
                              d_rel_tol,
                              d_abs_tol.value_or(-1.0),
                              static_cast<double>(d_max_its),
                              static_cast<double>(three_phase) };
    double controls_max[6];
    std::copy(std::begin(controls_min), std::end(controls_min), controls_max);
    IBTK_MPI::minReduction(controls_min, 6);
    IBTK_MPI::maxReduction(controls_max, 6);
    if (!std::equal(std::begin(controls_min), std::end(controls_min), controls_max))
    {
        fail("inconsistent correction controls across ranks");
    }

    Pointer<PatchHierarchy<NDIM>> hierarchy = integrator->getPatchHierarchy();
    Pointer<HierarchyMathOps> ops = integrator->getHierarchyMathOps();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellVariable<NDIM, double>> variable = d_ls_container.getLevelSetVariable();
    const int phi_idx = var_db->mapVariableAndContextToIndex(variable, integrator->getNewContext());
    const int psi_idx = three_phase ? var_db->mapVariableAndContextToIndex(d_ls_container.getLevelSetVariable(1),
                                                                           integrator->getNewContext()) :
                                      -1;
    const int weight_idx = ops->getCellWeightPatchDescriptorIndex();
    const int finest_ln = hierarchy->getFinestLevelNumber();
    double lower = std::numeric_limits<double>::infinity();
    double upper = -std::numeric_limits<double>::infinity();
    double local_capacity = 0.0;
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(phi_idx) || (three_phase && !level->checkAllocated(psi_idx)))
        {
            fail("NEW level-set data are not allocated");
        }
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> phi = patch->getPatchData(phi_idx);
            Pointer<CellData<NDIM, double>> psi = three_phase ? patch->getPatchData(psi_idx) : nullptr;
            Pointer<CellData<NDIM, double>> weight = patch->getPatchData(weight_idx);
            Pointer<CartesianPatchGeometry<NDIM>> geom = patch->getPatchGeometry();
            double cell_volume = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                cell_volume *= geom->getDx()[d];
            }
            const double alpha = ncells * std::pow(cell_volume, 1.0 / NDIM);
            if (!std::isfinite(alpha) || alpha <= 0.0)
            {
                fail("smoothing half-width is not representable");
            }
            for (Box<NDIM>::Iterator i(patch->getBox()); i; i++)
            {
                const double value = (*phi)(i());
                const double dv = (*weight)(i());
                if (!std::isfinite(value) || !std::isfinite(dv) || dv < 0.0 ||
                    (three_phase && !std::isfinite((*psi)(i()))))
                {
                    fail("nonfinite level-set data or invalid composite weight");
                }
                const double fluid_fraction = three_phase ? IBTK::smooth_heaviside((*psi)(i()), alpha) : 1.0;
                local_capacity += fluid_fraction * dv;
                if (fluid_fraction * dv > 0.0)
                {
                    lower = std::min(lower, -alpha - value);
                    upper = std::max(upper, alpha - value);
                }
            }
        }
    }
    capacity = IBTK_MPI::sumReduction(local_capacity);
    if (!std::isfinite(capacity) || capacity < 0.0)
    {
        fail("invalid fluid capacity");
    }
    const double abs_tol = d_abs_tol.value_or(64.0 * std::numeric_limits<double>::epsilon() * capacity);
    const double tolerance = abs_tol + d_rel_tol * std::abs(d_vol_target);
    if (!std::isfinite(tolerance))
    {
        fail("volume tolerance is not representable");
    }
    if (d_vol_target < -tolerance || (d_vol_target > capacity && d_vol_target - capacity > tolerance))
    {
        fail("target volume is outside the fluid capacity");
    }
    const auto evaluate = [&](const int idx, const double q)
    {
        return three_phase ? compute_heaviside_integrals(ops, idx, psi_idx, ncells, q) :
                             compute_heaviside_integrals(ops, idx, ncells, q);
    };
    const int phase = three_phase ? 1 : 0;
    const double orientation = three_phase ? 1.0 : -1.0;
    std::vector<double> integrals = evaluate(phi_idx, 0.0);
    residual = integrals[phase] - d_vol_target;
    if (!std::isfinite(residual))
    {
        fail("nonfinite phase volume");
    }
    double q = 0.0;
    if (std::abs(residual) > tolerance)
    {
        lower = std::nextafter(IBTK_MPI::minReduction(lower), -std::numeric_limits<double>::infinity());
        upper = std::nextafter(IBTK_MPI::maxReduction(upper), std::numeric_limits<double>::infinity());
        if (!std::isfinite(lower) || !std::isfinite(upper) || !(lower < upper))
        {
            fail("cannot form a finite correction bracket");
        }
        const double internal_target = std::max(0.0, std::min(capacity, d_vol_target));
        const double v_lower = evaluate(phi_idx, lower)[phase];
        const double v_upper = evaluate(phi_idx, upper)[phase];
        if (!std::isfinite(v_lower) || !std::isfinite(v_upper) || orientation * (v_upper - v_lower) < 0.0 ||
            orientation * (v_lower - internal_target) > tolerance ||
            orientation * (v_upper - internal_target) < -tolerance)
        {
            fail("invalid correction bracket volumes");
        }
        if (std::abs(v_lower - d_vol_target) <= tolerance)
        {
            q = lower;
            residual = v_lower - d_vol_target;
        }
        else if (std::abs(v_upper - d_vol_target) <= tolerance)
        {
            q = upper;
            residual = v_upper - d_vol_target;
        }
        while (std::abs(residual) > tolerance && trials < d_max_its)
        {
            if (q > lower && q < upper)
            {
                if (orientation * (integrals[phase] - internal_target) < 0.0)
                {
                    lower = q;
                }
                else
                {
                    upper = q;
                }
            }
            double candidate = 0.5 * lower + 0.5 * upper;
            const double derivative = orientation * integrals.back();
            if (std::isfinite(derivative) && derivative != 0.0)
            {
                const double newton = q - (integrals[phase] - internal_target) / derivative;
                if (std::isfinite(newton) && newton > lower && newton < upper && newton != q &&
                    std::abs(newton - q) <= 0.5 * upper - 0.5 * lower)
                {
                    candidate = newton;
                }
            }
            if (!(candidate > lower && candidate < upper) || candidate == q)
            {
                fail("correction stagnated before satisfying the volume tolerance");
            }
            q = candidate;
            integrals = evaluate(phi_idx, q);
            ++trials;
            residual = integrals[phase] - d_vol_target;
            if (!std::isfinite(residual))
            {
                fail("nonfinite correction trial volume");
            }
        }
        if (std::abs(residual) > tolerance)
        {
            fail("correction trial budget exhausted");
        }

        const int candidate_idx = var_db->registerClonedPatchDataIndex(variable, phi_idx);
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            level->allocatePatchData(candidate_idx, new_time);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> original = patch->getPatchData(phi_idx);
                Pointer<CellData<NDIM, double>> candidate = patch->getPatchData(candidate_idx);
                for (Box<NDIM>::Iterator i(patch->getBox()); i; i++)
                {
                    const double value = (*original)(i()) + q;
                    if (!std::isfinite(value))
                    {
                        fail("nonfinite candidate level-set value");
                    }
                    (*candidate)(i()) = value;
                }
            }
        }
        residual = evaluate(candidate_idx, 0.0)[phase] - d_vol_target;
        if (!std::isfinite(residual) || std::abs(residual) > tolerance)
        {
            fail("candidate field does not satisfy the volume tolerance");
        }
        HierarchyCellDataOpsReal<NDIM, double> data_ops(hierarchy, 0, finest_ln);
        data_ops.copyData(phi_idx, candidate_idx, true);
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(candidate_idx);
        }
        var_db->removePatchDataIndex(candidate_idx);
    }
    d_q = q;
    d_time = new_time;
    if (d_enable_logging)
    {
        plog << d_object_name << "::correctVolume(): phase = " << (three_phase ? "three-phase liquid" : "two-phase gas")
             << ", target = " << d_vol_target << ", capacity = " << capacity << ", residual = " << residual
             << ", tolerance = " << tolerance << ", trials = " << trials << ", shift = " << q << '\n';
    }
}

void
SetLSProperties::setLSData(int ls_idx,
                           SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                           const int integrator_step,
                           const double current_time,
                           const bool initial_time,
                           const bool regrid_time)
{
    // If at the regrid time, force reinitialization
    d_ls_ops->setReinitializeLSData(regrid_time);
    d_ls_ops->initializeLSData(ls_idx, hier_math_ops, integrator_step, current_time, initial_time);

    return;
} // setLSData

void
fixMassLoss2PhaseFlows(double /*current_time*/,
                       double new_time,
                       bool /*skip_synchronize_new_state_data*/,
                       int /*num_cycles*/,
                       void* ctx)
{
    LevelSetMassLossFixer* mass_fixer = static_cast<LevelSetMassLossFixer*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(mass_fixer);
#endif
    mass_fixer->correctVolume(new_time, false);
} // fixMassLoss2PhaseFlows

void
fixMassLoss3PhaseFlows(double /*current_time*/,
                       double new_time,
                       bool /*skip_synchronize_new_state_data*/,
                       int /*num_cycles*/,
                       void* ctx)
{
    LevelSetMassLossFixer* mass_fixer = static_cast<LevelSetMassLossFixer*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(mass_fixer);
#endif
    mass_fixer->correctVolume(new_time, true);
} // fixMassLoss3PhaseFlows

std::vector<double>
computeHeavisideIntegrals2PhaseFlows(const LevelSetContainer& lsc)
{
    const double ncells = lsc.getInterfaceHalfWidth();
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = lsc.getAdvDiffHierarchyIntegrator();
    Pointer<PatchHierarchy<NDIM>> patch_hier = adv_diff_integrator->getPatchHierarchy();
    Pointer<HierarchyMathOps> hier_math_ops = adv_diff_integrator->getHierarchyMathOps();

    // NOTE: In practice the level set mass is computed after integrating the hierarchy. Hence the application time
    // would be the new time and the variable context would be the current context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_idx =
        var_db->mapVariableAndContextToIndex(lsc.getLevelSetVariable(), adv_diff_integrator->getCurrentContext());

    std::vector<double> integrals = compute_heaviside_integrals(hier_math_ops, ls_idx, ncells);

    return integrals;

} // computeHeavisideIntegrals2PhaseFlows

std::vector<double>
computeHeavisideIntegrals3PhaseFlows(const LevelSetContainer& lsc)
{
    const double ncells = lsc.getInterfaceHalfWidth();
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = lsc.getAdvDiffHierarchyIntegrator();
    Pointer<PatchHierarchy<NDIM>> patch_hier = adv_diff_integrator->getPatchHierarchy();
    Pointer<HierarchyMathOps> hier_math_ops = adv_diff_integrator->getHierarchyMathOps();

    // NOTE: In practice the level set mass is computed after integrating the hierarchy. Hence the application time
    // would be the new time and the variable context would be the current context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int fluid_ls_idx =
        var_db->mapVariableAndContextToIndex(lsc.getLevelSetVariable(0), adv_diff_integrator->getCurrentContext());
    const int solid_ls_idx =
        var_db->mapVariableAndContextToIndex(lsc.getLevelSetVariable(1), adv_diff_integrator->getCurrentContext());

    std::vector<double> integrals = compute_heaviside_integrals(hier_math_ops, fluid_ls_idx, solid_ls_idx, ncells);

    return integrals;
} // computeHeavisideIntegrals3PhaseFlows

void
setLSDataPatchHierarchy(int ls_idx,
                        Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                        const int integrator_step,
                        const double current_time,
                        const bool initial_time,
                        const bool regrid_time,
                        void* ctx)
{
    // Set the density from the level set information
    SetLSProperties* ptr_SetLSProperties = static_cast<SetLSProperties*>(ctx);
    ptr_SetLSProperties->setLSData(ls_idx, hier_math_ops, integrator_step, current_time, initial_time, regrid_time);

    return;
} // setLSDataPatchHierarchy
////////////////////////////// PROTECTED ///////////////////////////////////////

void
LevelSetMassLossFixer::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }

    d_vol_init = db->getDouble("vol_init");
    d_vol_target = db->getDouble("vol_target");

    LevelSetContainer& ls_container = getLevelSetContainer();
    ls_container.setInterfaceHalfWidth(db->getDouble("ncells"));

    return;
} // getFromRestart

void
LevelSetMassLossFixer::getFromInput(Pointer<Database> input_db)
{
    d_enable_logging = input_db->getBoolWithDefault("enable_logging", false);
    d_interval = input_db->getIntegerWithDefault("correction_interval", 1);
    d_max_its = input_db->getIntegerWithDefault("max_its", 64);
    d_rel_tol = input_db->getDoubleWithDefault("rel_tol", 1e-12);
    if (input_db->keyExists("abs_tol"))
    {
        d_abs_tol = input_db->getDouble("abs_tol");
    }

    LevelSetContainer& ls_container = getLevelSetContainer();
    ls_container.setInterfaceHalfWidth(input_db->getDoubleWithDefault("half_width", 1.0));

    return;
} // getFromInput

} // namespace LevelSetUtilities

} // namespace IBAMR
