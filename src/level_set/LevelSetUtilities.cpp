// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2024 by the IBAMR developers
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

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/LevelSetUtilities.h"

#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/ibtk_utilities.h"

#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/RestartManager.h"

#include "ibamr/app_namespaces.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
std::vector<double>
compute_heaviside_integrals(Pointer<HierarchyMathOps> hier_math_ops, int phi_idx, double ncells)
{
    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
    Pointer<PatchHierarchy<NDIM> > patch_hier = hier_math_ops->getPatchHierarchy();

    const int hier_finest_ln = patch_hier->getFinestLevelNumber();
    double vol_phase1 = 0.0;
    double vol_phase2 = 0.0;
    double integral_delta = 0.0;
    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_idx);
            Pointer<CellData<NDIM, double> > wgt_data = patch->getPatchData(wgt_cc_idx);

            // Get grid spacing information
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            double cell_size = 1.0;
            for (int d = 0; d < NDIM; ++d) cell_size *= patch_dx[d];
            cell_size = std::pow(cell_size, 1.0 / static_cast<double>(NDIM));
            const double alpha = ncells * cell_size;

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double phi = (*phi_data)(ci);
                const double dv = (*wgt_data)(ci);

                // smoothed delta and Heaviside functions
                const double h_phi = IBTK::smooth_heaviside(phi, alpha);
                const double h_prime = IBTK::smooth_delta(phi, alpha);

                vol_phase1 += (1.0 - h_phi) * dv;
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
compute_heaviside_integrals(Pointer<HierarchyMathOps> hier_math_ops, int phi_idx, int psi_idx, double ncells)
{
    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
    Pointer<PatchHierarchy<NDIM> > patch_hier = hier_math_ops->getPatchHierarchy();

    const int hier_finest_ln = patch_hier->getFinestLevelNumber();
    double vol_phase1 = 0.0;
    double vol_phase2 = 0.0;
    double vol_phase3 = 0.0;
    double integral_delta = 0.0;
    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_idx);
            Pointer<CellData<NDIM, double> > psi_data = patch->getPatchData(psi_idx);
            Pointer<CellData<NDIM, double> > wgt_data = patch->getPatchData(wgt_cc_idx);

            // Get grid spacing information
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            double cell_size = 1.0;
            for (int d = 0; d < NDIM; ++d) cell_size *= patch_dx[d];
            cell_size = std::pow(cell_size, 1.0 / static_cast<double>(NDIM));
            const double alpha = ncells * cell_size;

            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());

                const double phi = (*phi_data)(ci);
                const double psi = (*psi_data)(ci);
                const double dv = (*wgt_data)(ci);

                // smoothed delta and Heaviside functions
                const double h_phi = IBTK::smooth_heaviside(phi, alpha);
                const double h_phi_prime = IBTK::smooth_delta(phi, alpha);
                const double h_psi = IBTK::smooth_heaviside(psi, alpha);

                vol_phase1 += (1.0 - h_phi) * h_psi * dv;
                vol_phase2 += h_phi * h_psi * dv;
                vol_phase3 += (1.0 - h_psi) * dv;
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
tagLSCells(Pointer<BasePatchHierarchy<NDIM> > hierarchy,
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
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, int> > tags_data = patch->getPatchData(tag_index);
        Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_idx);

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
                                             std::vector<Pointer<CellVariable<NDIM, double> > > ls_vars,
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

    bool is_from_restart = RestartManager::getManager()->isFromRestart();
    if (is_from_restart)
    {
        getFromRestart();
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
    const LevelSetContainer& ls_container = mass_fixer->getLevelSetContainer();
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = ls_container.getAdvDiffHierarchyIntegrator();
    const int integrator_step = adv_diff_integrator->getIntegratorStep();
    const int mass_correction_interval = mass_fixer->getCorrectionInterval();

    if (integrator_step % mass_correction_interval != 0) return;

    Pointer<PatchHierarchy<NDIM> > patch_hier = adv_diff_integrator->getPatchHierarchy();
    Pointer<HierarchyMathOps> hier_math_ops = adv_diff_integrator->getHierarchyMathOps();

    const int hier_finest_ln = patch_hier->getFinestLevelNumber();
    const double vol_target = mass_fixer->getTargetVolume();
    const double ncells = ls_container.getInterfaceHalfWidth();

    // NOTE: In practice the level set mass loss would be fixed during the postprocess integrate hierarchy stage.
    // Hence the application time would be the new time and the variable context would be the new context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_idx =
        var_db->mapVariableAndContextToIndex(ls_container.getLevelSetVariable(), adv_diff_integrator->getNewContext());

    // Carry out the Newton iterations
    double rel_error = 1.0e12;
    int current_iter = 0;

    const double min_rel_error = mass_fixer->getErrorRelTolerance();
    const int max_its = mass_fixer->getMaxIterations();

    double q = 0.0;
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_ops(patch_hier, 0, hier_finest_ln);
    while (rel_error > min_rel_error && current_iter < max_its)
    {
        std::vector<double> integrals = compute_heaviside_integrals(hier_math_ops, ls_idx, ncells);
        const double& vol_phase1 = integrals[0]; // Target the gas volume
        const double& integral_delta = integrals[2];

        rel_error = std::abs(vol_phase1 / vol_target - 1.0);

        if (mass_fixer->enableLogging())
        {
            plog << "fixMassLoss2PhaseFlows():: current iter = " << current_iter << " , rel error  = " << rel_error
                 << std::endl;
        }

        const double delta_q = (vol_phase1 - vol_target) / integral_delta;
        hier_cc_ops.addScalar(ls_idx, ls_idx, delta_q);

        q += delta_q;
        current_iter += 1;
    }

    // For logging purposes.
    mass_fixer->setLagrangeMultiplier(q);
    mass_fixer->setTime(new_time);

    return;
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
    const LevelSetContainer& ls_container = mass_fixer->getLevelSetContainer();
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = ls_container.getAdvDiffHierarchyIntegrator();
    const int integrator_step = adv_diff_integrator->getIntegratorStep();
    const int mass_correction_interval = mass_fixer->getCorrectionInterval();

    if (integrator_step % mass_correction_interval != 0) return;

    Pointer<PatchHierarchy<NDIM> > patch_hier = adv_diff_integrator->getPatchHierarchy();
    Pointer<HierarchyMathOps> hier_math_ops = adv_diff_integrator->getHierarchyMathOps();

    const int hier_finest_ln = patch_hier->getFinestLevelNumber();
    const double vol_target = mass_fixer->getTargetVolume();
    const double ncells = ls_container.getInterfaceHalfWidth();

    // NOTE: In practice the level set mass loss would be fixed during the postprocess integrate hierarchy stage.
    // Hence the application time would be the new time and the variable context would be the new context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int fluid_ls_idx =
        var_db->mapVariableAndContextToIndex(ls_container.getLevelSetVariable(0), adv_diff_integrator->getNewContext());
    const int solid_ls_idx =
        var_db->mapVariableAndContextToIndex(ls_container.getLevelSetVariable(1), adv_diff_integrator->getNewContext());

    // Carry out the Newton iterations
    double rel_error = 1.0e12;
    int current_iter = 0;

    const double min_rel_error = mass_fixer->getErrorRelTolerance();
    const int max_its = mass_fixer->getMaxIterations();

    double q = 0.0;
    HierarchyCellDataOpsReal<NDIM, double> hier_cc_ops(patch_hier, 0, hier_finest_ln);
    while (rel_error > min_rel_error && current_iter < max_its)
    {
        std::vector<double> integrals = compute_heaviside_integrals(hier_math_ops, fluid_ls_idx, solid_ls_idx, ncells);
        const double& vol_phase2 = integrals[1]; // Target the liquid volume
        const double& integral_delta = integrals[3];

        rel_error = std::abs(vol_phase2 / vol_target - 1.0);

        if (mass_fixer->enableLogging())
        {
            plog << "fixMassLoss3PhaseFlows():: current iter = " << current_iter << " , rel error  = " << rel_error
                 << std::endl;
        }

        const double delta_q = -(vol_phase2 - vol_target) / integral_delta;
        hier_cc_ops.addScalar(fluid_ls_idx, fluid_ls_idx, delta_q);

        q += delta_q;
        current_iter += 1;
    }

    // For logging purposes.
    mass_fixer->setLagrangeMultiplier(q);
    mass_fixer->setTime(new_time);

    return;
} // fixMassLoss3PhaseFlows

std::vector<double>
computeHeavisideIntegrals2PhaseFlows(const LevelSetContainer& lsc)
{
    const double ncells = lsc.getInterfaceHalfWidth();
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = lsc.getAdvDiffHierarchyIntegrator();
    Pointer<PatchHierarchy<NDIM> > patch_hier = adv_diff_integrator->getPatchHierarchy();
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
    Pointer<PatchHierarchy<NDIM> > patch_hier = adv_diff_integrator->getPatchHierarchy();
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
LevelSetMassLossFixer::getFromInput(Pointer<Database> input_db)
{
    d_enable_logging = input_db->getBoolWithDefault("enable_logging", false);
    d_interval = input_db->getIntegerWithDefault("correction_interval", 1);
    d_max_its = input_db->getIntegerWithDefault("max_its", 4);
    d_rel_tol = input_db->getDoubleWithDefault("rel_tol", 1e-12);

    LevelSetContainer& ls_container = getLevelSetContainer();
    ls_container.setInterfaceHalfWidth(input_db->getDoubleWithDefault("half_width", 1.0));

    return;
} // getFromInput

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

} // namespace LevelSetUtilities

} // namespace IBAMR
