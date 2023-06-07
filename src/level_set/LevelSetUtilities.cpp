// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2020 by the IBAMR developers
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
#include "ibamr/LSInitStrategy.h"
#include "ibamr/LevelSetUtilities.h"

#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IBTK_MPI.h"
#include "ibtk/ibtk_utilities.h"

#include "tbox/Database.h"
#include "tbox/RestartManager.h"

#include "ibamr/app_namespaces.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{

/////////////////////////////// PUBLIC ///////////////////////////////////////

LevelSetUtilities::LevelSetMassLossFixer::LevelSetMassLossFixer(std::string object_name,
                                                                Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator,
                                                                Pointer<CellVariable<NDIM, double> > ls_var,
                                                                double ncells,
                                                                bool register_for_restart)
    : LevelSetContainer(adv_diff_integrator, ls_var, ncells),
      d_object_name(std::move(object_name)),
      d_registered_for_restart(register_for_restart)
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    bool is_from_restart = RestartManager::getManager()->isFromRestart();
    if (is_from_restart)
    {
        getFromRestart();
    }
    return;
} // LevelSetMassLossFixer

LevelSetUtilities::LevelSetMassLossFixer::~LevelSetMassLossFixer()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
        d_registered_for_restart = false;
    }
    return;
} // ~LevelSetMassLossFixer

void
LevelSetUtilities::LevelSetMassLossFixer::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_vol_init", d_vol_init);
} // putToDatabase

void
LevelSetUtilities::LevelSetMassLossFixer::setInitialVolume(double v0)
{
    if (RestartManager::getManager()->isFromRestart()) return;

    d_vol_init = v0;
    return;
} // setInitialVolume

void
LevelSetUtilities::SetLSProperties::setLSData(int ls_idx,
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

/////////////////////////////// STATIC ///////////////////////////////////////

void
LevelSetUtilities::fixLevelSetMassLoss(double /*current_time*/, double new_time, int /*cycle_num*/, void* ctx)
{
    LevelSetMassLossFixer* mass_fixer = static_cast<LevelSetMassLossFixer*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(mass_fixer);
#endif

    const double v0 = mass_fixer->getInitialVolume();
    const double ncells = mass_fixer->getInterfaceHalfWidth();
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = mass_fixer->getAdvDiffHierarchyIntegrator();

    Pointer<PatchHierarchy<NDIM> > patch_hier = adv_diff_integrator->getPatchHierarchy();
    Pointer<HierarchyMathOps> hier_math_ops = adv_diff_integrator->getHierarchyMathOps();

    // NOTE: In practice the level set mass loss would be fixed after advection. Hence the application time
    // would be the new time and the variable context would be the new context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_idx =
        var_db->mapVariableAndContextToIndex(mass_fixer->getLevelSetVariable(), adv_diff_integrator->getNewContext());
    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();

    const int hier_finest_ln = patch_hier->getFinestLevelNumber();
    double vol_integral = 0.0;
    double vol_interface_integral = 0.0;
    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(ls_idx);
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

                vol_integral += (1.0 - h_phi) * dv;
                vol_interface_integral += h_prime * dv;
            }
        }
    }
    const double nmr = IBTK_MPI::sumReduction(vol_integral) - v0;
    const double dnr = IBTK_MPI::sumReduction(vol_interface_integral);

    // Compute the spatially constant Lagrange multiplier and adjust the level set function
    const double q = nmr / dnr;

    HierarchyCellDataOpsReal<NDIM, double> hier_cc_ops(patch_hier, 0, hier_finest_ln);
    hier_cc_ops.addScalar(ls_idx, ls_idx, q);

    // For logging purposes.
    mass_fixer->setLagrangeMultiplier(q);
    mass_fixer->setTime(new_time);

    return;
} // fixLevelSetMassLoss

std::pair<double, double>
LevelSetUtilities::computeIntegralHeavisideFcns(double /*current_time*/,
                                                double /*new_time*/,
                                                int /*cycle_num*/,
                                                void* ctx)
{
    LevelSetContainer* lsc = static_cast<LevelSetContainer*>(ctx);
#if !defined(NDEBUG)
    TBOX_ASSERT(lsc);
#endif
    const double ncells = lsc->getInterfaceHalfWidth();

    Pointer<AdvDiffHierarchyIntegrator> adv_diff_integrator = lsc->getAdvDiffHierarchyIntegrator();
    Pointer<PatchHierarchy<NDIM> > patch_hier = adv_diff_integrator->getPatchHierarchy();
    Pointer<HierarchyMathOps> hier_math_ops = adv_diff_integrator->getHierarchyMathOps();

    // NOTE: In practice the level set mass is computed after integrating the hierarchy. Hence the application time
    // would be the new time and the variable context would be the current context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int ls_idx =
        var_db->mapVariableAndContextToIndex(lsc->getLevelSetVariable(), adv_diff_integrator->getCurrentContext());
    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();

    const int hier_finest_ln = patch_hier->getFinestLevelNumber();
    double vol_integral = 0.0;
    double vol_integral_complement = 0.0;
    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();

            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(ls_idx);
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

                // smoothed Heaviside function
                const double h_phi = IBTK::smooth_heaviside(phi, alpha);

                vol_integral += h_phi * dv;
                vol_integral_complement += (1.0 - h_phi) * dv;
            }
        }
    }
    const double H1 = IBTK_MPI::sumReduction(vol_integral);
    const double H2 = IBTK_MPI::sumReduction(vol_integral_complement);

    return std::make_pair(H1, H2);

} // computeIntegralHeavisideFcns

void
LevelSetUtilities::setLSDataPatchHierarchy(int ls_idx,
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
LevelSetUtilities::LevelSetMassLossFixer::getFromRestart()
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

    d_vol_init = db->getDouble("d_vol_init");
    return;
}

} // namespace IBAMR