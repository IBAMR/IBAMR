// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
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

#include "ibamr/BrinkmanPenalizationMethod.h"

#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/ibtk_utilities.h"

#include <string>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC //////////////////////////////////////
BrinkmanPenalizationMethod::BrinkmanPenalizationMethod(std::string object_name,
                                                       Pointer<HierarchyIntegrator> time_integrator,
                                                       Pointer<Database> input_db,
                                                       bool register_for_restart)
    : BrinkmanPenalizationStrategy(std::move(object_name), register_for_restart), d_time_integrator(time_integrator)
{
    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    return;
} // BrinkmanPenalizationMethod

void
BrinkmanPenalizationMethod::registerSolidLevelSet(int ls_current_idx,
                                                  int ls_new_idx,
                                                  int ls_scratch_idx,
                                                  SAMRAI::solv::RobinBcCoefStrategy<NDIM>* ls_bc_coef)
{
    d_ls_current_idx = ls_current_idx;
    d_ls_new_idx = ls_new_idx;
    d_ls_scratch_idx = ls_scratch_idx;
    d_ls_bc_coef = ls_bc_coef;
    return;
} // registerSolidLevelSet

void
BrinkmanPenalizationMethod::computeBrinkmanVelocity(int u_idx, double time, int /*cycle_num*/)
{
    // Ghost fill the level set values.
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_time_integrator->getPatchHierarchy();
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent transaction_comp(d_ls_scratch_idx,
                                                       d_ls_new_idx,
                                                       "CONSERVATIVE_LINEAR_REFINE",
                                                       false,
                                                       "CONSERVATIVE_COARSEN",
                                                       "LINEAR",
                                                       false,
                                                       d_ls_bc_coef);

    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(transaction_comp, patch_hierarchy);
    hier_bdry_fill->fillData(time);

    // Set the rigid body velocity in u_idx
    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double h_min = *(std::min_element(patch_dx, patch_dx + NDIM));
        double vol_cell = 1.0;
        for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
        const double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

        Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(d_ls_scratch_idx);
        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                const double phi_lower = (*phi_data)(s_i.toCell(0));
                const double phi_upper = (*phi_data)(s_i.toCell(1));
                const double phi = 0.5 * (phi_lower + phi_upper);
                const double Hphi = IBTK::smooth_heaviside(phi, alpha);

                if (phi <= alpha)
                {
                    const double solid_velocity = 0.0;
                    (*u_data)(s_i) = solid_velocity * d_penalty_factor * (1.0 - Hphi);
                }
            }
        }
    }

    return;
} // computeBrinkmanVelocity

void
BrinkmanPenalizationMethod::demarcateBrinkmanZone(int u_idx, double time, int /*cycle_num*/)
{
    // Ghost fill the level set values.
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = d_time_integrator->getPatchHierarchy();
    typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent transaction_comp(d_ls_scratch_idx,
                                                       d_ls_new_idx,
                                                       "CONSERVATIVE_LINEAR_REFINE",
                                                       false,
                                                       "CONSERVATIVE_COARSEN",
                                                       "LINEAR",
                                                       false,
                                                       d_ls_bc_coef);

    Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
    hier_bdry_fill->initializeOperatorState(transaction_comp, patch_hierarchy);
    hier_bdry_fill->fillData(time);

    int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
    for (PatchLevel<NDIM>::Iterator p(finest_level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = finest_level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
        const double* const patch_dx = patch_geom->getDx();
        const double h_min = *(std::min_element(patch_dx, patch_dx + NDIM));
        double vol_cell = 1.0;
        for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
        const double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

        Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(d_ls_scratch_idx);
        Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_idx);

        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
            {
                SideIndex<NDIM> s_i(it(), axis, SideIndex<NDIM>::Lower);

                const double phi_lower = (*phi_data)(s_i.toCell(0));
                const double phi_upper = (*phi_data)(s_i.toCell(1));
                const double phi = 0.5 * (phi_lower + phi_upper);
                const double Hphi = IBTK::smooth_heaviside(phi, alpha);

                if (phi <= alpha)
                {
                    (*u_data)(s_i) = d_penalty_factor * (1.0 - Hphi);
                }
            }
        }
    }

    return;
} // demarcateBrinkmanZone

void
BrinkmanPenalizationMethod::putToDatabase(Pointer<Database> db)
{
    db->putDouble("penalty_factor", d_penalty_factor);
    db->putDouble("num_interface_cells", d_num_interface_cells);

    return;
} // putToDatabase

/////////////////////////////// PRIVATE //////////////////////////////////////

void
BrinkmanPenalizationMethod::getFromInput(Pointer<Database> input_db, bool is_from_restart)
{
    if (!is_from_restart)
    {
        if (input_db->keyExists("penalty_factor"))
        {
            d_penalty_factor = input_db->getDouble("penalty_factor");
        }

        if (input_db->keyExists("num_interface_cells"))
        {
            d_num_interface_cells = input_db->getDouble("num_interface_cells");
        }
    }
    return;
} // getFromInput

void
BrinkmanPenalizationMethod::getFromRestart()
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
    d_penalty_factor = db->getDouble("penalty_factor");
    d_num_interface_cells = db->getDouble("num_interface_cells");

    return;
} // getFromRestart

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
