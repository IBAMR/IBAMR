// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/StokesWaveGeneratorStrategy.h"
#include "ibamr/WaveGenerationFunctions.h"
#include "ibamr/WaveUtilities.h"

#include "Box.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "CellVariable.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "SideData.h"
#include "SideGeometry.h"
#include "SideIndex.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "tbox/Pointer.h"
#include "tbox/Utilities.h"

#include <cmath>
#include <string>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

namespace WaveGenerationFunctions
{
void
callStokesWaveRelaxationCallbackFunction(double /*current_time*/,
                                         double new_time,
                                         bool /*skip_synchronize_new_state_data*/,
                                         int /*num_cycles*/,
                                         void* ctx)
{
    auto stokes_wave_generator = static_cast<StokesWaveGeneratorStrategy*>(ctx);
    const double x_zone_start = stokes_wave_generator->d_wave_gen_data.d_x_zone_start;
    const double x_zone_end = stokes_wave_generator->d_wave_gen_data.d_x_zone_end;
    const double alpha = stokes_wave_generator->d_wave_gen_data.d_alpha;
    const double sign_gas = stokes_wave_generator->d_wave_gen_data.d_sign_gas_phase;
    const double depth = stokes_wave_generator->getWaterDepth();

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy =
        stokes_wave_generator->d_wave_gen_data.d_ins_hier_integrator->getPatchHierarchy();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int u_new_idx = var_db->mapVariableAndContextToIndex(
        stokes_wave_generator->d_wave_gen_data.d_ins_hier_integrator->getVelocityVariable(),
        stokes_wave_generator->d_wave_gen_data.d_ins_hier_integrator->getNewContext());

    Pointer<CellVariable<NDIM, double> > phi_cc_var = stokes_wave_generator->d_wave_gen_data.d_phi_var;
    if (!phi_cc_var)
    {
        TBOX_ERROR(
            "callFifthOrderStokesWaveRelaxationCallbackFunction(): Level "
            "set variable must be cell centered!");
    }
    int phi_new_idx = var_db->mapVariableAndContextToIndex(
        phi_cc_var, stokes_wave_generator->d_wave_gen_data.d_adv_diff_hier_integrator->getNewContext());
    static const int dir = NDIM == 2 ? 1 : 2;

    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_x_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();

            Pointer<SideData<NDIM, double> > u_data = patch->getPatchData(u_new_idx);

            // Compute a representative grid spacing
            double vol_cell = 1.0;
            for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
            auto beta = stokes_wave_generator->d_wave_gen_data.d_num_interface_cells *
                        std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch_box, axis)); it; it++)
                {
                    hier::Index<NDIM> i = it();
                    SideIndex<NDIM> i_side(i, axis, SideIndex<NDIM>::Lower);
                    double x_posn = patch_x_lower[0] + patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)));
                    const double shift_x = (axis == 0 ? 0.0 : 0.5);
                    x_posn += patch_dx[0] * shift_x;

                    double dir_posn =
                        patch_x_lower[dir] + patch_dx[dir] * (static_cast<double>(i(dir) - patch_lower(dir)));
                    const double shift_dir = (axis == dir ? 0.0 : 0.5);
                    dir_posn += patch_dx[dir] * shift_dir;

                    // Compute a numerical heaviside from the analytical wave elevation
                    const double z_plus_d = dir_posn;
                    const double eta = stokes_wave_generator->getSurfaceElevation(x_posn, new_time);
                    const double phi = -eta + (z_plus_d - depth);
                    double h_phi;
                    if (phi < -beta)
                        h_phi = 1.0;
                    else if (std::abs(phi) <= beta)
                        h_phi = 1.0 - (0.5 + 0.5 * phi / beta + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / beta));
                    else
                        h_phi = 0.0;

                    if (x_posn >= x_zone_start && x_posn <= x_zone_end)
                    {
                        const double xtilde = (x_posn - x_zone_start) / (x_zone_end - x_zone_start);
                        const double gamma = 1.0 - std::expm1(std::pow(xtilde, alpha)) / std::expm1(1.0);

                        if (axis == 0 || axis == (NDIM - 1))
                        {
                            const double target = h_phi * stokes_wave_generator->getVelocity(x_posn,
                                                                                             z_plus_d,
                                                                                             new_time,
                                                                                             /*comp*/ axis);
                            (*u_data)(i_side, 0) = (1.0 - gamma) * (*u_data)(i_side, 0) + gamma * target;
                        }
                    }
                }
            }
        }
    }

    // Modify the level set in the generation zone.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* const patch_dx = patch_geom->getDx();
            const double* const patch_x_lower = patch_geom->getXLower();
            const Box<NDIM>& patch_box = patch->getBox();
            const IntVector<NDIM>& patch_lower = patch_box.lower();

            Pointer<CellData<NDIM, double> > phi_data = patch->getPatchData(phi_new_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                hier::Index<NDIM> i = it();

                const double x_posn =
                    patch_x_lower[0] + patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + 0.5);

                const double dir_posn =
                    patch_x_lower[dir] + patch_dx[dir] * (static_cast<double>(i(dir) - patch_lower(dir)) + 0.5);

                if (x_posn >= x_zone_start && x_posn <= x_zone_end)
                {
                    const double eta = stokes_wave_generator->getSurfaceElevation(x_posn, new_time);
                    const double xtilde = (x_posn - x_zone_start) / (x_zone_end - x_zone_start);
                    const double gamma = 1.0 - std::expm1(std::pow(xtilde, alpha)) / std::expm1(1.0);

                    const double target = sign_gas * (-eta + dir_posn - depth);

                    (*phi_data)(i, 0) = (1.0 - gamma) * (*phi_data)(i, 0) + gamma * target;
                }
            }
        }
    }

    return;
} // callStokesWaveRelaxationCallbackFunction

} // namespace WaveGenerationFunctions

} // namespace IBAMR
