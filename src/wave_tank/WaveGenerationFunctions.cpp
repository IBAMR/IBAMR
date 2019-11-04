// Filename: WaveGenerationFunctions.h
// Created on 15 Oct 2019 by Amneet Bhalla
//
// Copyright (c) 2002-2019, Amneet Bhalla
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/INSVCStaggeredHierarchyIntegrator.h"
#include "ibamr/StokesWaveGeneratorStrategy.h"
#include "ibamr/WaveGenerationFunctions.h"
#include "ibamr/app_namespaces.h"

#include "CartesianGridGeometry.h"

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
        TBOX_ERROR("callFifthOrderStokesWaveRelaxationCallbackFunction(): Level set variable must be cell centered!");
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
                    Index<NDIM> i = it();
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
                        const double gamma = 1.0 - (exp(std::pow(xtilde, alpha)) - 1.0) / (exp(1.0) - 1.0);

                        if (axis == 0 || axis == (NDIM - 1))
                        {
                            const double target =
                                h_phi * stokes_wave_generator->getVelocity(x_posn, z_plus_d, new_time, /*comp*/ axis);
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
                Index<NDIM> i = it();

                const double x_posn =
                    patch_x_lower[0] + patch_dx[0] * (static_cast<double>(i(0) - patch_lower(0)) + 0.5);

                const double dir_posn =
                    patch_x_lower[dir] + patch_dx[dir] * (static_cast<double>(i(dir) - patch_lower(dir)) + 0.5);

                if (x_posn >= x_zone_start && x_posn <= x_zone_end)
                {
                    const double eta = stokes_wave_generator->getSurfaceElevation(x_posn, new_time);
                    const double xtilde = (x_posn - x_zone_start) / (x_zone_end - x_zone_start);
                    const double gamma = 1.0 - (exp(std::pow(xtilde, alpha)) - 1.0) / (exp(1.0) - 1.0);

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
