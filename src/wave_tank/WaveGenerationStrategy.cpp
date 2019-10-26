// Filename: WaveGenerationStrategy.h
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

#include "ibamr/WaveGenerationStrategy.h"
#include "CartesianGridGeometry.h"
#include "ibamr/app_namespaces.h"

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
    auto stokes_wave_generator = static_cast<StokesWaveGenerator*>(ctx);
    const double x_zone_start = stokes_wave_generator->d_x_zone_start;
    const double x_zone_end = stokes_wave_generator->d_x_zone_end;
    const double alpha = stokes_wave_generator->d_alpha;
    const double sign_gas = stokes_wave_generator->d_sign_gas_phase;
    const double depth = stokes_wave_generator->getWaterDepth();

    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = stokes_wave_generator->d_ins_hier_integrator->getPatchHierarchy();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    int u_new_idx =
        var_db->mapVariableAndContextToIndex(stokes_wave_generator->d_ins_hier_integrator->getVelocityVariable(),
                                             stokes_wave_generator->d_ins_hier_integrator->getNewContext());

    Pointer<CellVariable<NDIM, double> > phi_cc_var = stokes_wave_generator->d_phi_var;
    if (!phi_cc_var)
    {
        TBOX_ERROR("callFifthOrderStokesWaveRelaxationCallbackFunction(): Level set variable must be cell centered!");
    }
    int phi_new_idx = var_db->mapVariableAndContextToIndex(
        phi_cc_var, stokes_wave_generator->d_adv_diff_hier_integrator->getNewContext());
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
            auto beta =
                stokes_wave_generator->d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

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

/////////////////////////////// PUBLIC //////////////////////////////////////

StokesWaveGenerator::StokesWaveGenerator(const std::string& object_name, Pointer<Database> input_db)
{
    d_object_name = object_name;
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
    TBOX_ASSERT(input_db);
#endif

    // Get wave parameters.
    getFromInput(input_db);

    return;
} // StokesWaveGenerator

FifthOrderStokesWaveGenerator::FifthOrderStokesWaveGenerator(const std::string& object_name, Pointer<Database> input_db)
    : StokesWaveGenerator(object_name, input_db)
{
    initStokesCoefficients();
    return;
} // FifthOrderStokesWaveGenerator

double
FifthOrderStokesWaveGenerator::getSurfaceElevation(double x, double time) const
{
    // This computes eta for both shallow and deep water.
    const double& k = d_wave_number;
    const double ka = k * d_amplitude;
    const double omega = sqrt(d_gravity * k) * (d_C[0] + ka * ka * d_C[2] + ka * ka * ka * ka * d_C[4]);
    const double phase = k * x - omega * time;

    return (d_eta[0] * cos(phase) + d_eta[1] * cos(2 * phase) + d_eta[2] * cos(3 * phase) + d_eta[3] * cos(4 * phase) +
            d_eta[4] * cos(5 * phase)) /
           k;

} // getSurfaceElevation

double
FifthOrderStokesWaveGenerator::getVelocity(double x, double z_plus_d, double time, int comp_idx) const
{
    const double& k = d_wave_number;
    const double ka = k * d_amplitude;
    const double omega = sqrt(d_gravity * k) * (d_C[0] + ka * ka * d_C[2] + ka * ka * ka * ka * d_C[4]);
    const double phase = k * x - omega * time;

    if (!d_deep_water_limit)
    {
        double kzd = k * z_plus_d;
        if (comp_idx == 0)
        {
            return d_p[0] * cosh(kzd) * cos(phase) * k + d_p[1] * cosh(2 * kzd) * cos(2 * phase) * k * 2 +
                   d_p[2] * cosh(3 * kzd) * cos(3 * phase) * k * 3 + d_p[3] * cosh(4 * kzd) * cos(4 * phase) * k * 4 +
                   d_p[4] * cosh(5 * kzd) * cos(5 * phase) * k * 5;
        }

        if (comp_idx == 1)
        {
#if (NDIM == 2)
            return d_p[0] * sinh(kzd) * sin(phase) * k + d_p[1] * sinh(2 * kzd) * sin(2 * phase) * k * 2 +
                   d_p[2] * sinh(3 * kzd) * sin(3 * phase) * k * 3 + d_p[3] * sinh(4 * kzd) * sin(4 * phase) * k * 4 +
                   d_p[4] * sinh(5 * kzd) * sin(5 * phase) * k * 5;
#elif (NDIM == 3)
            return 0;
#endif
        }
#if (NDIM == 3)
        if (comp_idx == 2)
        {
            return d_p[0] * sinh(kzd) * sin(phase) * k + d_p[1] * sinh(2 * kzd) * sin(2 * phase) * k * 2 +
                   d_p[2] * sinh(3 * kzd) * sin(3 * phase) * k * 3 + d_p[3] * sinh(4 * kzd) * sin(4 * phase) * k * 4 +
                   d_p[4] * sinh(5 * kzd) * sin(5 * phase) * k * 5;
        }
#endif
    }
    else
    {
        double kz = k * (z_plus_d - d_depth);
        if (comp_idx == 0)
        {
            return d_p[0] * exp(kz) * cos(phase) * k + d_p[1] * exp(2 * kz) * cos(2 * phase) * k * 2 +
                   d_p[2] * exp(3 * kz) * cos(3 * phase) * k * 3;
        }

        if (comp_idx == 1)
        {
#if (NDIM == 2)
            return d_p[0] * exp(kz) * sin(phase) * k + d_p[1] * exp(2 * kz) * sin(2 * phase) * k * 2 +
                   d_p[2] * exp(3 * kz) * sin(3 * phase) * k * 3;
#elif (NDIM == 3)

            return 0;
#endif
        }

#if (NDIM == 3)
        if (comp_idx == 2)
        {
            return d_p[0] * exp(kz) * sin(phase) * k + d_p[1] * exp(2 * kz) * sin(2 * phase) * k * 2 +
                   d_p[2] * exp(3 * kz) * sin(3 * phase) * k * 3;
        }
#endif
    }

    return std::numeric_limits<double>::signaling_NaN();
} // getVelocity

/////////////////////////////// PRIVATE //////////////////////////////////////
void
StokesWaveGenerator::getFromInput(Pointer<Database> input_db)
{
    Pointer<Database> wave_db = input_db;
    if (input_db->isDatabase("wave_parameters_db"))
    {
        wave_db = input_db->getDatabase("wave_parameters_db");
    }

    d_depth = wave_db->getDouble("depth");
    d_omega = wave_db->getDoubleWithDefault("omega", std::numeric_limits<double>::quiet_NaN());
    d_gravity = wave_db->getDouble("gravitational_constant");
    d_wave_number = wave_db->getDouble("wave_number");
    d_amplitude = wave_db->getDouble("amplitude");
    d_num_interface_cells = wave_db->getDouble("num_interface_cells");
    d_deep_water_limit = wave_db->getBoolWithDefault("deep_water_limit", d_deep_water_limit);

    return;
} // getFromInput

void
FifthOrderStokesWaveGenerator::initStokesCoefficients()
{
    const double kd = d_wave_number * d_depth;
    const double ka = d_wave_number * d_amplitude;
    const double ka2 = ka * ka;
    const double ka3 = ka2 * ka;
    const double ka4 = ka3 * ka;
    const double ka5 = ka4 * ka;

    if (d_deep_water_limit)
    {
        // Terms to calculate velocity.
        d_A[1][1] = 1;
        d_A[3][1] = -0.5;
        d_A[4][2] = 0.5;
        d_A[5][1] = -37.0 / 24;
        d_A[5][3] = 1.0 / 12;

        // Terms to calculate eta.
        d_B[2][2] = 0.5;
        d_B[3][1] = -3 / 8.0;
        d_B[4][2] = 1 / 3.0;
        d_B[4][4] = 1 / 3.0;
        d_B[5][3] = 99 / 128.0;
        d_B[5][5] = 125 / 384.0;

        // Wave Dispersion
        d_C[0] = 1;
        d_C[2] = 0.5;
        d_C[4] = 1 / 8.0;

        // Velocity coefficients.
        const double alpha = sqrt(d_gravity / pow(d_wave_number, 3));
        d_p[0] = alpha * (ka * d_A[1][1] + ka3 * d_A[3][1] + ka5 * d_A[5][1]);
        d_p[1] = alpha * (ka4 * d_A[4][2]);
        d_p[2] = alpha * (ka5 * d_A[5][3]);
    }
    else
    {
        double S = 1 / cosh(2 * kd);

        // Terms to calculate velocity.
        d_A[1][1] = 1 / sinh(kd);
        d_A[2][2] = 3 * S * S / (2 * (1 - S) * (1 - S));
        d_A[3][1] = (-4 - 20 * S + 10 * S * S - 13 * S * S * S) / (8 * sinh(kd) * pow((1 - S), 3));
        d_A[3][3] = (-2 * S * S + 11 * S * S * S) / (8 * sinh(kd) * pow((1 - S), 3));
        d_A[4][2] = (12 * S - 14 * S * S - 264 * S * S * S - 45 * S * S * S * S - 13 * S * S * S * S * S) /
                    (24 * pow((1 - S), 5));
        d_A[4][4] = (10 * S * S * S - 174 * S * S * S * S + 291 * pow(S, 5) + 278 * pow(S, 6)) /
                    (48 * (3 + 2 * S) * pow((1 - S), 5));
        d_A[5][1] = (-1184 + 32 * S + 13232 * S * S + 21712 * S * S * S + 20940 * S * S * S * S + 12554 * pow(S, 5) -
                     500 * pow(S, 6) - 3341 * pow(S, 7) - 670 * pow(S, 8)) /
                    (64 * sinh(kd) * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));
        d_A[5][3] = (4 * S + 105 * S * S + 198 * S * S * S - 1376 * S * S * S * S - 1302 * pow(S, 5) - 117 * pow(S, 6) +
                     58 * pow(S, 7)) /
                    (32 * sinh(kd) * (3 + 2 * S) * pow((1 - S), 6));
        d_A[5][5] = (-6 * S * S * S + 272 * S * S * S * S - 1552 * pow(S, 5) + 852 * pow(S, 6) + 2029 * pow(S, 7) +
                     430 * pow(S, 8)) /
                    (64 * sinh(kd) * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));

        // Terms to calculate eta.
        d_B[2][2] = 1 / tanh(kd) * (1 + 2 * S) / (2 * (1 - S));
        d_B[3][1] = -3 * (1 + 3 * S + 3 * S * S + 2 * S * S * S) / (8 * pow((1 - S), 3));
        d_B[4][2] = 1 / tanh(kd) * (6 - 26 * S - 182 * S * S - 204 * S * S * S - 25 * pow(S, 4) + 26 * pow(S, 5)) /
                    (6 * (3 + 2 * S) * pow((1 - S), 4));
        d_B[4][4] = 1 / tanh(kd) * (24 + 92 * S + 122 * S * S + 66 * S * S * S + 67 * pow(S, 4) + 34 * pow(S, 5)) /
                    (24 * (3 + 2 * S) * pow((1 - S), 4));
        d_B[5][3] = 9 *
                    (132 + 17 * S - 2216 * S * S - 5897 * S * S * S - 6292 * pow(S, 4) - 2687 * pow(S, 5) +
                     194 * pow(S, 6) + 467 * pow(S, 7) + 82 * pow(S, 8)) /
                    (128 * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));
        d_B[5][5] = 5 *
                    (300 + 1579 * S + 3176 * S * S + 2949 * S * S * S + 1188 * pow(S, 4) + 675 * pow(S, 5) +
                     1326 * pow(S, 6) + 827 * pow(S, 7) + 130 * pow(S, 8)) /
                    (384 * (3 + 2 * S) * (4 + S) * pow((1 - S), 6));

        // Wave dispersion
        d_C[0] = sqrt(tanh(kd));
        d_C[2] = d_C[0] * (2 + 7 * S * S) / (4 * pow((1 - S), 2));
        d_C[4] = d_C[0] * (4 + 32 * S - 116 * S * S - 400 * S * S * S - 71 * pow(S, 4) + 146 * pow(S, 5)) /
                 (32 * pow((1 - S), 5));

        // Velocity coefficients.
        const double alpha = d_C[0] * sqrt(d_gravity / pow(d_wave_number, 3));
        d_p[0] = alpha * (ka * d_A[1][1] + ka3 * d_A[3][1] + ka5 * d_A[5][1]);
        d_p[1] = alpha * (ka2 * d_A[2][2] + ka4 * d_A[4][2]);
        d_p[2] = alpha * (ka3 * d_A[3][3] + ka5 * d_A[5][3]);
        d_p[3] = alpha * (ka4 * d_A[4][4]);
        d_p[4] = alpha * (ka5 * d_A[5][5]);
    }

    // Component wave eta.
    d_eta[0] = ka + ka3 * d_B[3][1] - ka5 * (d_B[5][3] + d_B[5][5]);
    d_eta[1] = ka2 * d_B[2][2] + ka4 * d_B[4][2];
    d_eta[2] = -ka3 * d_B[3][1] + ka5 * d_B[5][3];
    d_eta[3] = ka4 * d_B[4][4];
    d_eta[4] = ka5 * d_B[5][5];

    return;

} // initStokesCoefficients

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
