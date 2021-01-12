// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
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

#include "ibamr/StokesFifthOrderWaveBcCoef.h"

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "BoxArray.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Database.h"

#include <cmath>
#include <limits>
#include <utility>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int EXTENSIONS_FILLABLE = 128;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

StokesFifthOrderWaveBcCoef::StokesFifthOrderWaveBcCoef(std::string object_name,
                                                       const int comp_idx,
                                                       Pointer<Database> input_db,
                                                       Pointer<CartesianGridGeometry<NDIM> > grid_geom)
    : d_object_name(std::move(object_name)),
      d_comp_idx(comp_idx),
      d_muparser_bcs(d_object_name + "::muParser", input_db, grid_geom),
      d_grid_geom(grid_geom)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
    TBOX_ASSERT(input_db);
#endif
    // Get wave parameters.
    getFromInput(input_db);

    // Initialize Stokes coefficients.
    initStokesCoefficients();

    return;
} // StokesFifthOrderWaveBcCoef

StokesFifthOrderWaveBcCoef::~StokesFifthOrderWaveBcCoef()
{
    // intentionally left-blank.
    return;
} // ~StokesFifthOrderWaveBcCoef

void
StokesFifthOrderWaveBcCoef::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
                                       Pointer<ArrayData<NDIM, double> >& bcoef_data,
                                       Pointer<ArrayData<NDIM, double> >& gcoef_data,
                                       const Pointer<Variable<NDIM> >& variable,
                                       const Patch<NDIM>& patch,
                                       const BoundaryBox<NDIM>& bdry_box,
                                       double fill_time) const
{
    // Get pgeom info.
    const Box<NDIM>& patch_box = patch.getBox();
    const SAMRAI::hier::Index<NDIM>& patch_lower = patch_box.lower();
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch.getPatchGeometry();
    const double* const x_lower = pgeom->getXLower();
    const double* const dx = pgeom->getDx();

    // Compute a representative grid spacing
    double vol_cell = 1.0;
    for (int d = 0; d < NDIM; ++d) vol_cell *= dx[d];
    auto alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / static_cast<double>(NDIM));

    // Get physical domain box extents on patch box level.
    const SAMRAI::hier::Box<NDIM> domain_box =
        SAMRAI::hier::Box<NDIM>::refine(d_grid_geom->getPhysicalDomain()[0], pgeom->getRatio());
    const int dir = (NDIM == 2) ? 1 : (NDIM == 3) ? 2 : -1;

    // We take location index = 0, i.e., x_lower to be the wave inlet.
    const unsigned int location_index = bdry_box.getLocationIndex();
    if (location_index != 0)
    {
        d_muparser_bcs.setBcCoefs(acoef_data, bcoef_data, gcoef_data, variable, patch, bdry_box, fill_time);
    }
    else
    {
        const unsigned int bdry_normal_axis = location_index / 2;
        const Box<NDIM>& bc_coef_box = (acoef_data ? acoef_data->getBox() :
                                        bcoef_data ? bcoef_data->getBox() :
                                        gcoef_data ? gcoef_data->getBox() :
                                                     Box<NDIM>());
#if !defined(NDEBUG)
        TBOX_ASSERT(!acoef_data || bc_coef_box == acoef_data->getBox());
        TBOX_ASSERT(!bcoef_data || bc_coef_box == bcoef_data->getBox());
        TBOX_ASSERT(!gcoef_data || bc_coef_box == gcoef_data->getBox());
#endif

        for (Box<NDIM>::Iterator b(bc_coef_box); b; b++)
        {
            const SAMRAI::hier::Index<NDIM>& i = b();
            if (acoef_data) (*acoef_data)(i, 0) = 1.0;
            if (bcoef_data) (*bcoef_data)(i, 0) = 0.0;

            IBTK::Point dof_posn;
            dof_posn.fill(std::numeric_limits<double>::signaling_NaN());
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                if (d != bdry_normal_axis)
                {
                    dof_posn[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)) + 0.5);
                }
                else
                {
                    dof_posn[d] = x_lower[d] + dx[d] * (static_cast<double>(i(d) - patch_lower(d)));
                }
            }

            // Compute a numerical heaviside at the boundary from the analytical wave
            // elevation
            const double z_plus_d = dof_posn[dir];
            const double eta = getSurfaceElevation(dof_posn[0], fill_time);
            const double phi = -eta + (z_plus_d - d_depth);
            double h_phi;
            if (phi < -alpha)
                h_phi = 1.0;
            else if (std::abs(phi) <= alpha)
                h_phi = 1.0 - (0.5 + 0.5 * phi / alpha + 1.0 / (2.0 * M_PI) * std::sin(M_PI * phi / alpha));
            else
                h_phi = 0.0;

            if (gcoef_data)
            {
                (*gcoef_data)(i, 0) = h_phi * getVelocity(dof_posn[0], z_plus_d, fill_time);
            }
        }
    }
    return;
} // setBcCoefs

IntVector<NDIM>
StokesFifthOrderWaveBcCoef::numberOfExtensionsFillable() const
{
    return EXTENSIONS_FILLABLE;
} // numberOfExtensionsFillable

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////
void
StokesFifthOrderWaveBcCoef::getFromInput(Pointer<Database> input_db)
{
    Pointer<Database> wave_db = input_db->getDatabase("wave_parameters_db");
#if !defined(NDEBUG)
    TBOX_ASSERT(input_db->isDatabase("wave_parameters_db"));
#endif

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
StokesFifthOrderWaveBcCoef::initStokesCoefficients()
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

double
StokesFifthOrderWaveBcCoef::getSurfaceElevation(double x, double time) const
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
StokesFifthOrderWaveBcCoef::getVelocity(double x, double z_plus_d, double time) const
{
    const double& k = d_wave_number;
    const double ka = k * d_amplitude;
    const double omega = sqrt(d_gravity * k) * (d_C[0] + ka * ka * d_C[2] + ka * ka * ka * ka * d_C[4]);
    const double phase = k * x - omega * time;

    if (!d_deep_water_limit)
    {
        double kzd = k * z_plus_d;
        if (d_comp_idx == 0)
        {
            return d_p[0] * cosh(kzd) * cos(phase) * k + d_p[1] * cosh(2 * kzd) * cos(2 * phase) * k * 2 +
                   d_p[2] * cosh(3 * kzd) * cos(3 * phase) * k * 3 + d_p[3] * cosh(4 * kzd) * cos(4 * phase) * k * 4 +
                   d_p[4] * cosh(5 * kzd) * cos(5 * phase) * k * 5;
        }

        if (d_comp_idx == 1)
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
        if (d_comp_idx == 2)
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
        if (d_comp_idx == 0)
        {
            return d_p[0] * exp(kz) * cos(phase) * k + d_p[1] * exp(2 * kz) * cos(2 * phase) * k * 2 +
                   d_p[2] * exp(3 * kz) * cos(3 * phase) * k * 3;
        }

        if (d_comp_idx == 1)
        {
#if (NDIM == 2)
            return d_p[0] * exp(kz) * sin(phase) * k + d_p[1] * exp(2 * kz) * sin(2 * phase) * k * 2 +
                   d_p[2] * exp(3 * kz) * sin(3 * phase) * k * 3;
#elif (NDIM == 3)

            return 0;
#endif
        }

#if (NDIM == 3)
        if (d_comp_idx == 2)
        {
            return d_p[0] * exp(kz) * sin(phase) * k + d_p[1] * exp(2 * kz) * sin(2 * phase) * k * 2 +
                   d_p[2] * exp(3 * kz) * sin(3 * phase) * k * 3;
        }
#endif
    }

    return std::numeric_limits<double>::signaling_NaN();
} // getVelocity

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR
