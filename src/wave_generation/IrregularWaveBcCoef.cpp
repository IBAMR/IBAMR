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

#include "ibamr/IrregularWaveBcCoef.h"
#include "ibamr/RNG.h"

#include "ibtk/IBTK_MPI.h"

#include "ArrayData.h"
#include "BoundaryBox.h"
#include "Box.h"
#include "CartesianGridGeometry.h"
#include "CartesianPatchGeometry.h"
#include "Index.h"
#include "IntVector.h"
#include "Patch.h"
#include "tbox/Database.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <utility>

#include "ibamr/namespaces.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int EXTENSIONS_FILLABLE = 128;
static const unsigned SEED = 1234567;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IrregularWaveBcCoef::IrregularWaveBcCoef(std::string object_name,
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

    // Resize vectors according to number of wave components
    d_amplitude.resize(d_num_waves);
    d_wave_number.resize(d_num_waves);
    d_omega.resize(d_num_waves);
    d_phase.resize(d_num_waves);

    double delta_omega = std::abs(d_omega_end - d_omega_begin) / (d_num_waves - 1);
    double omega_s = 2 * M_PI / d_Ts;

    RNG::srandgen(SEED);

    // Calculating component waves
    {
        for (int i = 0; i < d_num_waves; i++)
        {
            double rn;
            RNG::genrand(&rn);
            d_phase[i] = 2 * M_PI * rn;

            d_omega[i] = d_omega_begin + i * delta_omega;

            // Using an approximate formula for the dispersion relationship, calculate
            // the wave number.
            // See Eqn. (5.4.22) in WAVES IN OCEANIC AND COASTAL WATERS by LEO H.
            // HOLTHUIJSEN.
            const double alpha = std::pow(d_omega[i], 2) * d_depth / d_gravity;
            const double beta = alpha * std::pow(tanh(alpha), -0.5);
            d_wave_number[i] = (alpha + std::pow(beta, 2) * std::pow(cosh(beta), -2)) /
                               (d_depth * (tanh(beta) + beta * std::pow(cosh(beta), -2)));

            double spectral_density = std::numeric_limits<double>::quiet_NaN();
            ;
            if (d_wave_spectrum == "JONSWAP")
            {
                double sigma, A, gamma = 3.3;
                ;
                if (d_omega[i] <= omega_s)
                    sigma = 0.07;
                else
                    sigma = 0.09;

                A = std::exp(-std::pow((d_omega[i] / omega_s - 1) / (sigma * std::sqrt(2)), 2));

                spectral_density = 320 * std::pow(d_Hs, 2) / (std::pow(d_Ts, 4) * std::pow(d_omega[i], 5)) *
                                   std::exp(-1950.0 / (std::pow(d_Ts * d_omega[i], 4))) * std::pow(gamma, A);
            }
            else if (d_wave_spectrum == "BRETSCHNEIDER")
            {
                spectral_density = 173 * std::pow(d_Hs, 2) / (std::pow(d_Ts, 4) * std::pow(d_omega[i], 5)) *
                                   std::exp(-692.0 / (std::pow(d_Ts * d_omega[i], 4)));
            }
            else
            {
                TBOX_ERROR(
                    "IrregularWaveBcCoef::IrregularWaveBcCoef(): Unknown wave "
                    "spectrum type "
                    << d_wave_spectrum
                    << " .This class supports only JONSWAP "
                       "and BRETSCHNEIDER wave spectra.");
            }

            d_amplitude[i] = std::sqrt(2.0 * spectral_density * delta_omega);
        }
    }

    if (!IBTK_MPI::getRank())
    {
        std::ofstream wave_stream;
        wave_stream.open("irregular_wave.txt", std::fstream::out);
        wave_stream.precision(10);

        for (int i = 0; i < d_num_waves; ++i)
        {
            wave_stream << d_amplitude[i] << "\t" << d_omega[i] << "\t" << d_wave_number[i] << "\t" << d_phase[i]
                        << std::endl;
        }

        wave_stream.close();
    }

    return;
} // IrregularWaveBcCoef

IrregularWaveBcCoef::~IrregularWaveBcCoef()
{
    // intentionally left-blank.
    return;
} // ~IrregularWaveBcCoef

void
IrregularWaveBcCoef::setBcCoefs(Pointer<ArrayData<NDIM, double> >& acoef_data,
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
IrregularWaveBcCoef::numberOfExtensionsFillable() const
{
    return EXTENSIONS_FILLABLE;
} // numberOfExtensionsFillable

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IrregularWaveBcCoef::getFromInput(Pointer<Database> input_db)
{
    Pointer<Database> wave_db = input_db->getDatabase("wave_parameters_db");
#if !defined(NDEBUG)
    TBOX_ASSERT(input_db->isDatabase("wave_parameters_db"));
#endif

    d_gravity = wave_db->getDouble("gravitational_constant");
    d_depth = wave_db->getDouble("depth");
    d_num_waves = wave_db->getIntegerWithDefault("num_waves", d_num_waves);
    d_omega_begin = wave_db->getDouble("omega_begin");
    d_omega_end = wave_db->getDouble("omega_end");
    d_Ts = wave_db->getDouble("significant_wave_period");
    d_Hs = wave_db->getDouble("significant_wave_height");
    d_wave_spectrum = wave_db->getString("wave_spectrum");
    d_num_interface_cells = wave_db->getDouble("num_interface_cells");

    return;
} // getFromInput

double
IrregularWaveBcCoef::getSurfaceElevation(double x, double time) const
{
    double eta = 0;
    for (int i = 0; i < d_num_waves; i++)
    {
        const double theta = d_wave_number[i] * x - d_omega[i] * time + d_phase[i];
        eta += d_amplitude[i] * cos(theta);
    }

    return eta;
} // getSurfaceElevation

double
IrregularWaveBcCoef::getVelocity(double x, double z_plus_d, double time) const
{
    if (d_comp_idx == 0)
    {
        double velocity_component = 0.0;
        for (int i = 0; i < d_num_waves; i++)
        {
            const double theta = d_wave_number[i] * x - d_omega[i] * time + d_phase[i];
            velocity_component += d_amplitude[i] * d_omega[i] * cosh(d_wave_number[i] * z_plus_d) * cos(theta) /
                                  sinh(d_wave_number[i] * d_depth);
        }
        return velocity_component;
    }
    if (d_comp_idx == 1)
    {
#if (NDIM == 2)
        double velocity_component = 0.0;
        for (int i = 0; i < d_num_waves; i++)
        {
            const double theta = d_wave_number[i] * x - d_omega[i] * time + d_phase[i];
            velocity_component += d_amplitude[i] * d_omega[i] * sinh(d_wave_number[i] * z_plus_d) * sin(theta) /
                                  sinh(d_wave_number[i] * d_depth);
        }
        return velocity_component;

#elif (NDIM == 3)
        return 0.0;
#endif
    }
#if (NDIM == 3)
    if (d_comp_idx == 2)
    {
        double velocity_component = 0.0;
        for (int i = 0; i < d_num_waves; i++)
        {
            const double theta = d_wave_number[i] * x - d_omega[i] * time + d_phase[i];
            velocity_component += d_amplitude[i] * d_omega[i] * sinh(d_wave_number[i] * z_plus_d) * sin(theta) /
                                  sinh(d_wave_number[i] * d_depth);
        }
        return velocity_component;
    }
#endif

    return std::numeric_limits<double>::signaling_NaN();
    ;
} // getVelocity

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR
