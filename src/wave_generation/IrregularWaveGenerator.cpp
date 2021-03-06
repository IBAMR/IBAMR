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

#include "ibamr/IrregularWaveGenerator.h"
#include "ibamr/RNG.h"

#include "ibtk/IBTK_MPI.h"

#include "tbox/Database.h"
#include "tbox/Utilities.h"

#include <algorithm>
#include <cmath>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////
namespace
{
static const unsigned SEED = 1234567;
}
/////////////////////////////// PUBLIC ///////////////////////////////////////

IrregularWaveGenerator::IrregularWaveGenerator(const std::string& object_name, Pointer<Database> input_db)
    : StokesWaveGeneratorStrategy(object_name, input_db)
{
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
    for (int i = 0; i < d_num_waves; i++)
    {
        double rn;
        RNG::genrand(&rn);
        d_phase[i] = 2 * M_PI * rn;

        d_omega[i] = d_omega_begin + i * delta_omega;

        // Using an approximate formula for the dispersion relationship, calculate the wave number.
        // See Eqn. (5.4.22) in WAVES IN OCEANIC AND COASTAL WATERS by LEO H. HOLTHUIJSEN.
        const double alpha = std::pow(d_omega[i], 2) * d_depth / d_gravity;
        const double beta = alpha * std::pow(tanh(alpha), -0.5);
        d_wave_number[i] = (alpha + std::pow(beta, 2) * std::pow(cosh(beta), -2)) /
                           (d_depth * (tanh(beta) + beta * std::pow(cosh(beta), -2)));

        double spectral_density = std::numeric_limits<double>::quiet_NaN();

        if (d_wave_spectrum == "JONSWAP")
        {
            double sigma, A, gamma = 3.3;

            sigma = d_omega[i] <= omega_s ? 0.07 : 0.09;

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
            TBOX_ERROR("IrregularWaveBcCoef::IrregularWaveBcCoef(): Unknown wave spectrum type "
                       << d_wave_spectrum << " .This class supports only JONSWAP and BRETSCHNEIDER wave spectra.");
        }
        d_amplitude[i] = std::sqrt(2.0 * spectral_density * delta_omega);
    }

    return;
} // IrregularWaveGenerator

double
IrregularWaveGenerator::getSurfaceElevation(const double x, const double time) const
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
IrregularWaveGenerator::getVelocity(const double x, const double z_plus_d, const double time, const int comp_idx) const
{
    if (comp_idx == 0)
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
    if (comp_idx == 1)
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
    if (comp_idx == 2)
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

} // getVelocity

void
IrregularWaveGenerator::printWaveData(ofstream& ostream) const
{
    if (!IBTK_MPI::getRank())
    {
        for (int i = 0; i < d_num_waves; ++i)
        {
            ostream << d_amplitude[i] << "\t" << d_omega[i] << "\t" << d_wave_number[i] << "\t" << d_phase[i]
                    << std::endl;
        }
    }

    return;
} // getWaveData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IrregularWaveGenerator::getFromInput(Pointer<Database> input_db)
{
    Pointer<Database> wave_db = input_db;
    if (input_db->isDatabase("wave_parameters_db"))
    {
        wave_db = input_db->getDatabase("wave_parameters_db");
    }

    d_num_waves = wave_db->getIntegerWithDefault("num_waves", d_num_waves);
    d_omega_begin = wave_db->getDouble("omega_begin");
    d_omega_end = wave_db->getDouble("omega_end");
    d_Hs = wave_db->getDouble("significant_wave_height");
    d_Ts = wave_db->getDouble("significant_wave_period");
    d_wave_spectrum = wave_db->getString("wave_spectrum");

    return;
} // getFromInput

} // namespace IBAMR
