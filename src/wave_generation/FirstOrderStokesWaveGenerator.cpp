// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
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

#include "ibamr/FirstOrderStokesWaveGenerator.h"

#include "tbox/Database.h"

#include <cmath>
#include <limits>

#include "ibamr/app_namespaces.h" // IWYU pragma: keep

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FirstOrderStokesWaveGenerator::FirstOrderStokesWaveGenerator(const std::string& object_name, Pointer<Database> input_db)
    : StokesWaveGeneratorStrategy(object_name, input_db)
{
    return;
} // FirstOrderStokesWaveGenerator

double
FirstOrderStokesWaveGenerator::getSurfaceElevation(double x, double time) const
{
    // This computes eta for first order water waves.

    const double theta = d_wave_number * x - d_omega * time;

    return (d_amplitude * cos(theta));

} // getSurfaceElevation

double
FirstOrderStokesWaveGenerator::getVelocity(double x, double z_plus_d, double time, int comp_idx) const
{
    const double theta = d_wave_number * x - d_omega * time;

    if (comp_idx == 0)
    {
        double velocity_component = 0.0;

        velocity_component =
            d_amplitude * d_omega * cosh(d_wave_number * z_plus_d) * cos(theta) / sinh(d_wave_number * d_depth);
        return velocity_component;
    }
    if (comp_idx == 1)
    {
#if (NDIM == 2)
        double velocity_component = 0.0;

        velocity_component =
            d_amplitude * d_omega * sinh(d_wave_number * z_plus_d) * sin(theta) / sinh(d_wave_number * d_depth);
        return velocity_component;

#elif (NDIM == 3)
        return 0.0;
#endif
    }
#if (NDIM == 3)
    if (comp_idx == 2)
    {
        double velocity_component = 0.0;

        velocity_component =
            d_amplitude * d_omega * sinh(d_wave_number * z_plus_d) * sin(theta) / sinh(d_wave_number * d_depth);

        return velocity_component;
    }
#endif

    return std::numeric_limits<double>::signaling_NaN();

} // getVelocity

} // namespace IBAMR
