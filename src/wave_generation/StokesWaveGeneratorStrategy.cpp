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

#include "ibamr/StokesWaveGeneratorStrategy.h"

#include "tbox/Database.h"

#include <limits>

#include "ibamr/app_namespaces.h"

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StokesWaveGeneratorStrategy::StokesWaveGeneratorStrategy(const std::string& object_name, Pointer<Database> input_db)
{
    d_object_name = object_name;
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_object_name.empty());
    TBOX_ASSERT(input_db);
#endif

    // Get wave parameters.
    getFromInput(input_db);

    return;
} // StokesWaveGeneratorStrategy

double
StokesWaveGeneratorStrategy::getWaterDepth() const
{
    return d_depth;
} // getWaterDept

double
StokesWaveGeneratorStrategy::getWaveAngularFrequency() const
{
    return d_omega;
} // getWaveAngularFrequency

double
StokesWaveGeneratorStrategy::getWaveNumber() const
{
    return d_wave_number;
} // getWaveNumber

double
StokesWaveGeneratorStrategy::getWaveAmplitude() const
{
    return d_amplitude;
} // getWaveAmplitude

double
StokesWaveGeneratorStrategy::getGravity() const
{
    return d_gravity;
} // getGravity

/////////////////////////////// PRIVATE //////////////////////////////////////
void
StokesWaveGeneratorStrategy::getFromInput(Pointer<Database> input_db)
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
    d_wave_gen_data.d_num_interface_cells = wave_db->getDouble("num_interface_cells");
    d_deep_water_limit = wave_db->getBoolWithDefault("deep_water_limit", d_deep_water_limit);

    return;
} // getFromInput

} // namespace IBAMR
