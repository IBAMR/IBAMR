// Filename: StokesWaveGeneratorStrategy.cpp
// Created on 4 Nov 2019 by Amneet Bhalla
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

#include "ibamr/StokesWaveGeneratorStrategy.h"
#include "ibamr/app_namespaces.h"

#include "ibtk/ibtk_utilities.h"

#include "tbox/Database.h"

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
