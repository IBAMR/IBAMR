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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_WaveDampingFunctions
#define included_WaveDampingFunctions

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/WaveUtilities.h"

namespace IBAMR
{
/*!
 * A collection of post processing call back functions to be hooked into
 * IBAMR::HierarchyIntegrator class to employ a wave damping zone of some
 * prescribed length within the computational domain. Additional strategies
 * can be implemented as simple C functions within this collection.
 *
 * The instance of a WaveDampingStrategy and the particular damping function
 * should be registered with an appropriate hierarchy integrator via
 * registerPostprocessIntegrateHierarchyCallback.
 *
 * \param ctx is the pointer to WaveDampingData struct.
 */

namespace WaveDampingFunctions
{
void callRelaxationZoneCallbackFunction(double current_time,
                                        double new_time,
                                        bool skip_synchronize_new_state_data,
                                        int num_cycles,
                                        void* ctx);

void callConservedWaveAbsorbingCallbackFunction(double current_time,
                                                double new_time,
                                                bool skip_synchronize_new_state_data,
                                                int num_cycles,
                                                void* ctx);
} // namespace WaveDampingFunctions

} // namespace IBAMR

#endif // #ifndef included_WaveDampingFunctions
