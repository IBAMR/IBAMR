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

#ifndef included_WaveGenerationFunctions
#define included_WaveGenerationFunctions

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/WaveUtilities.h"

namespace IBAMR
{
/*!
 * A collection of post processing call back functions to be hooked into
 * IBAMR::HierarchyIntegrator class to employ a wave generation zone of some
 * prescribed length within the computational domain. Additional strategies
 * can be implemented as simple C functions within this collection.
 *
 * The instance of a WaveGenerationStrategy and the particular generation function
 * should be registered with an appropriate hierarchy integrator via
 * registerPostprocessIntegrateHierarchyCallback.
 *
 * \param ctx is the pointer to the relevant wave generator object.
 */

namespace WaveGenerationFunctions
{
void callStokesWaveRelaxationCallbackFunction(double current_time,
                                              double new_time,
                                              bool skip_synchronize_new_state_data,
                                              int num_cycles,
                                              void* ctx);
} // namespace WaveGenerationFunctions

} // namespace IBAMR

#endif // #ifndef included_WaveGenerationFunctions
