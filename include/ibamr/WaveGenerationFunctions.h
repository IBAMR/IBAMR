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

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_WaveGenerationFunctions
#define included_WaveGenerationFunctions

///////////////////////////// INCLUDES ///////////////////////////////////

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
