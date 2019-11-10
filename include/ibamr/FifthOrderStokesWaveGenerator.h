// Filename: FifthOrderStokesWaveGenerator.h
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

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_FifthOrderStokesWaveGenerator
#define included_FifthOrderStokesWaveGenerator

///////////////////////////// INCLUDES ///////////////////////////////////

#include "ibamr/StokesWaveGeneratorStrategy.h"

namespace IBAMR
{
/*!
 * \brief Class for generating fifth order Stokes wave using Fenton's (Stokes) wave theory.
 */
class FifthOrderStokesWaveGenerator : public StokesWaveGeneratorStrategy
{
public:
    FifthOrderStokesWaveGenerator(const std::string& object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * Get surface elevation at a specified horizontal position and time.
     */
    double getSurfaceElevation(double x, double time) const;

    /*!
     * Get velocity component at a specified position and time.
     */
    double getVelocity(double x, double z_plus_d, double time, int comp_idx) const;

private:
    /*!
     * Initialize Stokes coefficients.
     */
    void initStokesCoefficients();

    /*!
     * Stokes coefficients.
     */
    double d_A[6][6], d_B[6][6], d_C[5];
    double d_p[5], d_eta[5];
};

} // namespace IBAMR

#endif included_FifthOrderStokesWaveGenerator
