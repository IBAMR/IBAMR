// Filename: StokesWaveGeneratorStrategy.h
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

#ifndef included_StokesWaveGeneratorStrategy
#define included_StokesWaveGeneratorStrategy

///////////////////////////// INCLUDES ///////////////////////////////////
#include "ibamr/WaveUtilities.h"

#include <string>

namespace SAMRAI
{
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

namespace IBAMR
{
/*!
 * \brief Strategy class for generating Stokes wave.
 */
class StokesWaveGeneratorStrategy
{
public:
    /*!
     * \brief Constructor of the class.
     */
    StokesWaveGeneratorStrategy(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor of the class.
     */
    virtual ~StokesWaveGeneratorStrategy() = default;

    /*!
     * Get surface elevation at a specified horizontal position and time.
     */
    virtual double getSurfaceElevation(double x, double time) const = 0;

    /*!
     * Get velocity component at a specified position and time.
     *
     * \note We assume that the bottom left corner of the numerical wave tank
     * is at (0,0) in 2D or at (0,0,0) in 3D. The \param z_plus_d is the
     * \f$ y \f$ (in 2D) or the \f$ z \f$ (in 3D) coordinate of a point
     * in the numerical wave tank. In partcular \param z_plus_d is \em not with respect
     * to free-surface as generally taken in some textbooks.
     */
    virtual double getVelocity(double x, double z_plus_d, double time, int comp_idx) const = 0;

    /*
     * \brief Get Stokes wave parameters.
     */
    double getWaterDepth() const;

    double getWaveAngularFrequency() const;

    double getWaveNumber() const;

    double getWaveAmplitude() const;

    double getGravity() const;

    /*
     * \brief WaveGeneration data object.
     */
    WaveGenerationData d_wave_gen_data;

protected:
    /*!
     * Book-keeping.
     */
    std::string d_object_name;

    /*!
     * \brief Wave parameters.
     *
     * \param d_wave_number  : Wave number of dominant wave component [$2\pi/m$]
     * \param d_amplitude    : Amplitude of the dominant wave component [m]
     * \param d_depth        : Depth of water, from sea bed to still water level [m]
     * \param d_gravity      : Acceleration due to gravity [$m/s^2$]
     * \param d_omega        : Angular frequency [$2 \pi/s$] (optional)
     *
     * \NOTE The class calculates a more accurate value of omega from the expansion coefficients
     * and the provided value in not used.
     */
    double d_depth, d_omega, d_wave_number, d_amplitude, d_gravity;

    /*!
     * If we are calculating in deep water limit.
     */
    bool d_deep_water_limit = false;

    /*
     *  \brief Stokes wave generator data members.
     */
    /*!
     * Get wave parameters from input db.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

}; // StokesWaveGeneratorStrategy

} // namespace IBAMR

#endif // #ifndef included_StokesWaveGeneratorStrategy
