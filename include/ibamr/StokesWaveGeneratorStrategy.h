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

#ifndef included_StokesWaveGeneratorStrategy
#define included_StokesWaveGeneratorStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/WaveUtilities.h"

#include "tbox/Pointer.h"

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
