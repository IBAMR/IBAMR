// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
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

#ifndef included_IrregularWaveGenerator
#define included_IrregularWaveGenerator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/StokesWaveGeneratorStrategy.h"

#include "tbox/Pointer.h"

#include <fstream>
#include <limits>
#include <string>
#include <vector>

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
 * \brief Class for generating Irregular waves.
 */
class IrregularWaveGenerator : public StokesWaveGeneratorStrategy
{
public:
    IrregularWaveGenerator(const std::string& object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * Get surface elevation at a specified horizontal position and time.
     */
    double getSurfaceElevation(double x, double time) const override;

    /*!
     * Get velocity component at a specified position and time.
     */
    double getVelocity(double x, double z_plus_d, double time, int comp_idx) const override;

    /*!
     * Print the wave data.
     */
    void printWaveData(std::ofstream& ostream) const;

private:
    /*!
     * Get wave parameters from input db.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    ///
    /// Number of component waves with random phases to be generated (default = 50).
    ///
    int d_num_waves = 50;

    ///
    /// Significant wave height [length].
    ///
    double d_Hs = std::numeric_limits<double>::quiet_NaN();

    ///
    /// Significant wave period [time].
    ///
    double d_Ts = std::numeric_limits<double>::quiet_NaN();

    ///
    /// Lowest angular frequency in the spectrum [rad/time].
    ///
    double d_omega_begin = std::numeric_limits<double>::quiet_NaN();

    ///
    /// Highest angular frequency in the spectrum [rad/time].
    ///
    double d_omega_end = std::numeric_limits<double>::quiet_NaN();

    ///
    /// JONSWAP/Bretschneider wave spectrum.
    ///
    std::string d_wave_spectrum;

    ///
    /// Angular frequencies of component waves [rad/time].
    ///
    std::vector<double> d_omega;

    ///
    /// Wave number of component waves [2\f$ \pi \f$/length].
    ///
    std::vector<double> d_wave_number;

    ///
    /// Amplitude of component waves [length].
    ///
    std::vector<double> d_amplitude;

    ///
    /// Phase (random) of component waves [rad].
    ///
    std::vector<double> d_phase;
};

} // namespace IBAMR
#endif
