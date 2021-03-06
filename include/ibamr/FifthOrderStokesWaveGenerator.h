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

#ifndef included_FifthOrderStokesWaveGenerator
#define included_FifthOrderStokesWaveGenerator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/StokesWaveGeneratorStrategy.h"

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
    double getSurfaceElevation(double x, double time) const override;

    /*!
     * Get velocity component at a specified position and time.
     */
    double getVelocity(double x, double z_plus_d, double time, int comp_idx) const override;

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

#endif
