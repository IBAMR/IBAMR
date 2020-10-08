// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_FirstOrderStokesWaveGenerator
#define included_FirstOrderStokesWaveGenerator

///////////////////////////// INCLUDES ///////////////////////////////////

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
 * \brief Class for generating first order water waves based upon linear wave theory.
 */
class FirstOrderStokesWaveGenerator : public StokesWaveGeneratorStrategy
{
public:
    FirstOrderStokesWaveGenerator(const std::string& object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * Get surface elevation at a specified horizontal position and time.
     */
    double getSurfaceElevation(double x, double time) const override;

    /*!
     * Get velocity component at a specified position and time.
     */
    double getVelocity(double x, double z_plus_d, double time, int comp_idx) const override;
};

} // namespace IBAMR

#endif
