// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_RNG
#define included_IBAMR_RNG

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

namespace IBAMR
{
/*!
 * \brief Class RNG organizes functions that provide random-number generator
 * functionality.
 */
class RNG
{
public:
    static void srandgen(unsigned long seed);

    static void genrand(double* rn);

    static void genrandn(double* result);

    static void parallel_seed(int global_seed);

private:
    RNG() = delete;
    RNG(RNG&) = delete;
    ~RNG() = delete;
    RNG& operator=(RNG&) = delete;
};
} // namespace IBAMR

#endif //#ifndef included_IBAMR_RNG
