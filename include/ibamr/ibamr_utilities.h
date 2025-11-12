// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2025 by the IBAMR developers
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

#ifndef included_IBAMR_ibamr_utilities
#define included_IBAMR_ibamr_utilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "tbox/PIO.h"
#include "tbox/Pointer.h"
/////////////////////////////// MACRO DEFINITIONS ////////////////////////////

#define IBAMR_DO_ONCE(task)                                                                                            \
    do                                                                                                                 \
    {                                                                                                                  \
        static bool done = false;                                                                                      \
        if (done == false)                                                                                             \
        {                                                                                                              \
            task;                                                                                                      \
            done = true;                                                                                               \
        }                                                                                                              \
    } while (0);

namespace IBAMR
{
static const bool ENABLE_TIMERS = true;
}

#define IBAMR_TIMER_START(timer)                                                                                       \
    do                                                                                                                 \
    {                                                                                                                  \
        if (IBAMR::ENABLE_TIMERS) timer->start();                                                                      \
    } while (0);

#define IBAMR_TIMER_STOP(timer)                                                                                        \
    do                                                                                                                 \
    {                                                                                                                  \
        if (IBAMR::ENABLE_TIMERS) timer->stop();                                                                       \
    } while (0);

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace std
{
template <typename T>
struct less<SAMRAI::tbox::Pointer<T> >
{
    inline bool operator()(const SAMRAI::tbox::Pointer<T>& k1, const SAMRAI::tbox::Pointer<T>& k2) const
    {
        return k1.getPointer() < k2.getPointer();
    }
};
} // namespace std

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_ibamr_utilities
