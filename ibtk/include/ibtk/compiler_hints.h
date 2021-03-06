// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_compiler_hints
#define included_IBTK_compiler_hints

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

/////////////////////////////// MACRO DEFINITIONS ////////////////////////////

namespace IBTK
{
#if IBTK_HAVE_BUILTIN_EXPECT
#define UNLIKELY(c) __builtin_expect(!!(c), 0)
#define LIKELY(c) __builtin_expect(!!(c), 1)
#else
#define UNLIKELY(c) (c)
#define LIKELY(c) (c)
#endif

#if IBTK_HAVE_BUILTIN_PREFETCH
#define PREFETCH_READ_NTA(a) __builtin_prefetch((a), 0, 0)
#define PREFETCH_READ_NTA_BLOCK(a, n)                                                                                  \
    do                                                                                                                 \
    {                                                                                                                  \
        for (int k = 0; k < static_cast<int>((n)); ++k)                                                                \
        {                                                                                                              \
            PREFETCH_READ_NTA((a) + k);                                                                                \
        }                                                                                                              \
    } while (0)
#if (NDIM == 2)
#define PREFETCH_READ_NTA_NDIM_BLOCK(a)                                                                                \
    do                                                                                                                 \
    {                                                                                                                  \
        PREFETCH_READ_NTA((a));                                                                                        \
        PREFETCH_READ_NTA((a) + 1);                                                                                    \
    } while (0)
#endif
#if (NDIM == 3)
#define PREFETCH_READ_NTA_NDIM_BLOCK(a)                                                                                \
    do                                                                                                                 \
    {                                                                                                                  \
        PREFETCH_READ_NTA((a));                                                                                        \
        PREFETCH_READ_NTA((a) + 1);                                                                                    \
        PREFETCH_READ_NTA((a) + 2);                                                                                    \
    } while (0)
#endif
#define PREFETCH_WRITE_NTA(a) __builtin_prefetch((a), 1, 0)
#define PREFETCH_WRITE_NTA_BLOCK(a, n)                                                                                 \
    do                                                                                                                 \
    {                                                                                                                  \
        for (int k = 0; k < static_cast<int>((n)); ++k)                                                                \
        {                                                                                                              \
            PREFETCH_WRITE_NTA((a) + k);                                                                               \
        }                                                                                                              \
    } while (0)
#if (NDIM == 2)
#define PREFETCH_WRITE_NTA_NDIM_BLOCK(a)                                                                               \
    do                                                                                                                 \
    {                                                                                                                  \
        PREFETCH_WRITE_NTA((a));                                                                                       \
        PREFETCH_WRITE_NTA((a) + 1);                                                                                   \
    } while (0)
#endif
#if (NDIM == 3)
#define PREFETCH_WRITE_NTA_NDIM_BLOCK(a)                                                                               \
    do                                                                                                                 \
    {                                                                                                                  \
        PREFETCH_WRITE_NTA((a));                                                                                       \
        PREFETCH_WRITE_NTA((a) + 1);                                                                                   \
        PREFETCH_WRITE_NTA((a) + 2);                                                                                   \
    } while (0)
#endif
#else
#define PREFETCH_READ_NTA(a)
#define PREFETCH_READ_NTA_BLOCK(a, n)
#define PREFETCH_READ_NTA_NDIM_BLOCK(a)
#define PREFETCH_WRITE_NTA(a)
#define PREFETCH_WRITE_NTA_BLOCK(a, n)
#define PREFETCH_WRITE_NTA_NDIM_BLOCK(a)
#endif

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_compiler_hints
