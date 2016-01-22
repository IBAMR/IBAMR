// Filename: compiler_hints.h
// Created on 17 Dec 2009 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_compiler_hints
#define included_compiler_hints

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "IBTK_config.h"

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

#endif //#ifndef included_compiler_hints
