// Filename: ibtk_macros.h
// Created on 13 June 2019 by David Wells
//
// Copyright (c) 2019, Boyce Griffith
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
#ifndef included_IBTK_ibtk_macros
#define included_IBTK_ibtk_macros

#include <IBTK_config.h>

#ifdef IBTK_HAVE_PRAGMA_KEYWORD
// Prevent clang-format from doing strange things to this very long macro:
// clang-format off

// The first four warnings here should be left in that order: new warnings
// should be placed at the end.
#define IBTK_DISABLE_EXTRA_WARNINGS                             \
_Pragma("GCC diagnostic push")                                  \
_Pragma("GCC diagnostic ignored \"-Wunknown-pragmas\"")         \
_Pragma("GCC diagnostic ignored \"-Wpragmas\"")                 \
_Pragma("GCC diagnostic ignored \"-Wunknown-warning-option\"")  \
_Pragma("GCC diagnostic ignored \"-Wunknown-warning\"")         \
_Pragma("GCC diagnostic ignored \"-Wunused-variable\"")         \
_Pragma("GCC diagnostic ignored \"-Wignored-attributes\"")      \
_Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"") \
_Pragma("GCC diagnostic ignored \"-Wmisleading-indentation\"")  \
_Pragma("GCC diagnostic ignored \"-Wint-in-bool-context\"")     \
_Pragma("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")     \
_Pragma("GCC diagnostic ignored \"-Wunused-local-typedefs\"")   \
_Pragma("GCC diagnostic ignored \"-Wdeprecated-copy\"")         \
_Pragma("GCC diagnostic ignored \"-Wunneeded-internal-declaration\"")

#define IBTK_ENABLE_EXTRA_WARNINGS _Pragma("GCC diagnostic pop")

// clang-format on
#else

#define IBTK_DISABLE_EXTRA_WARNINGS
#define IBTK_ENABLE_EXTRA_WARNINGS

#endif // #ifdef IBTK_HAVE_PRAGMA_KEYWORD

#endif // #ifndef included_IBTK_ibtk_macros
