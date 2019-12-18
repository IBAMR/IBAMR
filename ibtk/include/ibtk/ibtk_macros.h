// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

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
_Pragma("GCC diagnostic ignored \"-Wunused-parameter\"")        \
_Pragma("GCC diagnostic ignored \"-Wunneeded-internal-declaration\"")

#define IBTK_ENABLE_EXTRA_WARNINGS _Pragma("GCC diagnostic pop")

// clang-format on
#else

#define IBTK_DISABLE_EXTRA_WARNINGS
#define IBTK_ENABLE_EXTRA_WARNINGS

#endif // #ifdef IBTK_HAVE_PRAGMA_KEYWORD

#endif // #ifndef included_IBTK_ibtk_macros
