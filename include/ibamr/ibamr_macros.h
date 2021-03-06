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

#ifndef included_IBAMR_ibamr_macros
#define included_IBAMR_ibamr_macros

#include <ibamr/config.h>

#define IBAMR_VERSION_GTE(major, minor, subminor)                                                                      \
    ((IBAMR_VERSION_MAJOR * 10000 + IBAMR_VERSION_MINOR * 100 + IBAMR_VERSION_SUBMINOR) >=                             \
     (major)*10000 + (minor)*100 + (subminor))

#endif // #ifndef included_IBAMR_ibamr_macros
