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

#ifndef included_IBAMR_config
#define included_IBAMR_config

#include <ibtk/config.h>

// Same as the IBTK include file - for compatibility we simply redefine all
// macros with IBAMR

#define IBAMR_VERSION_MAJOR IBTK_VERSION_MAJOR

#define IBAMR_VERSION_MINOR IBTK_VERSION_MINOR

#define IBAMR_VERSION_SUBMINOR IBTK_VERSION_SUBMINOR

#define IBAMR_VERSION_GTE(major, minor, subminor) IBTK_VERSION_GTE(major, minor, subminor)

#define IBAMR_HAVE_PRAGMA_KEYWORD IBTK_HAVE_PRAGMA_KEYWORD

#define IBAMR_HAVE_LIBMESH IBTK_HAVE_LIBMESH

#define IBAMR_HAVE_SILO IBTK_HAVE_SILO

#define IBAMR_FC_FUNC(name, NAME) IBTK_FC_FUNC(name, NAME)

#define IBAMR_FC_FUNC_(name, NAME) IBTK_FC_FUNC_(name, NAME)

#define IBAMR_DISABLE_EXTRA_WARNINGS IBTK_DISABLE_EXTRA_WARNINGS

#define IBAMR_ENABLE_EXTRA_WARNINGS IBTK_ENABLE_EXTRA_WARNINGS

#endif // included_IBAMR_config
