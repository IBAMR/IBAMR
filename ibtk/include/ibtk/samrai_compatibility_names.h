// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_samrai_compatibility_names
#define included_IBTK_samrai_compatibility_names

#include <ibtk/config.h>

#include <ibtk/samrai_compatibility_legacy_aliases.h>
#include <ibtk/samrai_compatibility_names.h>

#include <samrai_compatibility/samrai_compatibility_environment.h>
#include <samrai_compatibility/shim_aliases.h>

// Re-export SAMRAI compatibility aliases at global scope so converted code can
// refer to SAMRAI* names directly without an IBTK:: prefix.
using namespace IBTK::SAMRAICompatibility;

#endif // #ifndef included_IBTK_samrai_compatibility_names
