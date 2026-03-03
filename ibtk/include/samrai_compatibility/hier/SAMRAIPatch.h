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

#ifndef included_IBTK_samrai_compatibility_hier_SAMRAIPatch
#define included_IBTK_samrai_compatibility_hier_SAMRAIPatch

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/hier/Patch.h>)
#include <SAMRAI/hier/Patch.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_PATCH 1
#else
#include <Patch.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_PATCH 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_PATCH)
using SAMRAIPatch = SAMRAI::hier::Patch;
#else
using SAMRAIPatch = SAMRAI::hier::Patch<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_PATCH

#endif // #ifndef included_IBTK_samrai_compatibility_hier_SAMRAIPatch
