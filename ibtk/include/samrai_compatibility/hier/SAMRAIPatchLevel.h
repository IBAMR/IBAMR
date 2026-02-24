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

#ifndef included_IBTK_samrai_compatibility_hier_SAMRAIPatchLevel
#define included_IBTK_samrai_compatibility_hier_SAMRAIPatchLevel

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/hier/PatchLevel.h")
#include <SAMRAI/hier/PatchLevel.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_PATCH_LEVEL 1
#else
#include <PatchLevel.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_PATCH_LEVEL 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_PATCH_LEVEL)
using SAMRAIPatchLevel = SAMRAI::hier::PatchLevel;
#else
using SAMRAIPatchLevel = SAMRAI::hier::PatchLevel<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_PATCH_LEVEL

#endif // #ifndef included_IBTK_samrai_compatibility_hier_SAMRAIPatchLevel
