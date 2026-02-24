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

#ifndef included_IBTK_samrai_compatibility_hier_SAMRAIIndex
#define included_IBTK_samrai_compatibility_hier_SAMRAIIndex

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/hier/Index.h")
#include <SAMRAI/hier/Index.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_INDEX 1
#else
#include <Index.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_INDEX 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_INDEX)
using SAMRAIIndex = SAMRAI::hier::Index;
#else
using SAMRAIIndex = SAMRAI::hier::Index<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_INDEX

#endif // #ifndef included_IBTK_samrai_compatibility_hier_SAMRAIIndex
