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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAISideIndex
#define included_IBTK_samrai_compatibility_pdat_SAMRAISideIndex

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/pdat/SideIndex.h")
#include <SAMRAI/pdat/SideIndex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideIndex 1
#else
#include <SideIndex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideIndex 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideIndex)
using SAMRAISideIndex = SAMRAI::pdat::SideIndex;
#else
using SAMRAISideIndex = SAMRAI::pdat::SideIndex<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideIndex

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAISideIndex
