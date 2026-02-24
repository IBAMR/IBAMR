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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAICellOverlap
#define included_IBTK_samrai_compatibility_pdat_SAMRAICellOverlap

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/CellOverlap.h>)
#include <SAMRAI/pdat/CellOverlap.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellOverlap 1
#else
#include <CellOverlap.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellOverlap 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellOverlap)
using SAMRAICellOverlap = SAMRAI::pdat::CellOverlap;
#else
using SAMRAICellOverlap = SAMRAI::pdat::CellOverlap<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellOverlap

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAICellOverlap
