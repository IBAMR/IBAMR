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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeData
#define included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeData

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/pdat/EdgeData.h")
#include <SAMRAI/pdat/EdgeData.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeData 1
#else
#include <EdgeData.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeData 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class T>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeData)
using SAMRAIEdgeData = SAMRAI::pdat::EdgeData<T>;
#else
using SAMRAIEdgeData = SAMRAI::pdat::EdgeData<NDIM, T>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeData

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeData
