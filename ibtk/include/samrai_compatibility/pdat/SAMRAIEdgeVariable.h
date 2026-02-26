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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeVariable
#define included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeVariable

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/EdgeVariable.h>)
#include <SAMRAI/pdat/EdgeVariable.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeVariable 1
#else
#include <EdgeVariable.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeVariable 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class T>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeVariable)
using SAMRAIEdgeVariable = SAMRAI::pdat::EdgeVariable<T>;
#else
using SAMRAIEdgeVariable = SAMRAI::pdat::EdgeVariable<NDIM, T>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeVariable

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeVariable
