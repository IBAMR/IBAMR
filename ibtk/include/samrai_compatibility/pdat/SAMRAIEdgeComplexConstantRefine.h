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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeComplexConstantRefine
#define included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeComplexConstantRefine

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/EdgeComplexConstantRefine.h>)
#include <SAMRAI/pdat/EdgeComplexConstantRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeComplexConstantRefine 1
#else
#include <EdgeComplexConstantRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeComplexConstantRefine 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeComplexConstantRefine)
using SAMRAIEdgeComplexConstantRefine = SAMRAI::pdat::EdgeComplexConstantRefine;
#else
using SAMRAIEdgeComplexConstantRefine = SAMRAI::pdat::EdgeComplexConstantRefine<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeComplexConstantRefine

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeComplexConstantRefine
