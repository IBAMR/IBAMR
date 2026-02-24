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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeFloatLinearTimeInterpolateOp
#define included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeFloatLinearTimeInterpolateOp

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/EdgeFloatLinearTimeInterpolateOp.h>)
#include <SAMRAI/pdat/EdgeFloatLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeFloatLinearTimeInterpolateOp 1
#else
#include <EdgeFloatLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeFloatLinearTimeInterpolateOp 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeFloatLinearTimeInterpolateOp)
using SAMRAIEdgeFloatLinearTimeInterpolateOp = SAMRAI::pdat::EdgeFloatLinearTimeInterpolateOp;
#else
using SAMRAIEdgeFloatLinearTimeInterpolateOp = SAMRAI::pdat::EdgeFloatLinearTimeInterpolateOp<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeFloatLinearTimeInterpolateOp

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeFloatLinearTimeInterpolateOp
