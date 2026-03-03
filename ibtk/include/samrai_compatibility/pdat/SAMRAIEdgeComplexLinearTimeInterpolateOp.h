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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeComplexLinearTimeInterpolateOp
#define included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeComplexLinearTimeInterpolateOp

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/EdgeComplexLinearTimeInterpolateOp.h>)
#include <SAMRAI/pdat/EdgeComplexLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeComplexLinearTimeInterpolateOp 1
#else
#include <EdgeComplexLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeComplexLinearTimeInterpolateOp 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeComplexLinearTimeInterpolateOp)
using SAMRAIEdgeComplexLinearTimeInterpolateOp = SAMRAI::pdat::EdgeComplexLinearTimeInterpolateOp;
#else
using SAMRAIEdgeComplexLinearTimeInterpolateOp = SAMRAI::pdat::EdgeComplexLinearTimeInterpolateOp<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeComplexLinearTimeInterpolateOp

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeComplexLinearTimeInterpolateOp
