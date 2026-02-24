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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeDoubleLinearTimeInterpolateOp
#define included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeDoubleLinearTimeInterpolateOp

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/EdgeDoubleLinearTimeInterpolateOp.h>)
#include <SAMRAI/pdat/EdgeDoubleLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeDoubleLinearTimeInterpolateOp 1
#else
#include <EdgeDoubleLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeDoubleLinearTimeInterpolateOp 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeDoubleLinearTimeInterpolateOp)
using SAMRAIEdgeDoubleLinearTimeInterpolateOp = SAMRAI::pdat::EdgeDoubleLinearTimeInterpolateOp;
#else
using SAMRAIEdgeDoubleLinearTimeInterpolateOp = SAMRAI::pdat::EdgeDoubleLinearTimeInterpolateOp<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_EdgeDoubleLinearTimeInterpolateOp

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIEdgeDoubleLinearTimeInterpolateOp
