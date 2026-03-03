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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAINodeFloatLinearTimeInterpolateOp
#define included_IBTK_samrai_compatibility_pdat_SAMRAINodeFloatLinearTimeInterpolateOp

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/NodeFloatLinearTimeInterpolateOp.h>)
#include <SAMRAI/pdat/NodeFloatLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeFloatLinearTimeInterpolateOp 1
#else
#include <NodeFloatLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeFloatLinearTimeInterpolateOp 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeFloatLinearTimeInterpolateOp)
using SAMRAINodeFloatLinearTimeInterpolateOp = SAMRAI::pdat::NodeFloatLinearTimeInterpolateOp;
#else
using SAMRAINodeFloatLinearTimeInterpolateOp = SAMRAI::pdat::NodeFloatLinearTimeInterpolateOp<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeFloatLinearTimeInterpolateOp

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAINodeFloatLinearTimeInterpolateOp
