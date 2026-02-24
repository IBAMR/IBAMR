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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIOuterfaceFloatLinearTimeInterpolateOp
#define included_IBTK_samrai_compatibility_pdat_SAMRAIOuterfaceFloatLinearTimeInterpolateOp

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/pdat/OuterfaceFloatLinearTimeInterpolateOp.h")
#include <SAMRAI/pdat/OuterfaceFloatLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceFloatLinearTimeInterpolateOp 1
#else
#include <OuterfaceFloatLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceFloatLinearTimeInterpolateOp 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceFloatLinearTimeInterpolateOp)
using SAMRAIOuterfaceFloatLinearTimeInterpolateOp = SAMRAI::pdat::OuterfaceFloatLinearTimeInterpolateOp;
#else
using SAMRAIOuterfaceFloatLinearTimeInterpolateOp = SAMRAI::pdat::OuterfaceFloatLinearTimeInterpolateOp<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceFloatLinearTimeInterpolateOp

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIOuterfaceFloatLinearTimeInterpolateOp
