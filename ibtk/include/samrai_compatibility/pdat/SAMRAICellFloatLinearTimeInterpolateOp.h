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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAICellFloatLinearTimeInterpolateOp
#define included_IBTK_samrai_compatibility_pdat_SAMRAICellFloatLinearTimeInterpolateOp

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/pdat/CellFloatLinearTimeInterpolateOp.h")
#include <SAMRAI/pdat/CellFloatLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellFloatLinearTimeInterpolateOp 1
#else
#include <CellFloatLinearTimeInterpolateOp.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellFloatLinearTimeInterpolateOp 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellFloatLinearTimeInterpolateOp)
using SAMRAICellFloatLinearTimeInterpolateOp = SAMRAI::pdat::CellFloatLinearTimeInterpolateOp;
#else
using SAMRAICellFloatLinearTimeInterpolateOp = SAMRAI::pdat::CellFloatLinearTimeInterpolateOp<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellFloatLinearTimeInterpolateOp

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAICellFloatLinearTimeInterpolateOp
