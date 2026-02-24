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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataNormOpsComplex
#define included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataNormOpsComplex

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/math/PatchCellDataNormOpsComplex.h")
#include <SAMRAI/math/PatchCellDataNormOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataNormOpsComplex 1
#else
#include <PatchCellDataNormOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataNormOpsComplex 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataNormOpsComplex)
using SAMRAIPatchCellDataNormOpsComplex = SAMRAI::math::PatchCellDataNormOpsComplex;
#else
using SAMRAIPatchCellDataNormOpsComplex = SAMRAI::math::PatchCellDataNormOpsComplex<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataNormOpsComplex

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataNormOpsComplex
