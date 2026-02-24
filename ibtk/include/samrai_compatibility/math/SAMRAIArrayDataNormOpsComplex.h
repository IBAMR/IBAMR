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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIArrayDataNormOpsComplex
#define included_IBTK_samrai_compatibility_math_SAMRAIArrayDataNormOpsComplex

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/math/ArrayDataNormOpsComplex.h")
#include <SAMRAI/math/ArrayDataNormOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataNormOpsComplex 1
#else
#include <ArrayDataNormOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataNormOpsComplex 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataNormOpsComplex)
using SAMRAIArrayDataNormOpsComplex = SAMRAI::math::ArrayDataNormOpsComplex;
#else
using SAMRAIArrayDataNormOpsComplex = SAMRAI::math::ArrayDataNormOpsComplex<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataNormOpsComplex

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIArrayDataNormOpsComplex
