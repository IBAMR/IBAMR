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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyCellDataOpsComplex
#define included_IBTK_samrai_compatibility_math_SAMRAIHierarchyCellDataOpsComplex

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/math/HierarchyCellDataOpsComplex.h")
#include <SAMRAI/math/HierarchyCellDataOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyCellDataOpsComplex 1
#else
#include <HierarchyCellDataOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyCellDataOpsComplex 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyCellDataOpsComplex)
using SAMRAIHierarchyCellDataOpsComplex = SAMRAI::math::HierarchyCellDataOpsComplex;
#else
using SAMRAIHierarchyCellDataOpsComplex = SAMRAI::math::HierarchyCellDataOpsComplex<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyCellDataOpsComplex

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyCellDataOpsComplex
