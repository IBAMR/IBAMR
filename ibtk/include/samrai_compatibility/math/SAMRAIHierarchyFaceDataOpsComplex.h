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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyFaceDataOpsComplex
#define included_IBTK_samrai_compatibility_math_SAMRAIHierarchyFaceDataOpsComplex

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/math/HierarchyFaceDataOpsComplex.h")
#include <SAMRAI/math/HierarchyFaceDataOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsComplex 1
#else
#include <HierarchyFaceDataOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsComplex 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsComplex)
using SAMRAIHierarchyFaceDataOpsComplex = SAMRAI::math::HierarchyFaceDataOpsComplex;
#else
using SAMRAIHierarchyFaceDataOpsComplex = SAMRAI::math::HierarchyFaceDataOpsComplex<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsComplex

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyFaceDataOpsComplex
