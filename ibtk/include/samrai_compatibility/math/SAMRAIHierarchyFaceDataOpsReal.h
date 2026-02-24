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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyFaceDataOpsReal
#define included_IBTK_samrai_compatibility_math_SAMRAIHierarchyFaceDataOpsReal

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/math/HierarchyFaceDataOpsReal.h")
#include <SAMRAI/math/HierarchyFaceDataOpsReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsReal 1
#else
#include <HierarchyFaceDataOpsReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsReal 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class T>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsReal)
using SAMRAIHierarchyFaceDataOpsReal = SAMRAI::math::HierarchyFaceDataOpsReal<T>;
#else
using SAMRAIHierarchyFaceDataOpsReal = SAMRAI::math::HierarchyFaceDataOpsReal<NDIM, T>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsReal

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyFaceDataOpsReal
