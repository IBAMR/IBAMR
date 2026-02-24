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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyFaceDataOpsInteger
#define included_IBTK_samrai_compatibility_math_SAMRAIHierarchyFaceDataOpsInteger

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/math/HierarchyFaceDataOpsInteger.h>)
#include <SAMRAI/math/HierarchyFaceDataOpsInteger.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsInteger 1
#else
#include <HierarchyFaceDataOpsInteger.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsInteger 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsInteger)
using SAMRAIHierarchyFaceDataOpsInteger = SAMRAI::math::HierarchyFaceDataOpsInteger;
#else
using SAMRAIHierarchyFaceDataOpsInteger = SAMRAI::math::HierarchyFaceDataOpsInteger<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyFaceDataOpsInteger

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyFaceDataOpsInteger
