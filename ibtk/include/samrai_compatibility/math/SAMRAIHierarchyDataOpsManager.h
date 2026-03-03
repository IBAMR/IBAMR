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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyDataOpsManager
#define included_IBTK_samrai_compatibility_math_SAMRAIHierarchyDataOpsManager

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/math/HierarchyDataOpsManager.h>)
#include <SAMRAI/math/HierarchyDataOpsManager.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyDataOpsManager 1
#else
#include <HierarchyDataOpsManager.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyDataOpsManager 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyDataOpsManager)
using SAMRAIHierarchyDataOpsManager = SAMRAI::math::HierarchyDataOpsManager;
#else
using SAMRAIHierarchyDataOpsManager = SAMRAI::math::HierarchyDataOpsManager<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyDataOpsManager

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyDataOpsManager
