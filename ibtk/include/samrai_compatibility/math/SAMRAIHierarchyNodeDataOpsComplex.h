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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyNodeDataOpsComplex
#define included_IBTK_samrai_compatibility_math_SAMRAIHierarchyNodeDataOpsComplex

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/math/HierarchyNodeDataOpsComplex.h>)
#include <SAMRAI/math/HierarchyNodeDataOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyNodeDataOpsComplex 1
#else
#include <HierarchyNodeDataOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyNodeDataOpsComplex 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyNodeDataOpsComplex)
using SAMRAIHierarchyNodeDataOpsComplex = SAMRAI::math::HierarchyNodeDataOpsComplex;
#else
using SAMRAIHierarchyNodeDataOpsComplex = SAMRAI::math::HierarchyNodeDataOpsComplex<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HierarchyNodeDataOpsComplex

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIHierarchyNodeDataOpsComplex
