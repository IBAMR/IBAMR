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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchEdgeDataOpsComplex
#define included_IBTK_samrai_compatibility_math_SAMRAIPatchEdgeDataOpsComplex

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/math/PatchEdgeDataOpsComplex.h>)
#include <SAMRAI/math/PatchEdgeDataOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchEdgeDataOpsComplex 1
#else
#include <PatchEdgeDataOpsComplex.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchEdgeDataOpsComplex 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchEdgeDataOpsComplex)
using SAMRAIPatchEdgeDataOpsComplex = SAMRAI::math::PatchEdgeDataOpsComplex;
#else
using SAMRAIPatchEdgeDataOpsComplex = SAMRAI::math::PatchEdgeDataOpsComplex<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchEdgeDataOpsComplex

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchEdgeDataOpsComplex
