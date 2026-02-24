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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchNodeDataBasicOps
#define included_IBTK_samrai_compatibility_math_SAMRAIPatchNodeDataBasicOps

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/math/PatchNodeDataBasicOps.h>)
#include <SAMRAI/math/PatchNodeDataBasicOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchNodeDataBasicOps 1
#else
#include <PatchNodeDataBasicOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchNodeDataBasicOps 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class... Args>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchNodeDataBasicOps)
using SAMRAIPatchNodeDataBasicOps = SAMRAI::math::PatchNodeDataBasicOps<Args...>;
#else
using SAMRAIPatchNodeDataBasicOps = SAMRAI::math::PatchNodeDataBasicOps<NDIM, Args...>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchNodeDataBasicOps

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchNodeDataBasicOps
