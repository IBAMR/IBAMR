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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchSideDataBasicOps
#define included_IBTK_samrai_compatibility_math_SAMRAIPatchSideDataBasicOps

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/math/PatchSideDataBasicOps.h>)
#include <SAMRAI/math/PatchSideDataBasicOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataBasicOps 1
#else
#include <PatchSideDataBasicOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataBasicOps 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class... Args>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataBasicOps)
using SAMRAIPatchSideDataBasicOps = SAMRAI::math::PatchSideDataBasicOps<Args...>;
#else
using SAMRAIPatchSideDataBasicOps = SAMRAI::math::PatchSideDataBasicOps<NDIM, Args...>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataBasicOps

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchSideDataBasicOps
