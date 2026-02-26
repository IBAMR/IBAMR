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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchFaceDataBasicOps
#define included_IBTK_samrai_compatibility_math_SAMRAIPatchFaceDataBasicOps

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/math/PatchFaceDataBasicOps.h>)
#include <SAMRAI/math/PatchFaceDataBasicOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchFaceDataBasicOps 1
#else
#include <PatchFaceDataBasicOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchFaceDataBasicOps 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class... Args>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchFaceDataBasicOps)
using SAMRAIPatchFaceDataBasicOps = SAMRAI::math::PatchFaceDataBasicOps<Args...>;
#else
using SAMRAIPatchFaceDataBasicOps = SAMRAI::math::PatchFaceDataBasicOps<NDIM, Args...>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchFaceDataBasicOps

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchFaceDataBasicOps
