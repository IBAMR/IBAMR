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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchFaceDataOpsInteger
#define included_IBTK_samrai_compatibility_math_SAMRAIPatchFaceDataOpsInteger

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/math/PatchFaceDataOpsInteger.h>)
#include <SAMRAI/math/PatchFaceDataOpsInteger.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchFaceDataOpsInteger 1
#else
#include <PatchFaceDataOpsInteger.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchFaceDataOpsInteger 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchFaceDataOpsInteger)
using SAMRAIPatchFaceDataOpsInteger = SAMRAI::math::PatchFaceDataOpsInteger;
#else
using SAMRAIPatchFaceDataOpsInteger = SAMRAI::math::PatchFaceDataOpsInteger<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchFaceDataOpsInteger

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchFaceDataOpsInteger
