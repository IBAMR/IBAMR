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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchSideDataNormOpsReal
#define included_IBTK_samrai_compatibility_math_SAMRAIPatchSideDataNormOpsReal

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/math/PatchSideDataNormOpsReal.h")
#include <SAMRAI/math/PatchSideDataNormOpsReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataNormOpsReal 1
#else
#include <PatchSideDataNormOpsReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataNormOpsReal 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class... Args>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataNormOpsReal)
using SAMRAIPatchSideDataNormOpsReal = SAMRAI::math::PatchSideDataNormOpsReal<Args...>;
#else
using SAMRAIPatchSideDataNormOpsReal = SAMRAI::math::PatchSideDataNormOpsReal<NDIM, Args...>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataNormOpsReal

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchSideDataNormOpsReal
