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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchSideDataMiscellaneousOpsReal
#define included_IBTK_samrai_compatibility_math_SAMRAIPatchSideDataMiscellaneousOpsReal

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/math/PatchSideDataMiscellaneousOpsReal.h")
#include <SAMRAI/math/PatchSideDataMiscellaneousOpsReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataMiscellaneousOpsReal 1
#else
#include <PatchSideDataMiscellaneousOpsReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataMiscellaneousOpsReal 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class TYPE>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataMiscellaneousOpsReal)
using SAMRAIPatchSideDataMiscellaneousOpsReal = SAMRAI::math::PatchSideDataMiscellaneousOpsReal<TYPE>;
#else
using SAMRAIPatchSideDataMiscellaneousOpsReal = SAMRAI::math::PatchSideDataMiscellaneousOpsReal<NDIM, TYPE>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchSideDataMiscellaneousOpsReal

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchSideDataMiscellaneousOpsReal
