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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataBasicOps
#define included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataBasicOps

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/math/PatchCellDataBasicOps.h>)
#include <SAMRAI/math/PatchCellDataBasicOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataBasicOps 1
#else
#include <PatchCellDataBasicOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataBasicOps 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class... Args>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataBasicOps)
using SAMRAIPatchCellDataBasicOps = SAMRAI::math::PatchCellDataBasicOps<Args...>;
#else
using SAMRAIPatchCellDataBasicOps = SAMRAI::math::PatchCellDataBasicOps<NDIM, Args...>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataBasicOps

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataBasicOps
