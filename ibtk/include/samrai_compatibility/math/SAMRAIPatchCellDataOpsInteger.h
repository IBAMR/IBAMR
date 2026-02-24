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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataOpsInteger
#define included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataOpsInteger

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/math/PatchCellDataOpsInteger.h")
#include <SAMRAI/math/PatchCellDataOpsInteger.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataOpsInteger 1
#else
#include <PatchCellDataOpsInteger.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataOpsInteger 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataOpsInteger)
using SAMRAIPatchCellDataOpsInteger = SAMRAI::math::PatchCellDataOpsInteger;
#else
using SAMRAIPatchCellDataOpsInteger = SAMRAI::math::PatchCellDataOpsInteger<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataOpsInteger

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataOpsInteger
