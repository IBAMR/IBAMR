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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataMiscellaneousOpsReal
#define included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataMiscellaneousOpsReal

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/math/PatchCellDataMiscellaneousOpsReal.h>)
#include <SAMRAI/math/PatchCellDataMiscellaneousOpsReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataMiscellaneousOpsReal 1
#else
#include <PatchCellDataMiscellaneousOpsReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataMiscellaneousOpsReal 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class TYPE>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataMiscellaneousOpsReal)
using SAMRAIPatchCellDataMiscellaneousOpsReal = SAMRAI::math::PatchCellDataMiscellaneousOpsReal<TYPE>;
#else
using SAMRAIPatchCellDataMiscellaneousOpsReal = SAMRAI::math::PatchCellDataMiscellaneousOpsReal<NDIM, TYPE>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchCellDataMiscellaneousOpsReal

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIPatchCellDataMiscellaneousOpsReal
