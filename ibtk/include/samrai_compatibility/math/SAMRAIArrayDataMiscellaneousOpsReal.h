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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIArrayDataMiscellaneousOpsReal
#define included_IBTK_samrai_compatibility_math_SAMRAIArrayDataMiscellaneousOpsReal

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/math/ArrayDataMiscellaneousOpsReal.h>)
#include <SAMRAI/math/ArrayDataMiscellaneousOpsReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataMiscellaneousOpsReal 1
#else
#include <ArrayDataMiscellaneousOpsReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataMiscellaneousOpsReal 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class TYPE>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataMiscellaneousOpsReal)
using SAMRAIArrayDataMiscellaneousOpsReal = SAMRAI::math::ArrayDataMiscellaneousOpsReal<TYPE>;
#else
using SAMRAIArrayDataMiscellaneousOpsReal = SAMRAI::math::ArrayDataMiscellaneousOpsReal<NDIM, TYPE>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataMiscellaneousOpsReal

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIArrayDataMiscellaneousOpsReal
