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

#ifndef included_IBTK_samrai_compatibility_algs_SAMRAITimeRefinementLevelStrategy
#define included_IBTK_samrai_compatibility_algs_SAMRAITimeRefinementLevelStrategy

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/algs/TimeRefinementLevelStrategy.h>)
#include <SAMRAI/algs/TimeRefinementLevelStrategy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_TimeRefinementLevelStrategy 1
#else
#include <TimeRefinementLevelStrategy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_TimeRefinementLevelStrategy 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_TimeRefinementLevelStrategy)
using SAMRAITimeRefinementLevelStrategy = SAMRAI::algs::TimeRefinementLevelStrategy;
#else
using SAMRAITimeRefinementLevelStrategy = SAMRAI::algs::TimeRefinementLevelStrategy<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_TimeRefinementLevelStrategy

#endif // #ifndef included_IBTK_samrai_compatibility_algs_SAMRAITimeRefinementLevelStrategy
