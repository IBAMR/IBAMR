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

#ifndef included_IBTK_samrai_compatibility_algs_SAMRAITimeRefinementIntegrator
#define included_IBTK_samrai_compatibility_algs_SAMRAITimeRefinementIntegrator

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/algs/TimeRefinementIntegrator.h")
#include <SAMRAI/algs/TimeRefinementIntegrator.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_TimeRefinementIntegrator 1
#else
#include <TimeRefinementIntegrator.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_TimeRefinementIntegrator 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_TimeRefinementIntegrator)
using SAMRAITimeRefinementIntegrator = SAMRAI::algs::TimeRefinementIntegrator;
#else
using SAMRAITimeRefinementIntegrator = SAMRAI::algs::TimeRefinementIntegrator<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_TimeRefinementIntegrator

#endif // #ifndef included_IBTK_samrai_compatibility_algs_SAMRAITimeRefinementIntegrator
