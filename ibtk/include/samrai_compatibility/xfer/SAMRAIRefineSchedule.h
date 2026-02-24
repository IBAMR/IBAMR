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

#ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefineSchedule
#define included_IBTK_samrai_compatibility_xfer_SAMRAIRefineSchedule

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/xfer/RefineSchedule.h")
#include <SAMRAI/xfer/RefineSchedule.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineSchedule 1
#else
#include <RefineSchedule.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineSchedule 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineSchedule)
using SAMRAIRefineSchedule = SAMRAI::xfer::RefineSchedule;
#else
using SAMRAIRefineSchedule = SAMRAI::xfer::RefineSchedule<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineSchedule

#endif // #ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefineSchedule
