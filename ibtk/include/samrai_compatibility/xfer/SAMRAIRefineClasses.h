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

#ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefineClasses
#define included_IBTK_samrai_compatibility_xfer_SAMRAIRefineClasses

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/xfer/RefineClasses.h>)
#include <SAMRAI/xfer/RefineClasses.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineClasses 1
#else
#include <RefineClasses.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineClasses 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineClasses)
using SAMRAIRefineClasses = SAMRAI::xfer::RefineClasses;
#else
using SAMRAIRefineClasses = SAMRAI::xfer::RefineClasses<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineClasses

#endif // #ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefineClasses
