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

#ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefinePatchStrategy
#define included_IBTK_samrai_compatibility_xfer_SAMRAIRefinePatchStrategy

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/xfer/RefinePatchStrategy.h>)
#include <SAMRAI/xfer/RefinePatchStrategy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefinePatchStrategy 1
#else
#include <RefinePatchStrategy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefinePatchStrategy 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefinePatchStrategy)
using SAMRAIRefinePatchStrategy = SAMRAI::xfer::RefinePatchStrategy;
#else
using SAMRAIRefinePatchStrategy = SAMRAI::xfer::RefinePatchStrategy<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefinePatchStrategy

#endif // #ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefinePatchStrategy
