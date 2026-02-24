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

#ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefineAlgorithm
#define included_IBTK_samrai_compatibility_xfer_SAMRAIRefineAlgorithm

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/xfer/RefineAlgorithm.h")
#include <SAMRAI/xfer/RefineAlgorithm.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineAlgorithm 1
#else
#include <RefineAlgorithm.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineAlgorithm 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineAlgorithm)
using SAMRAIRefineAlgorithm = SAMRAI::xfer::RefineAlgorithm;
#else
using SAMRAIRefineAlgorithm = SAMRAI::xfer::RefineAlgorithm<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineAlgorithm

#endif // #ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefineAlgorithm
