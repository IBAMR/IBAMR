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

#ifndef included_IBTK_samrai_compatibility_xfer_SAMRAICoarsenAlgorithm
#define included_IBTK_samrai_compatibility_xfer_SAMRAICoarsenAlgorithm

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/xfer/CoarsenAlgorithm.h>)
#include <SAMRAI/xfer/CoarsenAlgorithm.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CoarsenAlgorithm 1
#else
#include <CoarsenAlgorithm.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CoarsenAlgorithm 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CoarsenAlgorithm)
using SAMRAICoarsenAlgorithm = SAMRAI::xfer::CoarsenAlgorithm;
#else
using SAMRAICoarsenAlgorithm = SAMRAI::xfer::CoarsenAlgorithm<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CoarsenAlgorithm

#endif // #ifndef included_IBTK_samrai_compatibility_xfer_SAMRAICoarsenAlgorithm
