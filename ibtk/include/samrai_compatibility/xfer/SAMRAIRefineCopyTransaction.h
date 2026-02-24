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

#ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefineCopyTransaction
#define included_IBTK_samrai_compatibility_xfer_SAMRAIRefineCopyTransaction

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/xfer/RefineCopyTransaction.h>)
#include <SAMRAI/xfer/RefineCopyTransaction.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineCopyTransaction 1
#else
#include <RefineCopyTransaction.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineCopyTransaction 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineCopyTransaction)
using SAMRAIRefineCopyTransaction = SAMRAI::xfer::RefineCopyTransaction;
#else
using SAMRAIRefineCopyTransaction = SAMRAI::xfer::RefineCopyTransaction<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineCopyTransaction

#endif // #ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefineCopyTransaction
