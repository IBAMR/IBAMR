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

#ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefineTimeTransaction
#define included_IBTK_samrai_compatibility_xfer_SAMRAIRefineTimeTransaction

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/xfer/RefineTimeTransaction.h>)
#include <SAMRAI/xfer/RefineTimeTransaction.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineTimeTransaction 1
#else
#include <RefineTimeTransaction.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineTimeTransaction 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineTimeTransaction)
using SAMRAIRefineTimeTransaction = SAMRAI::xfer::RefineTimeTransaction;
#else
using SAMRAIRefineTimeTransaction = SAMRAI::xfer::RefineTimeTransaction<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RefineTimeTransaction

#endif // #ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIRefineTimeTransaction
