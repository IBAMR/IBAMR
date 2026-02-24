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

#ifndef included_IBTK_samrai_compatibility_xfer_SAMRAICoarsenTransactionFactory
#define included_IBTK_samrai_compatibility_xfer_SAMRAICoarsenTransactionFactory

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/xfer/CoarsenTransactionFactory.h")
#include <SAMRAI/xfer/CoarsenTransactionFactory.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CoarsenTransactionFactory 1
#else
#include <CoarsenTransactionFactory.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CoarsenTransactionFactory 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CoarsenTransactionFactory)
using SAMRAICoarsenTransactionFactory = SAMRAI::xfer::CoarsenTransactionFactory;
#else
using SAMRAICoarsenTransactionFactory = SAMRAI::xfer::CoarsenTransactionFactory<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CoarsenTransactionFactory

#endif // #ifndef included_IBTK_samrai_compatibility_xfer_SAMRAICoarsenTransactionFactory
