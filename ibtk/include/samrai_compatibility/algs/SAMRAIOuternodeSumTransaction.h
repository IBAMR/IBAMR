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

#ifndef included_IBTK_samrai_compatibility_algs_SAMRAIOuternodeSumTransaction
#define included_IBTK_samrai_compatibility_algs_SAMRAIOuternodeSumTransaction

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/algs/OuternodeSumTransaction.h")
#include <SAMRAI/algs/OuternodeSumTransaction.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuternodeSumTransaction 1
#else
#include <OuternodeSumTransaction.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuternodeSumTransaction 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuternodeSumTransaction)
using SAMRAIOuternodeSumTransaction = SAMRAI::algs::OuternodeSumTransaction;
#else
using SAMRAIOuternodeSumTransaction = SAMRAI::algs::OuternodeSumTransaction<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuternodeSumTransaction

#endif // #ifndef included_IBTK_samrai_compatibility_algs_SAMRAIOuternodeSumTransaction
