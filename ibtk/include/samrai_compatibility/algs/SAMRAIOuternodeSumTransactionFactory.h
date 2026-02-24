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

#ifndef included_IBTK_samrai_compatibility_algs_SAMRAIOuternodeSumTransactionFactory
#define included_IBTK_samrai_compatibility_algs_SAMRAIOuternodeSumTransactionFactory

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/algs/OuternodeSumTransactionFactory.h>)
#include <SAMRAI/algs/OuternodeSumTransactionFactory.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuternodeSumTransactionFactory 1
#else
#include <OuternodeSumTransactionFactory.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuternodeSumTransactionFactory 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuternodeSumTransactionFactory)
using SAMRAIOuternodeSumTransactionFactory = SAMRAI::algs::OuternodeSumTransactionFactory;
#else
using SAMRAIOuternodeSumTransactionFactory = SAMRAI::algs::OuternodeSumTransactionFactory<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuternodeSumTransactionFactory

#endif // #ifndef included_IBTK_samrai_compatibility_algs_SAMRAIOuternodeSumTransactionFactory
