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

#ifndef included_IBTK_samrai_compatibility_algs_SAMRAIOuteredgeSumTransaction
#define included_IBTK_samrai_compatibility_algs_SAMRAIOuteredgeSumTransaction

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/algs/OuteredgeSumTransaction.h>)
#include <SAMRAI/algs/OuteredgeSumTransaction.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuteredgeSumTransaction 1
#else
#include <OuteredgeSumTransaction.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuteredgeSumTransaction 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuteredgeSumTransaction)
using SAMRAIOuteredgeSumTransaction = SAMRAI::algs::OuteredgeSumTransaction;
#else
using SAMRAIOuteredgeSumTransaction = SAMRAI::algs::OuteredgeSumTransaction<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuteredgeSumTransaction

#endif // #ifndef included_IBTK_samrai_compatibility_algs_SAMRAIOuteredgeSumTransaction
