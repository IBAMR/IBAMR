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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIArrayDataOperationUtilities
#define included_IBTK_samrai_compatibility_pdat_SAMRAIArrayDataOperationUtilities

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/ArrayDataOperationUtilities.h>)
#include <SAMRAI/pdat/ArrayDataOperationUtilities.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataOperationUtilities 1
#else
#include <ArrayDataOperationUtilities.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataOperationUtilities 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class TYPE, class OP>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataOperationUtilities)
using SAMRAIArrayDataOperationUtilities = SAMRAI::pdat::ArrayDataOperationUtilities<TYPE, OP>;
#else
using SAMRAIArrayDataOperationUtilities = SAMRAI::pdat::ArrayDataOperationUtilities<NDIM, TYPE, OP>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataOperationUtilities

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIArrayDataOperationUtilities
