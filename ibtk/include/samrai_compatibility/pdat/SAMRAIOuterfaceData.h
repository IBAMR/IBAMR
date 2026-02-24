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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIOuterfaceData
#define included_IBTK_samrai_compatibility_pdat_SAMRAIOuterfaceData

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/OuterfaceData.h>)
#include <SAMRAI/pdat/OuterfaceData.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceData 1
#else
#include <OuterfaceData.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceData 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class... Args>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceData)
using SAMRAIOuterfaceData = SAMRAI::pdat::OuterfaceData<Args...>;
#else
using SAMRAIOuterfaceData = SAMRAI::pdat::OuterfaceData<NDIM, Args...>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceData

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIOuterfaceData
