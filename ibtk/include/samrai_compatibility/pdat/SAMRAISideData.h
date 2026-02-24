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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAISideData
#define included_IBTK_samrai_compatibility_pdat_SAMRAISideData

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/SideData.h>)
#include <SAMRAI/pdat/SideData.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideData 1
#else
#include <SideData.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideData 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class T>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideData)
using SAMRAISideData = SAMRAI::pdat::SideData<T>;
#else
using SAMRAISideData = SAMRAI::pdat::SideData<NDIM, T>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideData

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAISideData
