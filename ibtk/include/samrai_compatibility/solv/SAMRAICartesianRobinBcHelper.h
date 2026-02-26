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

#ifndef included_IBTK_samrai_compatibility_solv_SAMRAICartesianRobinBcHelper
#define included_IBTK_samrai_compatibility_solv_SAMRAICartesianRobinBcHelper

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/solv/CartesianRobinBcHelper.h>)
#include <SAMRAI/solv/CartesianRobinBcHelper.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianRobinBcHelper 1
#else
#include <CartesianRobinBcHelper.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianRobinBcHelper 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianRobinBcHelper)
using SAMRAICartesianRobinBcHelper = SAMRAI::solv::CartesianRobinBcHelper;
#else
using SAMRAICartesianRobinBcHelper = SAMRAI::solv::CartesianRobinBcHelper<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianRobinBcHelper

#endif // #ifndef included_IBTK_samrai_compatibility_solv_SAMRAICartesianRobinBcHelper
