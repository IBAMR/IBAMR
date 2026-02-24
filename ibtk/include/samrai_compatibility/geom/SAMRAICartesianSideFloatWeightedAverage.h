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

#ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianSideFloatWeightedAverage
#define included_IBTK_samrai_compatibility_geom_SAMRAICartesianSideFloatWeightedAverage

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/geom/CartesianSideFloatWeightedAverage.h")
#include <SAMRAI/geom/CartesianSideFloatWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianSideFloatWeightedAverage 1
#else
#include <CartesianSideFloatWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianSideFloatWeightedAverage 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianSideFloatWeightedAverage)
using SAMRAICartesianSideFloatWeightedAverage = SAMRAI::geom::CartesianSideFloatWeightedAverage;
#else
using SAMRAICartesianSideFloatWeightedAverage = SAMRAI::geom::CartesianSideFloatWeightedAverage<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianSideFloatWeightedAverage

#endif // #ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianSideFloatWeightedAverage
