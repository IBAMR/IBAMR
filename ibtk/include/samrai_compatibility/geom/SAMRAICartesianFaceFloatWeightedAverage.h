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

#ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianFaceFloatWeightedAverage
#define included_IBTK_samrai_compatibility_geom_SAMRAICartesianFaceFloatWeightedAverage

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/geom/CartesianFaceFloatWeightedAverage.h>)
#include <SAMRAI/geom/CartesianFaceFloatWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianFaceFloatWeightedAverage 1
#else
#include <CartesianFaceFloatWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianFaceFloatWeightedAverage 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianFaceFloatWeightedAverage)
using SAMRAICartesianFaceFloatWeightedAverage = SAMRAI::geom::CartesianFaceFloatWeightedAverage;
#else
using SAMRAICartesianFaceFloatWeightedAverage = SAMRAI::geom::CartesianFaceFloatWeightedAverage<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianFaceFloatWeightedAverage

#endif // #ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianFaceFloatWeightedAverage
