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

#ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianCellDoubleWeightedAverage
#define included_IBTK_samrai_compatibility_geom_SAMRAICartesianCellDoubleWeightedAverage

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/geom/CartesianCellDoubleWeightedAverage.h>)
#include <SAMRAI/geom/CartesianCellDoubleWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianCellDoubleWeightedAverage 1
#else
#include <CartesianCellDoubleWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianCellDoubleWeightedAverage 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianCellDoubleWeightedAverage)
using SAMRAICartesianCellDoubleWeightedAverage = SAMRAI::geom::CartesianCellDoubleWeightedAverage;
#else
using SAMRAICartesianCellDoubleWeightedAverage = SAMRAI::geom::CartesianCellDoubleWeightedAverage<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianCellDoubleWeightedAverage

#endif // #ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianCellDoubleWeightedAverage
