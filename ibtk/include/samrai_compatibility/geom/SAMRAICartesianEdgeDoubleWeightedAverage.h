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

#ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianEdgeDoubleWeightedAverage
#define included_IBTK_samrai_compatibility_geom_SAMRAICartesianEdgeDoubleWeightedAverage

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/geom/CartesianEdgeDoubleWeightedAverage.h>)
#include <SAMRAI/geom/CartesianEdgeDoubleWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeDoubleWeightedAverage 1
#else
#include <CartesianEdgeDoubleWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeDoubleWeightedAverage 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeDoubleWeightedAverage)
using SAMRAICartesianEdgeDoubleWeightedAverage = SAMRAI::geom::CartesianEdgeDoubleWeightedAverage;
#else
using SAMRAICartesianEdgeDoubleWeightedAverage = SAMRAI::geom::CartesianEdgeDoubleWeightedAverage<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeDoubleWeightedAverage

#endif // #ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianEdgeDoubleWeightedAverage
