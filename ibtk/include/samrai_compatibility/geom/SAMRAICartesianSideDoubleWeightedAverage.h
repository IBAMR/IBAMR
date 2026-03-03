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

#ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianSideDoubleWeightedAverage
#define included_IBTK_samrai_compatibility_geom_SAMRAICartesianSideDoubleWeightedAverage

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/geom/CartesianSideDoubleWeightedAverage.h>)
#include <SAMRAI/geom/CartesianSideDoubleWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianSideDoubleWeightedAverage 1
#else
#include <CartesianSideDoubleWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianSideDoubleWeightedAverage 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianSideDoubleWeightedAverage)
using SAMRAICartesianSideDoubleWeightedAverage = SAMRAI::geom::CartesianSideDoubleWeightedAverage;
#else
using SAMRAICartesianSideDoubleWeightedAverage = SAMRAI::geom::CartesianSideDoubleWeightedAverage<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianSideDoubleWeightedAverage

#endif // #ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianSideDoubleWeightedAverage
