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

#ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianOuterfaceFloatWeightedAverage
#define included_IBTK_samrai_compatibility_geom_SAMRAICartesianOuterfaceFloatWeightedAverage

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/geom/CartesianOuterfaceFloatWeightedAverage.h>)
#include <SAMRAI/geom/CartesianOuterfaceFloatWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianOuterfaceFloatWeightedAverage 1
#else
#include <CartesianOuterfaceFloatWeightedAverage.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianOuterfaceFloatWeightedAverage 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianOuterfaceFloatWeightedAverage)
using SAMRAICartesianOuterfaceFloatWeightedAverage = SAMRAI::geom::CartesianOuterfaceFloatWeightedAverage;
#else
using SAMRAICartesianOuterfaceFloatWeightedAverage = SAMRAI::geom::CartesianOuterfaceFloatWeightedAverage<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianOuterfaceFloatWeightedAverage

#endif // #ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianOuterfaceFloatWeightedAverage
