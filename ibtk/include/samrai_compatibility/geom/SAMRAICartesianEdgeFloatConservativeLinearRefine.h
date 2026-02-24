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

#ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianEdgeFloatConservativeLinearRefine
#define included_IBTK_samrai_compatibility_geom_SAMRAICartesianEdgeFloatConservativeLinearRefine

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/geom/CartesianEdgeFloatConservativeLinearRefine.h")
#include <SAMRAI/geom/CartesianEdgeFloatConservativeLinearRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeFloatConservativeLinearRefine 1
#else
#include <CartesianEdgeFloatConservativeLinearRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeFloatConservativeLinearRefine 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeFloatConservativeLinearRefine)
using SAMRAICartesianEdgeFloatConservativeLinearRefine = SAMRAI::geom::CartesianEdgeFloatConservativeLinearRefine;
#else
using SAMRAICartesianEdgeFloatConservativeLinearRefine = SAMRAI::geom::CartesianEdgeFloatConservativeLinearRefine<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeFloatConservativeLinearRefine

#endif // #ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianEdgeFloatConservativeLinearRefine
