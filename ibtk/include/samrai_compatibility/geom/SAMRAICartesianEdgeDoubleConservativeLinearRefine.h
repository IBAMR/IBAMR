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

#ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianEdgeDoubleConservativeLinearRefine
#define included_IBTK_samrai_compatibility_geom_SAMRAICartesianEdgeDoubleConservativeLinearRefine

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/geom/CartesianEdgeDoubleConservativeLinearRefine.h>)
#include <SAMRAI/geom/CartesianEdgeDoubleConservativeLinearRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeDoubleConservativeLinearRefine 1
#else
#include <CartesianEdgeDoubleConservativeLinearRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeDoubleConservativeLinearRefine 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeDoubleConservativeLinearRefine)
using SAMRAICartesianEdgeDoubleConservativeLinearRefine = SAMRAI::geom::CartesianEdgeDoubleConservativeLinearRefine;
#else
using SAMRAICartesianEdgeDoubleConservativeLinearRefine =
    SAMRAI::geom::CartesianEdgeDoubleConservativeLinearRefine<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianEdgeDoubleConservativeLinearRefine

#endif // #ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianEdgeDoubleConservativeLinearRefine
