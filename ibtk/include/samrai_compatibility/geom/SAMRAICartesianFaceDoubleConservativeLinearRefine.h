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

#ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianFaceDoubleConservativeLinearRefine
#define included_IBTK_samrai_compatibility_geom_SAMRAICartesianFaceDoubleConservativeLinearRefine

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/geom/CartesianFaceDoubleConservativeLinearRefine.h>)
#include <SAMRAI/geom/CartesianFaceDoubleConservativeLinearRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianFaceDoubleConservativeLinearRefine 1
#else
#include <CartesianFaceDoubleConservativeLinearRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianFaceDoubleConservativeLinearRefine 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianFaceDoubleConservativeLinearRefine)
using SAMRAICartesianFaceDoubleConservativeLinearRefine = SAMRAI::geom::CartesianFaceDoubleConservativeLinearRefine;
#else
using SAMRAICartesianFaceDoubleConservativeLinearRefine =
    SAMRAI::geom::CartesianFaceDoubleConservativeLinearRefine<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CartesianFaceDoubleConservativeLinearRefine

#endif // #ifndef included_IBTK_samrai_compatibility_geom_SAMRAICartesianFaceDoubleConservativeLinearRefine
