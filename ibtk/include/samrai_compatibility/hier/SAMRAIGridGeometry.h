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

#ifndef included_IBTK_samrai_compatibility_hier_SAMRAIGridGeometry
#define included_IBTK_samrai_compatibility_hier_SAMRAIGridGeometry

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/geom/GridGeometry.h")
#include <SAMRAI/geom/GridGeometry.h>
#define IBTK_SAMRAI_COMPAT_GRID_GEOM_IN_GEOM_NS 1
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_GridGeometry 1
#else
#include <GridGeometry.h>
#define IBTK_SAMRAI_COMPAT_GRID_GEOM_IN_GEOM_NS 0
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_GridGeometry 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_GRID_GEOM_IN_GEOM_NS)
using SAMRAIGridGeometry = SAMRAI::geom::GridGeometry;
#else
using SAMRAIGridGeometry = SAMRAI::hier::GridGeometry<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_GRID_GEOM_IN_GEOM_NS
#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_GridGeometry

#endif // #ifndef included_IBTK_samrai_compatibility_hier_SAMRAIGridGeometry
