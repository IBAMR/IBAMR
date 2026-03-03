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

#ifndef included_IBTK_samrai_compatibility_mesh_SAMRAIGriddingAlgorithm
#define included_IBTK_samrai_compatibility_mesh_SAMRAIGriddingAlgorithm

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/mesh/GriddingAlgorithm.h>)
#include <SAMRAI/mesh/GriddingAlgorithm.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_GriddingAlgorithm 1
#else
#include <GriddingAlgorithm.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_GriddingAlgorithm 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_GriddingAlgorithm)
using SAMRAIGriddingAlgorithm = SAMRAI::mesh::GriddingAlgorithm;
#else
using SAMRAIGriddingAlgorithm = SAMRAI::mesh::GriddingAlgorithm<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_GriddingAlgorithm

#endif // #ifndef included_IBTK_samrai_compatibility_mesh_SAMRAIGriddingAlgorithm
