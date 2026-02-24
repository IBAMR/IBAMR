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

#ifndef included_IBTK_samrai_compatibility_mesh_SAMRAIBoxGeneratorStrategy
#define included_IBTK_samrai_compatibility_mesh_SAMRAIBoxGeneratorStrategy

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/mesh/BoxGeneratorStrategy.h>)
#include <SAMRAI/mesh/BoxGeneratorStrategy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxGeneratorStrategy 1
#else
#include <BoxGeneratorStrategy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxGeneratorStrategy 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxGeneratorStrategy)
using SAMRAIBoxGeneratorStrategy = SAMRAI::mesh::BoxGeneratorStrategy;
#else
using SAMRAIBoxGeneratorStrategy = SAMRAI::mesh::BoxGeneratorStrategy<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxGeneratorStrategy

#endif // #ifndef included_IBTK_samrai_compatibility_mesh_SAMRAIBoxGeneratorStrategy
