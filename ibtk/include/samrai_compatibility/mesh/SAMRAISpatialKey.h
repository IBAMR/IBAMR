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

#ifndef included_IBTK_samrai_compatibility_mesh_SAMRAISpatialKey
#define included_IBTK_samrai_compatibility_mesh_SAMRAISpatialKey

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/mesh/SpatialKey.h>)
#include <SAMRAI/mesh/SpatialKey.h>
#else
#include <SpatialKey.h>
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
using SAMRAISpatialKey = SAMRAI::mesh::SpatialKey;
} // namespace SAMRAICompatibility
} // namespace IBTK

#endif // #ifndef included_IBTK_samrai_compatibility_mesh_SAMRAISpatialKey
