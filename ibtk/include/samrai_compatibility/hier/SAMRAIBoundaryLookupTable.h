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

#ifndef included_IBTK_samrai_compatibility_hier_SAMRAIBoundaryLookupTable
#define included_IBTK_samrai_compatibility_hier_SAMRAIBoundaryLookupTable

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/hier/BoundaryLookupTable.h")
#include <SAMRAI/hier/BoundaryLookupTable.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoundaryLookupTable 1
#else
#include <BoundaryLookupTable.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoundaryLookupTable 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoundaryLookupTable)
using SAMRAIBoundaryLookupTable = SAMRAI::hier::BoundaryLookupTable;
#else
using SAMRAIBoundaryLookupTable = SAMRAI::hier::BoundaryLookupTable<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoundaryLookupTable

#endif // #ifndef included_IBTK_samrai_compatibility_hier_SAMRAIBoundaryLookupTable
