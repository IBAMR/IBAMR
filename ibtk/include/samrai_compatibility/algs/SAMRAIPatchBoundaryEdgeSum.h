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

#ifndef included_IBTK_samrai_compatibility_algs_SAMRAIPatchBoundaryEdgeSum
#define included_IBTK_samrai_compatibility_algs_SAMRAIPatchBoundaryEdgeSum

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/algs/PatchBoundaryEdgeSum.h")
#include <SAMRAI/algs/PatchBoundaryEdgeSum.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchBoundaryEdgeSum 1
#else
#include <PatchBoundaryEdgeSum.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchBoundaryEdgeSum 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchBoundaryEdgeSum)
using SAMRAIPatchBoundaryEdgeSum = SAMRAI::algs::PatchBoundaryEdgeSum;
#else
using SAMRAIPatchBoundaryEdgeSum = SAMRAI::algs::PatchBoundaryEdgeSum<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchBoundaryEdgeSum

#endif // #ifndef included_IBTK_samrai_compatibility_algs_SAMRAIPatchBoundaryEdgeSum
