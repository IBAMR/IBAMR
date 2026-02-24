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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIOuteredgeGeometry
#define included_IBTK_samrai_compatibility_pdat_SAMRAIOuteredgeGeometry

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/pdat/OuteredgeGeometry.h")
#include <SAMRAI/pdat/OuteredgeGeometry.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuteredgeGeometry 1
#else
#include <OuteredgeGeometry.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuteredgeGeometry 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuteredgeGeometry)
using SAMRAIOuteredgeGeometry = SAMRAI::pdat::OuteredgeGeometry;
#else
using SAMRAIOuteredgeGeometry = SAMRAI::pdat::OuteredgeGeometry<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuteredgeGeometry

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIOuteredgeGeometry
