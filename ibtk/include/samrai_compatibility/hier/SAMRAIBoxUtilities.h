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

#ifndef included_IBTK_samrai_compatibility_hier_SAMRAIBoxUtilities
#define included_IBTK_samrai_compatibility_hier_SAMRAIBoxUtilities

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/hier/BoxUtilities.h>)
#include <SAMRAI/hier/BoxUtilities.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxUtilities 1
#else
#include <BoxUtilities.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxUtilities 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxUtilities)
using SAMRAIBoxUtilities = SAMRAI::hier::BoxUtilities;
#else
using SAMRAIBoxUtilities = SAMRAI::hier::BoxUtilities<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxUtilities

#endif // #ifndef included_IBTK_samrai_compatibility_hier_SAMRAIBoxUtilities
