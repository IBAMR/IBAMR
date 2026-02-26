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

#ifndef included_IBTK_samrai_compatibility_hier_SAMRAIBasePatchHierarchy
#define included_IBTK_samrai_compatibility_hier_SAMRAIBasePatchHierarchy

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/hier/PatchHierarchy.h>)
#include <SAMRAI/hier/PatchHierarchy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BasePatchHierarchy 1
#else
#include <BasePatchHierarchy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BasePatchHierarchy 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BasePatchHierarchy)
using SAMRAIBasePatchHierarchy = SAMRAI::hier::PatchHierarchy;
#else
using SAMRAIBasePatchHierarchy = SAMRAI::hier::BasePatchHierarchy<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BasePatchHierarchy

#endif // #ifndef included_IBTK_samrai_compatibility_hier_SAMRAIBasePatchHierarchy
