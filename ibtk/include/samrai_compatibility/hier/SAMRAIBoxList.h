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

#ifndef included_IBTK_samrai_compatibility_hier_SAMRAIBoxList
#define included_IBTK_samrai_compatibility_hier_SAMRAIBoxList

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/hier/BoxContainer.h")
#include <SAMRAI/hier/BoxContainer.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxList 1
#else
#include <BoxList.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxList 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxList)
using SAMRAIBoxList = SAMRAI::hier::BoxContainer;
#else
using SAMRAIBoxList = SAMRAI::hier::BoxList<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_BoxList

#endif // #ifndef included_IBTK_samrai_compatibility_hier_SAMRAIBoxList
