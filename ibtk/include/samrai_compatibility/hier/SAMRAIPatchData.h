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

#ifndef included_IBTK_samrai_compatibility_hier_SAMRAIPatchData
#define included_IBTK_samrai_compatibility_hier_SAMRAIPatchData

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/hier/PatchData.h")
#include <SAMRAI/hier/PatchData.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchData 1
#else
#include <PatchData.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchData 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchData)
using SAMRAIPatchData = SAMRAI::hier::PatchData;
#else
using SAMRAIPatchData = SAMRAI::hier::PatchData<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_PatchData

#endif // #ifndef included_IBTK_samrai_compatibility_hier_SAMRAIPatchData
