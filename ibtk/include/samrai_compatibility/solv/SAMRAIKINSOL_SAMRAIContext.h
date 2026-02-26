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

#ifndef included_IBTK_samrai_compatibility_solv_SAMRAIKINSOL_SAMRAIContext
#define included_IBTK_samrai_compatibility_solv_SAMRAIKINSOL_SAMRAIContext

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/solv/KINSOL_SAMRAIContext.h>)
#include <SAMRAI/solv/KINSOL_SAMRAIContext.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_KINSOL_SAMRAIContext 1
#else
#include <KINSOL_SAMRAIContext.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_KINSOL_SAMRAIContext 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_KINSOL_SAMRAIContext)
using SAMRAIKINSOL_SAMRAIContext = SAMRAI::solv::KINSOL_SAMRAIContext;
#else
using SAMRAIKINSOL_SAMRAIContext = SAMRAI::solv::KINSOL_SAMRAIContext<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_KINSOL_SAMRAIContext

#endif // #ifndef included_IBTK_samrai_compatibility_solv_SAMRAIKINSOL_SAMRAIContext
