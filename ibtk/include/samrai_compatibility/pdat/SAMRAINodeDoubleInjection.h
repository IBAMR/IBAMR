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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAINodeDoubleInjection
#define included_IBTK_samrai_compatibility_pdat_SAMRAINodeDoubleInjection

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/NodeDoubleInjection.h>)
#include <SAMRAI/pdat/NodeDoubleInjection.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeDoubleInjection 1
#else
#include <NodeDoubleInjection.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeDoubleInjection 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeDoubleInjection)
using SAMRAINodeDoubleInjection = SAMRAI::pdat::NodeDoubleInjection;
#else
using SAMRAINodeDoubleInjection = SAMRAI::pdat::NodeDoubleInjection<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeDoubleInjection

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAINodeDoubleInjection
