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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAINodeData
#define included_IBTK_samrai_compatibility_pdat_SAMRAINodeData

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/NodeData.h>)
#include <SAMRAI/pdat/NodeData.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeData 1
#else
#include <NodeData.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeData 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class T>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeData)
using SAMRAINodeData = SAMRAI::pdat::NodeData<T>;
#else
using SAMRAINodeData = SAMRAI::pdat::NodeData<NDIM, T>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_NodeData

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAINodeData
