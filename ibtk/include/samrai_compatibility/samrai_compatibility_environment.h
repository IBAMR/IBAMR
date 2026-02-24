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

#ifndef included_IBTK_samrai_compatibility_environment
#define included_IBTK_samrai_compatibility_environment

#include <ibtk/config.h>

#include <samrai_compatibility/hier/SAMRAIBasePatchHierarchy.h>
#include <samrai_compatibility/hier/SAMRAIBasePatchLevel.h>
#include <samrai_compatibility/hier/SAMRAIBoxList.h>
#include <samrai_compatibility/tbox/SAMRAIPointer.h>

namespace IBTK
{
namespace SAMRAICompatibility
{
}
} // namespace IBTK

// SAMRAI 4 drops or renames several legacy SAMRAI 2 symbols used throughout
// IBAMR. Reintroduce a minimal compatibility surface in the original SAMRAI
// namespaces so old call sites compile with fewer edits.
#if defined(SAMRAI_VERSION_MAJOR) && (SAMRAI_VERSION_MAJOR < 100)

namespace SAMRAI
{
namespace tbox
{
template <class T>
using Pointer = IBTK::SAMRAICompatibility::SAMRAIPointer<T>;

template <class T>
using ConstPointer = IBTK::SAMRAICompatibility::SAMRAIPointer<const T>;
} // namespace tbox

namespace hier
{
template <int DIM = NDIM>
using BasePatchHierarchy = IBTK::SAMRAICompatibility::SAMRAIBasePatchHierarchy<DIM>;

template <int DIM = NDIM>
using BasePatchLevel = IBTK::SAMRAICompatibility::SAMRAIBasePatchLevel<DIM>;

template <int DIM = NDIM>
using BoxList = IBTK::SAMRAICompatibility::SAMRAIBoxList<DIM>;
} // namespace hier
} // namespace SAMRAI

#endif

#endif // #ifndef included_IBTK_samrai_compatibility_environment
