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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIIndexDataFactory
#define included_IBTK_samrai_compatibility_pdat_SAMRAIIndexDataFactory

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/IndexDataFactory.h>)
#include <SAMRAI/pdat/IndexDataFactory.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_IndexDataFactory 1
#else
#include <IndexDataFactory.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_IndexDataFactory 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class T, class BoxGeometry>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_IndexDataFactory)
using SAMRAIIndexDataFactory = SAMRAI::pdat::IndexDataFactory<T, BoxGeometry>;
#else
using SAMRAIIndexDataFactory = SAMRAI::pdat::IndexDataFactory<NDIM, T, BoxGeometry>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_IndexDataFactory

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIIndexDataFactory
