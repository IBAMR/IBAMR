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

#ifndef included_IBTK_samrai_compatibility_math_SAMRAIArrayDataBasicOps
#define included_IBTK_samrai_compatibility_math_SAMRAIArrayDataBasicOps

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/math/ArrayDataBasicOps.h")
#include <SAMRAI/math/ArrayDataBasicOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataBasicOps 1
#else
#include <ArrayDataBasicOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataBasicOps 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class... Args>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataBasicOps)
using SAMRAIArrayDataBasicOps = SAMRAI::math::ArrayDataBasicOps<Args...>;
#else
using SAMRAIArrayDataBasicOps = SAMRAI::math::ArrayDataBasicOps<NDIM, Args...>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_ArrayDataBasicOps

#endif // #ifndef included_IBTK_samrai_compatibility_math_SAMRAIArrayDataBasicOps
