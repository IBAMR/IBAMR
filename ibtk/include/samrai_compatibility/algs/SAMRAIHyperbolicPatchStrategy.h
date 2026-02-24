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

#ifndef included_IBTK_samrai_compatibility_algs_SAMRAIHyperbolicPatchStrategy
#define included_IBTK_samrai_compatibility_algs_SAMRAIHyperbolicPatchStrategy

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/algs/HyperbolicPatchStrategy.h>)
#include <SAMRAI/algs/HyperbolicPatchStrategy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HyperbolicPatchStrategy 1
#else
#include <HyperbolicPatchStrategy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HyperbolicPatchStrategy 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HyperbolicPatchStrategy)
using SAMRAIHyperbolicPatchStrategy = SAMRAI::algs::HyperbolicPatchStrategy;
#else
using SAMRAIHyperbolicPatchStrategy = SAMRAI::algs::HyperbolicPatchStrategy<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_HyperbolicPatchStrategy

#endif // #ifndef included_IBTK_samrai_compatibility_algs_SAMRAIHyperbolicPatchStrategy
