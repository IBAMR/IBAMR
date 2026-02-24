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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAISideComplexConstantRefine
#define included_IBTK_samrai_compatibility_pdat_SAMRAISideComplexConstantRefine

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/pdat/SideComplexConstantRefine.h")
#include <SAMRAI/pdat/SideComplexConstantRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideComplexConstantRefine 1
#else
#include <SideComplexConstantRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideComplexConstantRefine 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideComplexConstantRefine)
using SAMRAISideComplexConstantRefine = SAMRAI::pdat::SideComplexConstantRefine;
#else
using SAMRAISideComplexConstantRefine = SAMRAI::pdat::SideComplexConstantRefine<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SideComplexConstantRefine

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAISideComplexConstantRefine
