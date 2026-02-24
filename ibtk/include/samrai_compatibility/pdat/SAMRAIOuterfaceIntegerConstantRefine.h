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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIOuterfaceIntegerConstantRefine
#define included_IBTK_samrai_compatibility_pdat_SAMRAIOuterfaceIntegerConstantRefine

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/pdat/OuterfaceIntegerConstantRefine.h")
#include <SAMRAI/pdat/OuterfaceIntegerConstantRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceIntegerConstantRefine 1
#else
#include <OuterfaceIntegerConstantRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceIntegerConstantRefine 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceIntegerConstantRefine)
using SAMRAIOuterfaceIntegerConstantRefine = SAMRAI::pdat::OuterfaceIntegerConstantRefine;
#else
using SAMRAIOuterfaceIntegerConstantRefine = SAMRAI::pdat::OuterfaceIntegerConstantRefine<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_OuterfaceIntegerConstantRefine

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIOuterfaceIntegerConstantRefine
