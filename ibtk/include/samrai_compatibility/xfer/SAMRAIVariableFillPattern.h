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

#ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIVariableFillPattern
#define included_IBTK_samrai_compatibility_xfer_SAMRAIVariableFillPattern

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/xfer/VariableFillPattern.h")
#include <SAMRAI/xfer/VariableFillPattern.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_VariableFillPattern 1
#else
#include <VariableFillPattern.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_VariableFillPattern 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_VariableFillPattern)
using SAMRAIVariableFillPattern = SAMRAI::xfer::VariableFillPattern;
#else
using SAMRAIVariableFillPattern = SAMRAI::xfer::VariableFillPattern<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_VariableFillPattern

#endif // #ifndef included_IBTK_samrai_compatibility_xfer_SAMRAIVariableFillPattern
