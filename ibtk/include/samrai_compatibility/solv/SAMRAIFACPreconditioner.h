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

#ifndef included_IBTK_samrai_compatibility_solv_SAMRAIFACPreconditioner
#define included_IBTK_samrai_compatibility_solv_SAMRAIFACPreconditioner

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/solv/FACPreconditioner.h")
#include <SAMRAI/solv/FACPreconditioner.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_FACPreconditioner 1
#else
#include <FACPreconditioner.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_FACPreconditioner 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_FACPreconditioner)
using SAMRAIFACPreconditioner = SAMRAI::solv::FACPreconditioner;
#else
using SAMRAIFACPreconditioner = SAMRAI::solv::FACPreconditioner<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_FACPreconditioner

#endif // #ifndef included_IBTK_samrai_compatibility_solv_SAMRAIFACPreconditioner
