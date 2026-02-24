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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIFaceDoubleConstantRefine
#define included_IBTK_samrai_compatibility_pdat_SAMRAIFaceDoubleConstantRefine

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE("SAMRAI/pdat/FaceDoubleConstantRefine.h")
#include <SAMRAI/pdat/FaceDoubleConstantRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_FaceDoubleConstantRefine 1
#else
#include <FaceDoubleConstantRefine.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_FaceDoubleConstantRefine 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_FaceDoubleConstantRefine)
using SAMRAIFaceDoubleConstantRefine = SAMRAI::pdat::FaceDoubleConstantRefine;
#else
using SAMRAIFaceDoubleConstantRefine = SAMRAI::pdat::FaceDoubleConstantRefine<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_FaceDoubleConstantRefine

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAIFaceDoubleConstantRefine
