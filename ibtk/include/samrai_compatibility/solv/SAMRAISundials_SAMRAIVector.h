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

#ifndef included_IBTK_samrai_compatibility_solv_SAMRAISundials_SAMRAIVector
#define included_IBTK_samrai_compatibility_solv_SAMRAISundials_SAMRAIVector

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/solv/Sundials_SAMRAIVector.h>)
#include <SAMRAI/solv/Sundials_SAMRAIVector.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_Sundials_SAMRAIVector 1
#else
#include <Sundials_SAMRAIVector.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_Sundials_SAMRAIVector 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_Sundials_SAMRAIVector)
using SAMRAISundials_SAMRAIVector = SAMRAI::solv::Sundials_SAMRAIVector;
#else
using SAMRAISundials_SAMRAIVector = SAMRAI::solv::Sundials_SAMRAIVector<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_Sundials_SAMRAIVector

#endif // #ifndef included_IBTK_samrai_compatibility_solv_SAMRAISundials_SAMRAIVector
