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

#ifndef included_IBTK_samrai_compatibility_solv_SAMRAISAMRAIVectorReal
#define included_IBTK_samrai_compatibility_solv_SAMRAISAMRAIVectorReal

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/solv/SAMRAIVectorReal.h>)
#include <SAMRAI/solv/SAMRAIVectorReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SAMRAIVECTORREAL 1
#else
#include <SAMRAIVectorReal.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SAMRAIVECTORREAL 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class T>
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SAMRAIVECTORREAL)
using SAMRAISAMRAIVectorReal = SAMRAI::solv::SAMRAIVectorReal<T>;
#else
using SAMRAISAMRAIVectorReal = SAMRAI::solv::SAMRAIVectorReal<NDIM, T>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_SAMRAIVECTORREAL

#endif // #ifndef included_IBTK_samrai_compatibility_solv_SAMRAISAMRAIVectorReal
