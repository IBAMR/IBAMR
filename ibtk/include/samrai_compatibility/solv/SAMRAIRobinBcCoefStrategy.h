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

#ifndef included_IBTK_samrai_compatibility_solv_SAMRAIRobinBcCoefStrategy
#define included_IBTK_samrai_compatibility_solv_SAMRAIRobinBcCoefStrategy

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/solv/RobinBcCoefStrategy.h>)
#include <SAMRAI/solv/RobinBcCoefStrategy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RobinBcCoefStrategy 1
#else
#include <RobinBcCoefStrategy.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RobinBcCoefStrategy 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RobinBcCoefStrategy)
using SAMRAIRobinBcCoefStrategy = SAMRAI::solv::RobinBcCoefStrategy;
#else
using SAMRAIRobinBcCoefStrategy = SAMRAI::solv::RobinBcCoefStrategy<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_RobinBcCoefStrategy

#endif // #ifndef included_IBTK_samrai_compatibility_solv_SAMRAIRobinBcCoefStrategy
