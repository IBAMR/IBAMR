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

#ifndef included_IBTK_samrai_compatibility_solv_SAMRAICellPoissonFACSolver
#define included_IBTK_samrai_compatibility_solv_SAMRAICellPoissonFACSolver

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/solv/CellPoissonFACSolver.h>)
#include <SAMRAI/solv/CellPoissonFACSolver.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellPoissonFACSolver 1
#else
#include <CellPoissonFACSolver.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellPoissonFACSolver 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellPoissonFACSolver)
using SAMRAICellPoissonFACSolver = SAMRAI::solv::CellPoissonFACSolver;
#else
using SAMRAICellPoissonFACSolver = SAMRAI::solv::CellPoissonFACSolver<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellPoissonFACSolver

#endif // #ifndef included_IBTK_samrai_compatibility_solv_SAMRAICellPoissonFACSolver
