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

#ifndef included_IBTK_samrai_compatibility_solv_SAMRAICellPoissonFACOps
#define included_IBTK_samrai_compatibility_solv_SAMRAICellPoissonFACOps

#include <ibtk/config.h>

#include <samrai_compatibility/samrai_compatibility_detect.h>

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/solv/CellPoissonFACOps.h>)
#include <SAMRAI/solv/CellPoissonFACOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellPoissonFACOps 1
#else
#include <CellPoissonFACOps.h>
#define IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellPoissonFACOps 0
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
#if (IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellPoissonFACOps)
using SAMRAICellPoissonFACOps = SAMRAI::solv::CellPoissonFACOps;
#else
using SAMRAICellPoissonFACOps = SAMRAI::solv::CellPoissonFACOps<NDIM>;
#endif
} // namespace SAMRAICompatibility
} // namespace IBTK

#undef IBTK_SAMRAI_COMPAT_UNTEMPLATED_DIM_CellPoissonFACOps

#endif // #ifndef included_IBTK_samrai_compatibility_solv_SAMRAICellPoissonFACOps
