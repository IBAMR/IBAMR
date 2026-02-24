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

#ifndef included_IBTK_samrai_compatibility_pdat_SAMRAICopyOperation
#define included_IBTK_samrai_compatibility_pdat_SAMRAICopyOperation

#include <ibtk/config.h>

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/pdat/CopyOperation.h>)
#include <SAMRAI/pdat/CopyOperation.h>
#else
#include <CopyOperation.h>
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
template <class TYPE>
using SAMRAICopyOperation = SAMRAI::pdat::CopyOperation<TYPE>;
} // namespace SAMRAICompatibility
} // namespace IBTK

#endif // #ifndef included_IBTK_samrai_compatibility_pdat_SAMRAICopyOperation
