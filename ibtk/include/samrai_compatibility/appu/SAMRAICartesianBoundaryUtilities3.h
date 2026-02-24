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

#ifndef included_IBTK_samrai_compatibility_appu_SAMRAICartesianBoundaryUtilities3
#define included_IBTK_samrai_compatibility_appu_SAMRAICartesianBoundaryUtilities3

#include "samrai_compatibility/samrai_compatibility_detect.h"

#if IBTK_SAMRAI_HAS_INCLUDE(<SAMRAI/appu/CartesianBoundaryUtilities3.h>)
#include <SAMRAI/appu/CartesianBoundaryUtilities3.h>
#else
#include <CartesianBoundaryUtilities3.h>
#endif

namespace IBTK
{
namespace SAMRAICompatibility
{
using SAMRAICartesianBoundaryUtilities3 = SAMRAI::appu::CartesianBoundaryUtilities3;
} // namespace SAMRAICompatibility
} // namespace IBTK

#endif // #ifndef included_IBTK_samrai_compatibility_appu_SAMRAICartesianBoundaryUtilities3
