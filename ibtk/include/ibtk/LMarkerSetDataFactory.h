// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_LMarkerSetDataFactory
#define included_IBTK_LMarkerSetDataFactory

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LMarker.h"
#include "ibtk/LSetDataFactory.h"

/////////////////////////////// TYPEDEFS /////////////////////////////////////

namespace IBTK
{
using LMarkerSetDataFactory = LSetDataFactory<LMarker>;
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LMarkerSetDataFactory
