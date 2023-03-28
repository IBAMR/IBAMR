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

#ifndef included_IBTK_LNodeSetDataIterator
#define included_IBTK_LNodeSetDataIterator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LNode.h"
#include "ibtk/LSetDataIterator.h"

/////////////////////////////// TYPEDEFS /////////////////////////////////////

namespace IBTK
{
using LNodeSetDataIterator = LSetDataIterator<LNode>;
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_LNodeSetDataIterator
