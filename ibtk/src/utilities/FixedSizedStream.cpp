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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibtk/FixedSizedStream.h"

#include <cstring>
#include <vector>

#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

FixedSizedStream::FixedSizedStream(const int bytes) : d_buffer_size(bytes), d_buffer(d_buffer_size)
{
    // intentionally blank
    return;
} // FixedSizedStream

FixedSizedStream::FixedSizedStream(const void* const buffer, const int bytes)
    : d_buffer_size(bytes), d_buffer(d_buffer_size)
{
    std::memcpy(static_cast<void*>(&d_buffer[0]), buffer, bytes);
    return;
} // FixedSizedStream

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
