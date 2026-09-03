// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_IBKernel_inl
#define included_IBTK_IBKernel_inl

namespace IBTK
{
inline bool
IBKernel::operator==(const IBKernel& other) const
{
    return d_name[0] == other.d_name[0] && d_name[1] == other.d_name[1];
}

inline bool
IBKernel::operator!=(const IBKernel& other) const
{
    return !(*this == other);
}
} // namespace IBTK
#endif
