// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Some linkers (like on macOS) complain if there is an empty static library. For
// backwards compatibility reasons we want to keep defining libIBAMR.a, even
// though nothing is in it - hence give it a single symbol to placate the linker.

namespace IBAMR
{
int
add_42(int a)
{
    return a + 42;
}
} // namespace IBAMR
