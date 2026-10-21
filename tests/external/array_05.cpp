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

// Check that tbox::MessageStream accepts zero-length operations on full and
// empty streams.

#include <tbox/MessageStream.h>

#include <fstream>

#include <ibtk/app_namespaces.h>

int
main()
{
    std::ofstream output("output");

    const double values[2] = { 1.0, 2.0 };
    tbox::MessageStream full(sizeof(values), tbox::MessageStream::Write);
    full.pack(values, 2);
    full.pack(values, 0);
    output << "bytes in a full stream after packing no values: " << full.getCurrentSize() << '\n';

    tbox::MessageStream empty(0, tbox::MessageStream::Write);
    empty.getBufferStart();
    output << "bytes in an empty stream: " << empty.getCurrentSize() << '\n';
}
