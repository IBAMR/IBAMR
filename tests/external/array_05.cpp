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

// Check that tbox::Array and tbox::MessageStream accept sizes of zero and
// below.

#include <tbox/Array.h>
#include <tbox/MessageStream.h>

#include <fstream>

#include <ibtk/app_namespaces.h>

int
main()
{
    std::ofstream output("output");

    tbox::Array<double> negative(-1);
    output << "size of an array of -1 elements: " << negative.size() << '\n';

    tbox::Array<double> resized(3);
    resized.resizeArray(-1);
    output << "size after resizing to -1 elements: " << resized.size() << '\n';

    const double values[2] = { 1.0, 2.0 };
    tbox::MessageStream full(sizeof(values), tbox::MessageStream::Write);
    full.pack(values, 2);
    full.pack(values, 0);
    output << "bytes in a full stream after packing no values: " << full.getCurrentSize() << '\n';

    tbox::MessageStream empty(0, tbox::MessageStream::Write);
    empty.getBufferStart();
    output << "bytes in an empty stream: " << empty.getCurrentSize() << '\n';
}
