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

// Check that a tbox::Array that cannot be allocated stops the program with an
// error.

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>

#include <tbox/Array.h>
#include <tbox/Logger.h>

#include <limits>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
// INT_MAX elements of 1 MiB need more memory than any system has.
struct Megabyte
{
    char bytes[1 << 20];
};
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<Logger::Appender> abort_append(new TestAppender());
    Logger::getInstance()->setAbortAppender(abort_append);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "array_04.log");

    tbox::Array<Megabyte> blocks(std::numeric_limits<int>::max());
    blocks[0].bytes[0] = 1;
    pout << "allocated " << blocks.size() << " blocks\n";
}
