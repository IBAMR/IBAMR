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

// Check that a negative tbox::Array size stops with a message.

#include <ibtk/IBTKInit.h>

#include <tbox/Array.h>
#include <tbox/Logger.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

int
main(int argc, char** argv)
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    tbox::PIO::logOnlyNodeZero("output");
    tbox::Pointer<tbox::Logger::Appender> abort_append(new TestAppender());
    tbox::Logger::getInstance()->setAbortAppender(abort_append);
    tbox::Array<double> negative(-1);
    tbox::plog << "accepted\n";
}
