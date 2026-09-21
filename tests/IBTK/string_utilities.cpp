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

#include <ibtk/IBTKInit.h>
#include <ibtk/string_utilities.h>

#include <tbox/PIO.h>

#include <mpi.h>

#include <string>

using IBTK::equals_ignore_case;

static_assert(equals_ignore_case("", ""));
static_assert(equals_ignore_case("Mixed-Case", "mIXED-cASE"));
static_assert(!equals_ignore_case("mixed", "mixe"));
static_assert(!equals_ignore_case("mixed", "mixef"));
// '@' and '`' differ by the same offset as 'A' and 'a' but are not letters.
static_assert(!equals_ignore_case("@", "`"));
// The second bytes of these two-byte UTF-8 characters differ by the same offset.
static_assert(!equals_ignore_case("\xC3\x89", "\xC3\xA9"));

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    SAMRAI::tbox::PIO::logOnlyNodeZero("output");
    const std::string mixed = "Mixed-Case";
    SAMRAI::tbox::plog << "equal = " << equals_ignore_case(mixed, "mixed-case") << '\n';
    return 0;
}
