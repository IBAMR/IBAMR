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
#include <tbox/Utilities.h>

#include <mpi.h>

#include <string>

using IBTK::equals_ignore_case;

static_assert(equals_ignore_case("", ""));
static_assert(equals_ignore_case("BLAS-LAPACK", "blas-lapack"));
static_assert(equals_ignore_case("Multiplicative", "mULTIPLICATIVE"));
static_assert(!equals_ignore_case("additive", "additiv"));
static_assert(!equals_ignore_case("additive", "additivf"));
// Only letters change case: '@' and '`' neighbor the ranges of the letters.
static_assert(!equals_ignore_case("@", "`"));
static_assert(!equals_ignore_case("_", "-"));

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    SAMRAI::tbox::PIO::logOnlyNodeZero("output");
    const std::string mixed = "Eigen-Schur-Complement";
    TBOX_ASSERT(equals_ignore_case(mixed, "eigen-schur-complement"));
    TBOX_ASSERT(!equals_ignore_case(mixed, "eigen-schur-complements"));
    return 0;
}
