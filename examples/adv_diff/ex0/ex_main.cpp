// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <vector>

#include "example.cpp"

int
main(int argc, char** argv)
{
    std::vector<double> Q_err;
    run_example(argc, argv, Q_err);
    return 0;
}
