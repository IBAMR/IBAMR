// ---------------------------------------------------------------------
//
// Copyright (c) 2025 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Ensure that we can staticly allocate and deallocate tbox::Arrays before and
// after some patches to SAMRAI.

#include <tbox/Array.h>

#include <fstream>

#include <ibamr/app_namespaces.h>

static tbox::Array<double> s_array_1;
static tbox::Array<tbox::Array<double> > s_array_2(10);
static tbox::Array<tbox::Array<double> > s_array_3(100);
static tbox::Array<tbox::Array<double> > s_array_4(1000);
static tbox::Array<tbox::Array<double> > s_array_5(10000);

int
main()
{
    std::ofstream output("output");

    output << "array 1 size: " << s_array_1.size() << '\n';
    output << "array 2 size: " << s_array_2.size() << '\n';
    output << "array 3 size: " << s_array_3.size() << '\n';
    output << "array 4 size: " << s_array_4.size() << '\n';
    output << "array 5 size: " << s_array_5.size() << '\n';
}
