// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/SAMRAIManager.h>

// Set up application namespace declarations
#include <boost/core/demangle.hpp>

#include <fstream>
#include <iostream>
#include <utility>

#include <ibamr/app_namespaces.h>

// Verify that the correspondence between MPI and C++ types in IBTK_MPI works.

template <typename T>
void
size_check()
{
    const std::string t_str = boost::core::demangle(typeid(T).name());
    int size = 0;
    TBOX_ASSERT(MPI_Type_size(IBTK::mpi_type_id(T{}), &size) == MPI_SUCCESS);
    plog << '\n';
    plog << "size of MPI type = " << size << '\n';
    plog << "sizeof(" << t_str << ") = " << sizeof(T) << '\n';
}

template <typename U, typename T>
void
size_check_pair()
{
    // Check the size against the official C type too
    struct
    {
        U a;
        T b;
    } x;
    using pair_type = std::pair<U, T>;
    static_assert(sizeof(pair_type) == sizeof(x), "unequal sizes");
    int size = 0;
    TBOX_ASSERT(MPI_Type_size(IBTK::mpi_type_id(pair_type{}), &size) == MPI_SUCCESS);

    const std::string u_str = boost::core::demangle(typeid(U).name());
    const std::string t_str = boost::core::demangle(typeid(T).name());
    plog << '\n';
    plog << "size of MPI type = " << size << '\n';
    plog << "sizeof(std::pair<" << u_str << ", " << t_str << ">) = " << sizeof(pair_type) << '\n';
    plog << "sizeof(struct {" << u_str << " a; " << t_str << " b;}) = " << sizeof(x) << '\n';
}

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");

    size_check<double>();
    size_check<int>();
    size_check<float>();
    size_check<char>();

    size_check_pair<float, int>();
    // size_check_pair<long, int>();
    size_check_pair<double, int>();
    size_check_pair<int, int>();
    // size_check_pair<long double, int>();
    plog << "OK\n";
} // main
