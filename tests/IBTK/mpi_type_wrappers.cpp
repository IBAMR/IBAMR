// Copyright (c) 2019, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/SAMRAIManager.h>
#include <tbox/SAMRAI_MPI.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

#include <boost/core/demangle.hpp>

#include <fstream>
#include <iostream>
#include <utility>

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
    MPI_Init(&argc, &argv);
    SAMRAI_MPI::setCommunicator(MPI_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

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
    MPI_Finalize();
} // main
