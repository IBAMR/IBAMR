// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <tbox/PIO.h>

#include <mpi.h>

#include <string>
#include <vector>

// A utility function that prints @p out to plog by sending each string to
// rank 0.
inline void
print_strings_on_plog_0(const std::string& out)
{
    using namespace SAMRAI::tbox;
    const int n_nodes = SAMRAI_MPI::getNodes();
    std::vector<unsigned long> string_sizes(n_nodes);

    const unsigned long size = out.size();
    int ierr = MPI_Gather(
        &size, 1, MPI_UNSIGNED_LONG, string_sizes.data(), 1, MPI_UNSIGNED_LONG, 0, SAMRAI_MPI::getCommunicator());
    TBOX_ASSERT(ierr == 0);

    // MPI_Gatherv would be more efficient, but this just a test so its
    // not too important
    if (SAMRAI_MPI::getRank() == 0)
    {
        plog << out;
        for (int r = 1; r < n_nodes; ++r)
        {
            std::string input;
            input.resize(string_sizes[r]);
            ierr =
                MPI_Recv(&input[0], string_sizes[r], MPI_CHAR, r, 0, SAMRAI_MPI::getCommunicator(), MPI_STATUS_IGNORE);
            TBOX_ASSERT(ierr == 0);
            plog << input;
        }
    }
    else
        MPI_Send(out.data(), size, MPI_CHAR, 0, 0, SAMRAI_MPI::getCommunicator());
}
