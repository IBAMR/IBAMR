// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/IBTK_MPI.h>

#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

#include <ostream>
#include <string>
#include <vector>

namespace IBTK
{
MPI_Comm IBTK_MPI::s_communicator = MPI_COMM_WORLD;

void
IBTK_MPI::setCommunicator(MPI_Comm communicator)
{
    s_communicator = communicator;
} // setCommunicator

MPI_Comm
IBTK_MPI::getCommunicator()
{
    return (s_communicator);
} // getCommunicator

MPI_Comm
IBTK_MPI::getSAMRAIWorld()
{
    return SAMRAI::tbox::SAMRAI_MPI::commWorld;
}

int
IBTK_MPI::getNodes()
{
    int nodes = 1;
    const int ierr = MPI_Comm_size(IBTK_MPI::getCommunicator(), &nodes);
    TBOX_ASSERT(ierr == MPI_SUCCESS);
    return nodes;
} // getNodes

int
IBTK_MPI::getRank()
{
    int node = 0;
    const int ierr = MPI_Comm_rank(IBTK_MPI::getCommunicator(), &node);
    TBOX_ASSERT(ierr == MPI_SUCCESS);
    return node;
} // getRank

void
IBTK_MPI::barrier()
{
    const int ierr = MPI_Barrier(IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == MPI_SUCCESS);
} // barrier

void
IBTK_MPI::allToOneSumReduction(int* x, const int n, const int root)
{
    if (getNodes() > 1)
    {
        if (IBTK_MPI::getRank() == root)
        {
            const int ierr = MPI_Reduce(MPI_IN_PLACE, x, n, MPI_INT, MPI_SUM, root, IBTK_MPI::getCommunicator());
            TBOX_ASSERT(ierr == MPI_SUCCESS);
        }
        else
        {
            const int ierr = MPI_Reduce(x, x, n, MPI_INT, MPI_SUM, root, IBTK_MPI::getCommunicator());
            TBOX_ASSERT(ierr == MPI_SUCCESS);
        }
    }
} // allToOneSumReduction

void
IBTK_MPI::sendBytes(const void* buf, const int number_bytes, const int receiving_proc_number)
{
    const int ierr =
        MPI_Send((void*)buf, number_bytes, MPI_BYTE, receiving_proc_number, 0, IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == MPI_SUCCESS);
} // sendBytes

int
IBTK_MPI::recvBytes(void* buf, int number_bytes)
{
    int rval = 0;
    MPI_Status status;
    const int ierr =
        MPI_Recv(buf, number_bytes, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, IBTK_MPI::getCommunicator(), &status);
    TBOX_ASSERT(ierr == MPI_SUCCESS);

    rval = status.MPI_SOURCE;
    return rval;
} // recvBytes

//////////////////////////////////////  PRIVATE  ///////////////////////////////////////////////////

void
IBTK_MPI::allGatherSetup(int size_in, int size_out, std::vector<int>& rcounts, std::vector<int>& disps)
{
    int np = getNodes();
    rcounts.resize(np);
    disps.resize(np);

    /* figure out where where each processor's input will be placed */
    allGather(size_in, rcounts.data());

    disps[0] = 0;
    for (int p = 1; p < np; ++p)
    {
        disps[p] = disps[p - 1] + rcounts[p - 1];
    }

    /* verify that the x_out array is the appropriate size! */
    int c = 0;
    for (int x = 0; x < np; ++x)
    {
        c += rcounts[x];
    }
    if (c != size_out)
    {
        TBOX_ERROR("IBTK_MPI::allGatherSetup error..."
                   << "\n   size_out =" << size_out << "appears to be incorrect; "
                   << "should be: " << c << std::endl);
    }
} // allGatherSetup
} // namespace IBTK
