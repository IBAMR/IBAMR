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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/IBTK_MPI.h>

#include <tbox/SAMRAI_MPI.h>
#include <tbox/Utilities.h>

#include <ostream>
#include <string>
#include <vector>

namespace IBTK
{
IBTK_MPI::comm IBTK_MPI::s_communicator = MPI_COMM_WORLD;

void
IBTK_MPI::setCommunicator(IBTK_MPI::comm communicator)
{
    s_communicator = communicator;
} // setCommunicator

IBTK_MPI::comm
IBTK_MPI::getCommunicator()
{
    return (s_communicator);
} // getCommunicator

IBTK_MPI::comm
IBTK_MPI::getSAMRAIWorld()
{
#if SAMRAI_VERSION_MAJOR == 2
    return SAMRAI::tbox::SAMRAI_MPI::commWorld;
#else
    return SAMRAI::tbox::SAMRAI_MPI::getSAMRAIWorld()
#endif
}

int
IBTK_MPI::getNodes(IBTK_MPI::comm communicator)
{
    int nodes = 1;
    MPI_Comm_size(communicator, &nodes);
    return nodes;
} // getNodes

int
IBTK_MPI::getRank(IBTK_MPI::comm communicator)
{
    int node = 0;
    MPI_Comm_rank(communicator, &node);
    return node;
} // getRank

void
IBTK_MPI::barrier(IBTK_MPI::comm communicator)
{
    (void)MPI_Barrier(communicator);
} // barrier

void
IBTK_MPI::allToOneSumReduction(int* x, const int n, const int root, IBTK_MPI::comm communicator)
{
    if (getNodes(communicator) > 1)
    {
        if (IBTK_MPI::getRank() == root)
            MPI_Reduce(MPI_IN_PLACE, x, n, MPI_INT, MPI_SUM, root, communicator);
        else
            MPI_Reduce(x, x, n, MPI_INT, MPI_SUM, root, communicator);
    }
} // allToOneSumReduction

void
IBTK_MPI::sendBytes(const void* buf,
                    const int number_bytes,
                    const int receiving_proc_number,
                    IBTK_MPI::comm communicator)
{
    MPI_Send((void*)buf, number_bytes, MPI_BYTE, receiving_proc_number, 0, communicator);
} // sendBytes

int
IBTK_MPI::recvBytes(void* buf, int number_bytes, IBTK_MPI::comm communicator)
{
    int rval = 0;
    MPI_Status status;
    MPI_Recv(buf, number_bytes, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, communicator, &status);

    rval = status.MPI_SOURCE;
    return rval;
} // recvBytes

//////////////////////////////////////  PRIVATE  ///////////////////////////////////////////////////

void
IBTK_MPI::allGatherSetup(int size_in,
                         int size_out,
                         std::vector<int>& rcounts,
                         std::vector<int>& disps,
                         IBTK_MPI::comm communicator)
{
    int np = getNodes(communicator);
    rcounts.resize(np);
    disps.resize(np);

    /* figure out where where each processor's input will be placed */
    allGather(size_in, rcounts.data(), communicator);

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
