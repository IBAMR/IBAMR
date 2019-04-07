// Filename: IBTK_MPI.cpp
//
// Copyright (c) 2002-2017, Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "tbox/SAMRAI_MPI.h"

#include "ibtk/IBTK_MPI.h"
#include "ibtk/app_namespaces.h"
#include "tbox/Utilities.h"

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
    return SAMRAI_MPI::commWorld;
#else
    return SAMRAI_MPI::getSAMRAIWorld()
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

template <typename T>
T
IBTK_MPI::minReduction(T x, int* rank_of_min, IBTK_MPI::comm communicator)
{
    minReduction(&x, 1, rank_of_min, communicator);
    return x;
} // minReduction

template <typename T>
void
IBTK_MPI::minReduction(T* x, const int n, int* rank_of_min, IBTK_MPI::comm communicator)
{
    if (n == 0) return;
    if (rank_of_min == nullptr)
    {
        MPI_Allreduce(MPI_IN_PLACE, x, n, mpi_type_id(x[0]), MPI_MIN, communicator);
    }
    else
    {
        minMaxReduction(x, n, rank_of_min, MPI_MINLOC, communicator);
    }
} // minReduction

template <typename T>
T
IBTK_MPI::maxReduction(T x, int* rank_of_max, IBTK_MPI::comm communicator)
{
    maxReduction(&x, 1, rank_of_max, communicator);
    return x;
} // maxReduction

template <typename T>
void
IBTK_MPI::maxReduction(T* x, const int n, int* rank_of_max, IBTK_MPI::comm communicator)
{
    if (n == 0) return;
    if (rank_of_max == nullptr)
    {
        MPI_Allreduce(MPI_IN_PLACE, x, n, mpi_type_id(x[0]), MPI_MAX, communicator);
    }
    else
    {
        minMaxReduction(x, n, rank_of_max, MPI_MAXLOC, communicator);
    }
} // maxReduction

template <typename T>
T
IBTK_MPI::sumReduction(T x, IBTK_MPI::comm communicator)
{
    sumReduction(&x, 1, communicator);
    return x;
} // sumReduction

template <typename T>
void
IBTK_MPI::sumReduction(T* x, const int n, IBTK_MPI::comm communicator)
{
    if (n == 0 || getNodes(communicator) < 2) return;
    MPI_Allreduce(MPI_IN_PLACE, x, n, mpi_type_id(x[0]), MPI_SUM, communicator);
} // sumReduction

void
IBTK_MPI::allToOneSumReduction(int* x, const int n, const int root, IBTK_MPI::comm communicator)
{
    if (getNodes(communicator) > 1)
    {
        MPI_Reduce(MPI_IN_PLACE, x, n, MPI_INT, MPI_SUM, root, communicator);
    }
} // allToOneSumReduction

template <typename T>
T
IBTK_MPI::bcast(const T x, const int root, IBTK_MPI::comm communicator)
{
    bcast(&x, 1, root, communicator);
}

template <typename T>
void
IBTK_MPI::bcast(T* x, int& length, const int root, IBTK_MPI::comm communicator)
{
    if (getNodes(communicator) > 1)
    {
        MPI_Bcast(x, length, mpi_type_id(x[0]), root, communicator);
    }
} // bcast

template <typename T>
void
IBTK_MPI::send(const T* buf,
               const int length,
               const int receiving_proc_number,
               const bool send_length,
               int tag,
               IBTK_MPI::comm communicator)
{
    tag = (tag >= 0) ? tag : 0;
    int size = length;
    if (send_length)
    {
        MPI_Send(&size, 1, MPI_INT, receiving_proc_number, tag, communicator);
    }
    MPI_Send(buf, length, mpi_type_id(buf[0]), receiving_proc_number, tag, communicator);
} // send

template <typename T>
void
IBTK_MPI::recv(T* buf,
               int& length,
               const int sending_proc_number,
               const bool get_length,
               int tag,
               IBTK_MPI::comm communicator)
{
    MPI_Status status;
    tag = (tag >= 0) ? tag : 0;
    if (get_length)
    {
        MPI_Recv(&length, 1, MPI_INT, sending_proc_number, tag, communicator, &status);
    }
    MPI_Recv(buf, length, mpi_type_id(buf[0]), sending_proc_number, tag, communicator, &status);
} // recv

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

template <typename T>
void
IBTK_MPI::allGather(const T* x_in, int size_in, T* x_out, int size_out, IBTK_MPI::comm communicator)
{
    std::vector<int> rcounts, disps;
    allGatherSetup(size_in, size_out, rcounts, disps, communicator);

    MPI_Allgatherv(
        x_in, size_in, mpi_type_id(x_in[0]), x_out, rcounts.data(), disps.data(), mpi_type_id(x_in[0]), communicator);
} // allGather

template <typename T>
void
IBTK_MPI::allGather(T x_in, T* x_out, IBTK_MPI::comm communicator)
{
    MPI_Allgather(&x_in, 1, mpi_type_id(x_in), x_out, 1, mpi_type_id(x_in), communicator);
} // allGather
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

template <typename T>
void
IBTK_MPI::minMaxReduction(T* x, const int n, int* rank, MPI_Op op, IBTK_MPI::comm communicator)
{
    std::vector<std::pair<T, int> > recv(n);
    std::vector<std::pair<T, int> > send(n);
    for (int i = 0; i < n; ++i)
    {
        send[i].first = x[i];
        send[i].second = getRank(communicator);
    }
    MPI_Allreduce(send.data(), recv.data(), n, mpi_type_id(recv[0]), op, communicator);
    for (int i = 0; i < n; ++i)
    {
        x[i] = recv[i].first;
        rank[i] = send[i].second;
    }
} // minMaxReduction
} // namespace IBTK
