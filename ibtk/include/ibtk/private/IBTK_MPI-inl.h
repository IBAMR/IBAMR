// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#ifndef included_IBTK_MPI_inc
#define included_IBTK_MPI_inc

#include "ibtk/IBTK_MPI.h"

namespace IBTK
{
template <typename T>
inline T
IBTK_MPI::minReduction(T x, int* rank_of_min, IBTK_MPI::comm communicator)
{
    minReduction(&x, 1, rank_of_min, communicator);
    return x;
} // minReduction

template <typename T>
inline void
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
inline T
IBTK_MPI::maxReduction(T x, int* rank_of_max, IBTK_MPI::comm communicator)
{
    maxReduction(&x, 1, rank_of_max, communicator);
    return x;
} // maxReduction

template <typename T>
inline void
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
inline T
IBTK_MPI::sumReduction(T x, IBTK_MPI::comm communicator)
{
    sumReduction(&x, 1, communicator);
    return x;
} // sumReduction

template <typename T>
inline void
IBTK_MPI::sumReduction(T* x, const int n, IBTK_MPI::comm communicator)
{
    if (n == 0 || getNodes(communicator) < 2) return;
    MPI_Allreduce(MPI_IN_PLACE, x, n, mpi_type_id(x[0]), MPI_SUM, communicator);
} // sumReduction

template <typename T>
inline T
IBTK_MPI::bcast(const T x, const int root, IBTK_MPI::comm communicator)
{
    int size = 1;
    T temp_copy = x;
    bcast(&temp_copy, size, root, communicator);
    return x;
}

template <typename T>
inline void
IBTK_MPI::bcast(T* x, int& length, const int root, IBTK_MPI::comm communicator)
{
    if (getNodes(communicator) > 1)
    {
        MPI_Bcast(x, length, mpi_type_id(x[0]), root, communicator);
    }
} // bcast

template <typename T>
inline void
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
inline void
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

template <typename T>
inline void
IBTK_MPI::allGather(const T* x_in, int size_in, T* x_out, int size_out, IBTK_MPI::comm communicator)
{
    std::vector<int> rcounts, disps;
    allGatherSetup(size_in, size_out, rcounts, disps, communicator);

    MPI_Allgatherv(
        x_in, size_in, mpi_type_id(x_in[0]), x_out, rcounts.data(), disps.data(), mpi_type_id(x_in[0]), communicator);
} // allGather

template <typename T>
inline void
IBTK_MPI::allGather(T x_in, T* x_out, IBTK_MPI::comm communicator)
{
    MPI_Allgather(&x_in, 1, mpi_type_id(x_in), x_out, 1, mpi_type_id(x_in), communicator);
} // allGather

//////////////////////////////////////  PRIVATE  ///////////////////////////////////////////////////
template <typename T>
inline void
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

#endif
