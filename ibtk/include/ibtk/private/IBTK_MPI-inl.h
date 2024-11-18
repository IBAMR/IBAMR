// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_MPI_inc
#define included_IBTK_MPI_inc

#include <ibtk/config.h>

#include <ibtk/IBTK_MPI.h>

#include <tbox/Utilities.h>

namespace IBTK
{
template <typename T>
inline T
IBTK_MPI::minReduction(T x, int* rank_of_min)
{
    minReduction(&x, 1, rank_of_min);
    return x;
} // minReduction

template <typename T>
inline void
IBTK_MPI::minReduction(T* x, const int n, int* rank_of_min)
{
    if (rank_of_min == nullptr)
    {
        const int ierr = MPI_Allreduce(MPI_IN_PLACE, x, n, mpi_type_id(T{}), MPI_MIN, IBTK_MPI::getCommunicator());
        TBOX_ASSERT(ierr == MPI_SUCCESS);
    }
    else
    {
        minMaxReduction(x, n, rank_of_min, MPI_MINLOC);
    }
} // minReduction

template <typename T>
inline T
IBTK_MPI::maxReduction(T x, int* rank_of_max)
{
    maxReduction(&x, 1, rank_of_max);
    return x;
} // maxReduction

template <typename T>
inline void
IBTK_MPI::maxReduction(T* x, const int n, int* rank_of_max)
{
    if (rank_of_max == nullptr)
    {
        const int ierr = MPI_Allreduce(MPI_IN_PLACE, x, n, mpi_type_id(T{}), MPI_MAX, IBTK_MPI::getCommunicator());
        TBOX_ASSERT(ierr == MPI_SUCCESS);
    }
    else
    {
        minMaxReduction(x, n, rank_of_max, MPI_MAXLOC);
    }
} // maxReduction

template <typename T>
inline T
IBTK_MPI::sumReduction(T x)
{
    // This check is useful since it occurs earlier than the mpi_type_id() base
    // case failure (and it has a clearer error message)
    static_assert(!std::is_pointer<T>::value,
                  "This function cannot be instantiated for pointer types "
                  "since it does not make sense to sum pointers.");
    sumReduction(&x, 1);
    return x;
} // sumReduction

template <typename T>
inline void
IBTK_MPI::sumReduction(T* x, const int n)
{
    const int ierr = MPI_Allreduce(MPI_IN_PLACE, x, n, mpi_type_id(T{}), MPI_SUM, IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == MPI_SUCCESS);
} // sumReduction

template <typename T>
inline T
IBTK_MPI::bcast(const T x, const int root)
{
    int size = 1;
    T temp_copy = x;
    bcast(&temp_copy, size, root);
    return temp_copy;
}

template <typename T>
inline void
IBTK_MPI::bcast(T* x, int& length, const int root)
{
    const int ierr = MPI_Bcast(x, length, mpi_type_id(T{}), root, IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == MPI_SUCCESS);
} // bcast

template <typename T>
inline void
IBTK_MPI::send(const T* buf, const int length, const int receiving_proc_number, const bool send_length, int tag)
{
    TBOX_ASSERT(tag >= 0);
    if (send_length)
    {
        const int ierr = MPI_Send(&length, 1, MPI_INT, receiving_proc_number, tag, IBTK_MPI::getCommunicator());
        TBOX_ASSERT(ierr == MPI_SUCCESS);
    }
    const int ierr =
        MPI_Send(buf, length, mpi_type_id(buf[0]), receiving_proc_number, tag, IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == MPI_SUCCESS);
} // send

template <typename T>
inline void
IBTK_MPI::recv(T* buf, int& length, const int sending_proc_number, const bool get_length, int tag)
{
    TBOX_ASSERT(tag >= 0);
    MPI_Status status;
    if (get_length)
    {
        const int ierr = MPI_Recv(&length, 1, MPI_INT, sending_proc_number, tag, IBTK_MPI::getCommunicator(), &status);
        TBOX_ASSERT(ierr == MPI_SUCCESS);
    }
    const int ierr =
        MPI_Recv(buf, length, mpi_type_id(buf[0]), sending_proc_number, tag, IBTK_MPI::getCommunicator(), &status);
    TBOX_ASSERT(ierr == MPI_SUCCESS);
} // recv

template <typename T>
inline void
IBTK_MPI::allGather(const T* x_in, int size_in, T* x_out, int size_out)
{
    std::vector<int> rcounts, disps;
    allGatherSetup(size_in, size_out, rcounts, disps);

    const int ierr = MPI_Allgatherv(x_in,
                                    size_in,
                                    mpi_type_id(T{}),
                                    x_out,
                                    rcounts.data(),
                                    disps.data(),
                                    mpi_type_id(T{}),
                                    IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == MPI_SUCCESS);
} // allGather

template <typename T>
inline void
IBTK_MPI::allGather(T x_in, T* x_out)
{
    const int ierr = MPI_Allgather(&x_in, 1, mpi_type_id(T{}), x_out, 1, mpi_type_id(T{}), IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == MPI_SUCCESS);
} // allGather

//////////////////////////////////////  PRIVATE  ///////////////////////////////////////////////////
template <typename T>
inline void
IBTK_MPI::minMaxReduction(T* x, const int n, int* rank, MPI_Op op)
{
    std::vector<std::pair<T, int> > recv(n);
    std::vector<std::pair<T, int> > send(n);
    for (int i = 0; i < n; ++i)
    {
        send[i].first = x[i];
        send[i].second = getRank();
    }
    const int ierr =
        MPI_Allreduce(send.data(), recv.data(), n, mpi_type_id(std::pair<T, int>{}), op, IBTK_MPI::getCommunicator());
    TBOX_ASSERT(ierr == MPI_SUCCESS);
    for (int i = 0; i < n; ++i)
    {
        x[i] = recv[i].first;
        rank[i] = recv[i].second;
    }
} // minMaxReduction
} // namespace IBTK

#endif
