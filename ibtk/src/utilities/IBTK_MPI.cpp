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

#include <stdlib.h>
#include <string.h>

#include <string>

#include "tbox/SAMRAI_MPI.h"

#include "ibtk/IBTK_MPI.h"
#include "ibtk/app_namespaces.h"
#include "tbox/Utilities.h"

namespace IBTK
{
IBTK_MPI::comm IBTK_MPI::d_communicator = MPI_COMM_WORLD;

void
IBTK_MPI::setCommunicator(IBTK_MPI::comm communicator)
{
    d_communicator = communicator;
} // setCommunicator

IBTK_MPI::comm
IBTK_MPI::getCommunicator()
{
    return (d_communicator);
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
IBTK_MPI::getNodes(IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    int nodes = 1;
    MPI_Comm_size(communicator, &nodes);
    return nodes;
} // getNodes

int
IBTK_MPI::getRank(IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    int node = 0;
    MPI_Comm_rank(communicator, &node);
    return node;
} // getRank

void
IBTK_MPI::barrier(IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    (void)MPI_Barrier(communicator);
} // barrier

void
IBTK_MPI::minReduction(double* x, const int n, int* rank_of_min, IBTK_MPI::comm communicator)
{
    minReduction<double>(x, n, rank_of_min, communicator);
}

void
IBTK_MPI::minReduction(int* x, const int n, int* rank_of_min, IBTK_MPI::comm communicator)
{
    minReduction<int>(x, n, rank_of_min, communicator);
}

void
IBTK_MPI::minReduction(float* x, const int n, int* rank_of_min, IBTK_MPI::comm communicator)
{
    minReduction<float>(x, n, rank_of_min, communicator);
}

double
IBTK_MPI::minReduction(double x, int* rank_of_min, IBTK_MPI::comm communicator)
{
    minReduction<double>(&x, 1, rank_of_min, communicator);
    return x;
}

int
IBTK_MPI::minReduction(int x, int* rank_of_min, IBTK_MPI::comm communicator)
{
    minReduction<int>(&x, 1, rank_of_min, communicator);
    return x;
}

float
IBTK_MPI::minReduction(float x, int* rank_of_min, IBTK_MPI::comm communicator)
{
    minReduction<float>(&x, 1, rank_of_min, communicator);
    return x;
}

void
IBTK_MPI::maxReduction(double* x, const int n, int* rank_of_max, IBTK_MPI::comm communicator)
{
    maxReduction<double>(x, n, rank_of_max, communicator);
}

void
IBTK_MPI::maxReduction(int* x, const int n, int* rank_of_max, IBTK_MPI::comm communicator)
{
    maxReduction<int>(x, n, rank_of_max, communicator);
}

void
IBTK_MPI::maxReduction(float* x, const int n, int* rank_of_max, IBTK_MPI::comm communicator)
{
    maxReduction<float>(x, n, rank_of_max, communicator);
}

double
IBTK_MPI::maxReduction(double x, int* rank_of_max, IBTK_MPI::comm communicator)
{
    maxReduction<double>(&x, 1, rank_of_max, communicator);
    return x;
}

int
IBTK_MPI::maxReduction(int x, int* rank_of_max, IBTK_MPI::comm communicator)
{
    maxReduction<int>(&x, 1, rank_of_max, communicator);
    return x;
}

float
IBTK_MPI::maxReduction(float x, int* rank_of_max, IBTK_MPI::comm communicator)
{
    maxReduction<float>(&x, 1, rank_of_max, communicator);
    return x;
}

void
IBTK_MPI::sumReduction(double* x, const int n, IBTK_MPI::comm communicator)
{
    sumReduction<double>(x, n, communicator);
}

void
IBTK_MPI::sumReduction(int* x, const int n, IBTK_MPI::comm communicator)
{
    sumReduction<int>(x, n, communicator);
}

void
IBTK_MPI::sumReduction(float* x, const int n, IBTK_MPI::comm communicator)
{
    sumReduction<float>(x, n, communicator);
}

double
IBTK_MPI::sumReduction(double x, IBTK_MPI::comm communicator)
{
    sumReduction<double>(&x, 1, communicator);
    return x;
}

int
IBTK_MPI::sumReduction(int x, IBTK_MPI::comm communicator)
{
    sumReduction<int>(&x, 1, communicator);
    return x;
}

float
IBTK_MPI::sumReduction(float x, IBTK_MPI::comm communicator)
{
    sumReduction<float>(&x, 1, communicator);
    return x;
}

void
IBTK_MPI::allToOneSumReduction(int* x, const int n, const int root, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        int* send = new int[n];
        memcpy(send, x, n * sizeof(int));
        MPI_Reduce(send, x, n, MPI_INT, MPI_SUM, root, communicator);
        delete[] send;
    }
} // allToOneSumReduction

int
IBTK_MPI::bcast(const int x, const int root, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    int recv = x;
    if (getNodes(communicator) > 1)
    {
        (void)MPI_Bcast(&recv, 1, MPI_INT, root, communicator);
    }
    return (recv);
} // bcast

void
IBTK_MPI::bcast(int* x, int& length, const int root, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        (void)MPI_Bcast((void*)x, length, MPI_INT, root, communicator);
    }
} // bcast

void
IBTK_MPI::bcast(char* x, int& length, const int root, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        (void)MPI_Bcast((void*)x, length, MPI_BYTE, root, communicator);
    }
} // bcast

void
IBTK_MPI::send(const int* buf,
               const int length,
               const int receiving_proc_number,
               const bool send_length,
               int tag,
               IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    tag = (tag >= 0) ? tag : 0;
    int size = length;
    if (send_length)
    {
        MPI_Send(&size, 1, MPI_INT, receiving_proc_number, tag, communicator);
    }
    MPI_Send((void*)buf, length, MPI_INT, receiving_proc_number, tag, communicator);
} // send

void
IBTK_MPI::recv(int* buf,
               int& length,
               const int sending_proc_number,
               const bool get_length,
               int tag,
               IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    MPI_Status status;
    tag = (tag >= 0) ? tag : 0;
    if (get_length)
    {
        MPI_Recv(&length, 1, MPI_INT, sending_proc_number, tag, communicator, &status);
    }
    MPI_Recv((void*)buf, length, MPI_INT, sending_proc_number, tag, communicator, &status);
} // recv

void
IBTK_MPI::sendBytes(const void* buf,
                    const int number_bytes,
                    const int receiving_proc_number,
                    IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    MPI_Send((void*)buf, number_bytes, MPI_BYTE, receiving_proc_number, 0, communicator);
} // sendBytes

int
IBTK_MPI::recvBytes(void* buf, int number_bytes, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    int rval = 0;
    MPI_Status status;
    MPI_Recv(buf, number_bytes, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, communicator, &status);

    rval = status.MPI_SOURCE;
    return rval;
} // recvBytes

void
IBTK_MPI::allGather(const int* x_in,
                    int size_in,
                    int* x_out,
                    int size_out,
                    IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    int* rcounts = (int*)NULL;
    int* disps = (int*)NULL;
    allGatherSetup(size_in, size_out, rcounts, disps, communicator);

    MPI_Allgatherv((void*)x_in, size_in, MPI_INT, x_out, rcounts, disps, MPI_INT, communicator);

    if (rcounts)
    {
        delete[] rcounts;
    }
    if (disps)
    {
        delete[] disps;
    }
} // allGather

void
IBTK_MPI::allGather(const double* x_in,
                    int size_in,
                    double* x_out,
                    int size_out,
                    IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    int* rcounts = (int*)NULL;
    int* disps = (int*)NULL;
    allGatherSetup(size_in, size_out, rcounts, disps, communicator);

    MPI_Allgatherv((void*)x_in, size_in, MPI_DOUBLE, x_out, rcounts, disps, MPI_DOUBLE, communicator);

    if (rcounts)
    {
        delete[] rcounts;
    }
    if (disps)
    {
        delete[] disps;
    }
} // allGather

void
IBTK_MPI::allGather(double x_in, double* x_out, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    MPI_Allgather(&x_in, 1, MPI_DOUBLE, x_out, 1, MPI_DOUBLE, communicator);
} // allGather

void
IBTK_MPI::allGather(int x_in, int* x_out, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    MPI_Allgather(&x_in, 1, MPI_INT, x_out, 1, MPI_INT, communicator);
} // allGather
//////////////////////////////////////  PRIVATE  ///////////////////////////////////////////////////

template <typename T>
void
IBTK_MPI::minReduction(T* x, const int n, int* rank_of_min, IBTK_MPI::comm communicator)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (n == 0) return;
    if (rank_of_min == NULL)
    {
        std::vector<T> send(n);
        std::copy(x, x + n, send.begin());
        MPI_Allreduce(send.data(), x, n, mpi_type_id(send[0]), MPI_MIN, communicator);
    }
    else
    {
        // TODO: both this approach and SAMRAI's assume that we can define
        // byte-for-byte compatible types with, e.g.,
        // MPI_FLOAT_INT. Check this with a unit test!
        std::vector<std::pair<T, int> > recv(n);
        std::vector<std::pair<T, int> > send(n);
        for (int i = 0; i < n; ++i)
        {
            send[i].first = x[i];
            send[i].second = getRank(communicator);
        }
        MPI_Allreduce(send.data(), recv.data(), n, mpi_type_id(recv[0]), MPI_MINLOC, communicator);
        for (int i = 0; i < n; ++i)
        {
            x[i] = recv[i].first;
            rank_of_min[i] = send[i].second;
        }
    }
} // minReduction

template <typename T>
void
IBTK_MPI::maxReduction(T* x, const int n, int* rank_of_max, IBTK_MPI::comm communicator)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (n == 0) return;
    if (rank_of_max == NULL)
    {
        std::vector<T> send(n);
        std::copy(x, x + n, send.begin());
        MPI_Allreduce(send.data(), x, n, mpi_type_id(send[0]), MPI_MAX, communicator);
    }
    else
    {
        // TODO: both this approach and SAMRAI's assume that we can define
        // byte-for-byte compatible types with, e.g.,
        // MPI_FLOAT_INT. Check this with a unit test!
        std::vector<std::pair<T, int> > recv(n);
        std::vector<std::pair<T, int> > send(n);
        for (int i = 0; i < n; ++i)
        {
            send[i].first = x[i];
            send[i].second = getRank(communicator);
        }
        MPI_Allreduce(send.data(), recv.data(), n, mpi_type_id(recv[0]), MPI_MAXLOC, communicator);
        for (int i = 0; i < n; ++i)
        {
            x[i] = recv[i].first;
            rank_of_max[i] = send[i].second;
        }
    }
} // maxReduction

template <typename T>
void
IBTK_MPI::sumReduction(T* x, const int n, IBTK_MPI::comm communicator)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (n == 0 || getNodes(communicator) < 2) return;
    std::vector<T> send(n);
    std::copy(x, x + n, send.begin());
    MPI_Allreduce(send.data(), x, n, mpi_type_id(send[0]), MPI_SUM, communicator);
} // sumReduction

void
IBTK_MPI::allGatherSetup(int size_in,
                         int size_out,
                         int*& rcounts,
                         int*& disps,
                         IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    int np = getNodes(communicator);
    rcounts = new int[np];
    disps = new int[np];

    /* figure out where where each processor's input will be placed */
    allGather(size_in, rcounts, communicator);

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
