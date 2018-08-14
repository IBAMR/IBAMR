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

#include "ibtk/IBTK_MPI.h"
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

void
IBTK_MPI::barrier(IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    (void)MPI_Barrier(communicator);
} // barrier

double
IBTK_MPI::sumReduction(const double x, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    double recv = x;
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        double send = x;
        MPI_Allreduce(&send, &recv, 1, MPI_DOUBLE, MPI_SUM, communicator);
    }
    return (recv);
} // sumReduction

void
IBTK_MPI::sumReduction(double* x, const int n, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        double* send = new double[n];
        memcpy(send, x, n * sizeof(double));
        MPI_Allreduce(send, x, n, MPI_DOUBLE, MPI_SUM, communicator);
        delete[] send;
    }
} // sumReduction

float
IBTK_MPI::sumReduction(const float x, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    float recv = x;
    if (getNodes(communicator) > 1)
    {
        float send = x;
        MPI_Allreduce(&send, &recv, 1, MPI_FLOAT, MPI_SUM, communicator);
    }
    return (recv);
} // sumReduction

void
IBTK_MPI::sumReduction(float* x, const int n, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        float* send = new float[n];
        memcpy(send, x, n * sizeof(float));
        MPI_Allreduce(send, x, n, MPI_FLOAT, MPI_SUM, communicator);
        delete[] send;
    }
} // sumReduction

int
IBTK_MPI::sumReduction(const int x, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    int recv = x;
    if (getNodes(communicator) > 1)
    {
        int send = x;
        MPI_Allreduce(&send, &recv, 1, MPI_INT, MPI_SUM, communicator);
    }
    return (recv);
} // sumReduction

void
IBTK_MPI::sumReduction(int* x, const int n, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        int* send = new int[n];
        memcpy(send, x, n * sizeof(int));
        MPI_Allreduce(send, x, n, MPI_INT, MPI_SUM, communicator);
        delete[] send;
    }
} // sumReduction

double
IBTK_MPI::minReduction(const double x, int* rank_of_min, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    double rval = x;
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    /*
     * If a rank_of_min argument is provided, set it to the current
     * rank of the process.
     */
    if (rank_of_min != NULL)
    {
        *rank_of_min = getRank(communicator);
    }

    if (getNodes(communicator) > 1)
    {
        if (rank_of_min == NULL)
        {
            double send = x;
            MPI_Allreduce(&send, &rval, 1, MPI_DOUBLE, MPI_MIN, communicator);
        }
        else
        {
            DoubleIntStruct recv;
            DoubleIntStruct send;
            send.d = x;
            send.i = getRank();
            MPI_Allreduce(&send, &recv, 1, MPI_DOUBLE_INT, MPI_MINLOC, communicator);
            rval = recv.d;
            *rank_of_min = recv.i;
        }
    }
    return (rval);
} // minReduction

void
IBTK_MPI::minReduction(double* x, const int n, int* rank_of_min, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        if (rank_of_min == NULL)
        {
            double* send = new double[n];
            memcpy(send, x, n * sizeof(double));
            MPI_Allreduce(send, x, n, MPI_DOUBLE, MPI_MIN, communicator);
            delete[] send;
        }
        else
        {
            DoubleIntStruct* recv = new DoubleIntStruct[n];
            DoubleIntStruct* send = new DoubleIntStruct[n];
            for (int i = 0; i < n; ++i)
            {
                send[i].d = x[i];
                send[i].i = getRank(communicator);
            }
            MPI_Allreduce(send, recv, n, MPI_DOUBLE_INT, MPI_MINLOC, communicator);
            for (int i = 0; i < n; ++i)
            {
                x[i] = recv[i].d;
                rank_of_min[i] = send[i].i;
            }
            delete[] recv;
            delete[] send;
        }
    }
} // minReduction

float
IBTK_MPI::minReduction(const float x, int* rank_of_min, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    float rval = x;

    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    /*
     * If a rank_of_min argument is provided, set it to the current
     * rank of the process.
     */
    if (rank_of_min != NULL)
    {
        *rank_of_min = getRank(communicator);
    }

    if (getNodes(communicator) > 1)
    {
        if (rank_of_min == NULL)
        {
            float send = x;
            MPI_Allreduce(&send, &rval, 1, MPI_FLOAT, MPI_MIN, communicator);
        }
        else
        {
            FloatIntStruct recv;
            FloatIntStruct send;
            send.f = x;
            send.i = getRank(communicator);
            MPI_Allreduce(&send, &recv, 1, MPI_FLOAT_INT, MPI_MINLOC, communicator);
            rval = recv.f;
            *rank_of_min = recv.i;
        }
    }
    return (rval);
} // minReduction

void
IBTK_MPI::minReduction(float* x, const int n, int* rank_of_min, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        if (rank_of_min == NULL)
        {
            float* send = new float[n];
            memcpy(send, x, n * sizeof(float));
            MPI_Allreduce(send, x, n, MPI_FLOAT, MPI_MIN, communicator);
            delete[] send;
        }
        else
        {
            FloatIntStruct* recv = new FloatIntStruct[n];
            FloatIntStruct* send = new FloatIntStruct[n];
            for (int i = 0; i < n; ++i)
            {
                send[i].f = x[i];
                send[i].i = getRank(communicator);
            }
            MPI_Allreduce(send, recv, n, MPI_FLOAT_INT, MPI_MINLOC, communicator);
            for (int i = 0; i < n; ++i)
            {
                x[i] = recv[i].f;
                rank_of_min[i] = send[i].i;
            }
            delete[] recv;
            delete[] send;
        }
    }
} // minReduction

int
IBTK_MPI::minReduction(const int x, int* rank_of_min, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    int rval = x;

    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    /*
     * If a rank_of_min argument is provided, set it to the current
     * rank of the process.
     */
    if (rank_of_min != NULL)
    {
        *rank_of_min = getRank(communicator);
    }

    if (getNodes(communicator) > 1)
    {
        if (rank_of_min == NULL)
        {
            int send = x;
            MPI_Allreduce(&send, &rval, 1, MPI_INT, MPI_MIN, communicator);
        }
        else
        {
            IntIntStruct recv;
            IntIntStruct send;
            send.j = x;
            send.i = getRank(communicator);
            MPI_Allreduce(&send, &recv, 1, MPI_2INT, MPI_MINLOC, communicator);
            rval = recv.j;
            *rank_of_min = recv.i;
        }
    }
    return (rval);
} // minReduction

void
IBTK_MPI::minReduction(int* x, const int n, int* rank_of_min, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        if (rank_of_min == NULL)
        {
            int* send = new int[n];
            memcpy(send, x, n * sizeof(int));
            MPI_Allreduce(send, x, n, MPI_INT, MPI_MIN, communicator);
            delete[] send;
        }
        else
        {
            IntIntStruct* recv = new IntIntStruct[n];
            IntIntStruct* send = new IntIntStruct[n];
            for (int i = 0; i < n; ++i)
            {
                send[i].j = x[i];
                send[i].i = getRank(communicator);
            }
            MPI_Allreduce(send, recv, n, MPI_2INT, MPI_MINLOC, communicator);
            for (int i = 0; i < n; ++i)
            {
                x[i] = recv[i].j;
                rank_of_min[i] = send[i].i;
            }
            delete[] recv;
            delete[] send;
        }
    }
} // minReduction

double
IBTK_MPI::maxReduction(const double x, int* rank_of_max, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    double rval = x;

    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    /*
     * If a rank_of_max argument is provided, set it to the current
     * rank of the process.
     */
    if (rank_of_max != NULL)
    {
        *rank_of_max = getRank(communicator);
    }

    if (getNodes() > 1)
    {
        if (rank_of_max == NULL)
        {
            double send = x;
            MPI_Allreduce(&send, &rval, 1, MPI_DOUBLE, MPI_MAX, communicator);
        }
        else
        {
            DoubleIntStruct recv;
            DoubleIntStruct send;
            send.d = x;
            send.i = getRank(communicator);
            MPI_Allreduce(&send, &recv, 1, MPI_DOUBLE_INT, MPI_MAXLOC, communicator);
            rval = recv.d;
            *rank_of_max = recv.i;
        }
    }
    return (rval);
} // maxReduction

void
IBTK_MPI::maxReduction(double* x, const int n, int* rank_of_max, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        if (rank_of_max == NULL)
        {
            double* send = new double[n];
            memcpy(send, x, n * sizeof(double));
            MPI_Allreduce(send, x, n, MPI_DOUBLE, MPI_MAX, communicator);
            delete[] send;
        }
        else
        {
            DoubleIntStruct* recv = new DoubleIntStruct[n];
            DoubleIntStruct* send = new DoubleIntStruct[n];
            for (int i = 0; i < n; ++i)
            {
                send[i].d = x[i];
                send[i].i = getRank(communicator);
            }
            MPI_Allreduce(send, recv, n, MPI_DOUBLE_INT, MPI_MAXLOC, communicator);
            for (int i = 0; i < n; ++i)
            {
                x[i] = recv[i].d;
                rank_of_max[i] = send[i].i;
            }
            delete[] recv;
            delete[] send;
        }
    }
} // maxReduction

float
IBTK_MPI::maxReduction(const float x, int* rank_of_max, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    float rval = x;

    /*
     * If a rank_of_max argument is provided, set it to the current
     * rank of the process.
     */
    if (rank_of_max != NULL)
    {
        *rank_of_max = getRank(communicator);
    }

    if (getNodes(communicator) > 1)
    {
        if (rank_of_max == NULL)
        {
            float send = x;
            MPI_Allreduce(&send, &rval, 1, MPI_FLOAT, MPI_MAX, communicator);
        }
        else
        {
            FloatIntStruct recv;
            FloatIntStruct send;
            send.f = x;
            send.i = getRank(communicator);
            MPI_Allreduce(&send, &recv, 1, MPI_FLOAT_INT, MPI_MAXLOC, communicator);
            rval = recv.f;
            *rank_of_max = recv.i;
        }
    }
    return (rval);
} // maxReduction

void
IBTK_MPI::maxReduction(float* x, const int n, int* rank_of_max, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        if (rank_of_max == NULL)
        {
            float* send = new float[n];
            memcpy(send, x, n * sizeof(float));
            MPI_Allreduce(send, x, n, MPI_FLOAT, MPI_MAX, communicator);
            delete[] send;
        }
        else
        {
            FloatIntStruct* recv = new FloatIntStruct[n];
            FloatIntStruct* send = new FloatIntStruct[n];
            for (int i = 0; i < n; ++i)
            {
                send[i].f = x[i];
                send[i].i = getRank(communicator);
            }
            MPI_Allreduce(send, recv, n, MPI_FLOAT_INT, MPI_MAXLOC, communicator);
            for (int i = 0; i < n; ++i)
            {
                x[i] = recv[i].f;
                rank_of_max[i] = send[i].i;
            }
            delete[] recv;
            delete[] send;
        }
    }
} // maxReduction

int
IBTK_MPI::maxReduction(const int x, int* rank_of_max, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    int rval = x;

    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    /*
     * If a rank_of_max argument is provided, set it to the current
     * rank of the process.
     */
    if (rank_of_max != NULL)
    {
        *rank_of_max = getRank(communicator);
    }

    if (getNodes(communicator) > 1)
    {
        if (rank_of_max == NULL)
        {
            int send = x;
            MPI_Allreduce(&send, &rval, 1, MPI_INT, MPI_MAX, communicator);
        }
        else
        {
            IntIntStruct recv;
            IntIntStruct send;
            send.j = x;
            send.i = getRank(communicator);
            MPI_Allreduce(&send, &recv, 1, MPI_2INT, MPI_MAXLOC, communicator);
            rval = recv.j;
            *rank_of_max = recv.i;
        }
    }
    return (rval);
} // maxReduction

void
IBTK_MPI::maxReduction(int* x, const int n, int* rank_of_max, IBTK_MPI::comm communicator /* = MPI_COMM_NULL */)
{
    communicator = communicator == MPI_COMM_NULL ? d_communicator : communicator;
    if (getNodes(communicator) > 1)
    {
        if (rank_of_max == NULL)
        {
            int* send = new int[n];
            memcpy(send, x, n * sizeof(int));
            MPI_Allreduce(send, x, n, MPI_INT, MPI_MAX, communicator);
            delete[] send;
        }
        else
        {
            IntIntStruct* recv = new IntIntStruct[n];
            IntIntStruct* send = new IntIntStruct[n];
            for (int i = 0; i < n; ++i)
            {
                send[i].j = x[i];
                send[i].i = getRank(communicator);
            }
            MPI_Allreduce(send, recv, n, MPI_2INT, MPI_MAXLOC, communicator);
            for (int i = 0; i < n; ++i)
            {
                x[i] = recv[i].j;
                rank_of_max[i] = send[i].i;
            }
            delete[] recv;
            delete[] send;
        }
    }
} // maxReduction

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

} // namespace IBTK
