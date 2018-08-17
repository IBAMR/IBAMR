// Filename: IBTK_MPI.h
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

#ifndef included_IBTK_MPI
#define included_IBTK_MPI

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "SAMRAI_config.h"

#include "mpi.h"

#include "tbox/Utilities.h"
#include <boost/iterator/iterator_concepts.hpp>

namespace IBTK
{
/**
 * @brief Provides C++ wrapper around MPI routines.
 *
 * The IBTK_MPI struct provides simple interfaces to common MPI routines. All function calls allow a communicator to be
 * used. If a null communicator or none is present, the default communicator set by setCommunicator() is used.
 *
 * Note that this class is a utility class to group function calls in one
 * name space (all calls are to static functions).  Thus, you should never
 * attempt to instantiate a class of type MPI; simply call the functions
 * as static functions using the MPI::function(...) syntax.
 */

struct IBTK_MPI
{
    /**
     * MPI Types
     */
    typedef MPI_Comm comm;
    typedef MPI_Group group;
    typedef MPI_Request request;
    typedef MPI_Status status;

    /**
     * Set the communicator that is used for the MPI communication routines.
     * The default communicator is MPI_COMM_WORLD.
     */
    static void setCommunicator(IBTK_MPI::comm communicator);

    /**
     * Get the current MPI communicator.  The default communicator is
     * MPI_COMM_WORLD.
     */
    static IBTK_MPI::comm getCommunicator();

    /**
     * Get SAMRAI World communicator.
     */
    static IBTK_MPI::comm getSAMRAIWorld();

    /**
     * Return the processor rank (identifier) from 0 through the number of
     * processors minus one.
     */
    static int getRank(IBTK_MPI::comm communicator = MPI_COMM_NULL);

    /**
     * Return the number of processors (nodes).
     */
    static int getNodes(IBTK_MPI::comm communicator = MPI_COMM_NULL);

    /**
     * Perform a global barrier across all processors.
     */
    static void barrier(IBTK_MPI::comm communicator = MPI_COMM_NULL);

    inline static MPI_Datatype mpi_type_id(const double)
    {
        return MPI_DOUBLE;
    }

    inline static MPI_Datatype mpi_type_id(const std::pair<double, int>)
    {
        return MPI_DOUBLE_INT;
    }

    inline static MPI_Datatype mpi_type_id(const int)
    {
        return MPI_INT;
    }

    inline static MPI_Datatype mpi_type_id(const std::pair<int, int>)
    {
        return MPI_2INT;
    }

    inline static MPI_Datatype mpi_type_id(const float)
    {
        return MPI_FLOAT;
    }

    inline static MPI_Datatype mpi_type_id(const std::pair<float, int>)
    {
        return MPI_FLOAT_INT;
    }

    //@{
    /**
     * Perform a min reduction on a data structure of type double, int or float. Each processor contributes an array of
     * values and element-wise min is returned in the same array. If the rank_of_min is not null, the rank of which
     * processor the min is located is stored in the array.
     */
    static void
    minReduction(double* x, const int n = 1, int* rank_of_min = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static void
    minReduction(int* x, const int n = 1, int* rank_of_min = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static void
    minReduction(float* x, const int n = 1, int* rank_of_min = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static double minReduction(double x, int* rank_of_min = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static int minReduction(int x, int* rank_of_min = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static float minReduction(float x, int* rank_of_min = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    //@}

    //@{
    /**
     * Perform a max reduction on a data structure of type double, int or float. Each processor contributes an array of
     * values and element-wise max is returned in the same array. If the rank_of_max is not null, the rank of which
     * processor the max is located is stored in the array.
     */
    static void
    maxReduction(double* x, const int n = 1, int* rank_of_max = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static void
    maxReduction(int* x, const int n = 1, int* rank_of_max = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static void
    maxReduction(float* x, const int n = 1, int* rank_of_max = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static double maxReduction(double x, int* rank_of_max = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static int maxReduction(int x, int* rank_of_max = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static float maxReduction(float x, int* rank_of_max = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    //@}

    //@{
    /**
     * Perform a sum reduction on a data structure of type double, int or float. Each processor contributes an array of
     * values and element-wise sum is returned in the same array.
     */
    static void sumReduction(double* x, const int n = 1, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static void sumReduction(int* x, const int n = 1, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static void sumReduction(float* x, const int n = 1, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static double sumReduction(double x, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static int sumReduction(int x, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static float sumReduction(float x, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    //@}

    /**
     * Perform an all-to-one sum reduction on an integer array.
     * The final result is only available on the root processor.
     */
    static void
    allToOneSumReduction(int* x, const int n, const int root = 0, IBTK_MPI::comm communicator = MPI_COMM_NULL);

    /**
     * Broadcast integer from specified root process to all other processes.
     * All processes other than root, receive a copy of the integer value.
     */
    static int bcast(const int x, const int root, IBTK_MPI::comm communicator = MPI_COMM_NULL);

    /**
     * Broadcast integer array from specified root processor to all other
     * processors.  For the root processor, "array" and "length"
     * are treated as const.
     */
    static void bcast(int* x, int& length, const int root, IBTK_MPI::comm communicator = MPI_COMM_NULL);

    /**
     * Broadcast char array from specified root processor to all other
     * processors.  For the root processor, "array" and "length"
     * are treated as const.
     */
    static void bcast(char* x, int& length, const int root, IBTK_MPI::comm communicator = MPI_COMM_NULL);

    /*!
     * @brief This function sends an MPI message with an integer
     * array to another processer.
     *
     * If the receiving processor knows in advance the length
     * of the array, use "send_length = false;"  otherwise,
     * this processor will first send the length of the array,
     * then send the data.  This call must be paired  with a
     * matching call to IBTK_MPI::recv.
     *
     * @param buf Pointer to integer array buffer with length integers.
     * @param length Number of integers in buf that we want to send.
     * @param receiving_proc_number Receiving processor number.
     * @param send_length Optional boolean argument specifiying if
     * we first need to send a message with the array size.
     * Default value is true.
     * @param tag Optional integer argument specifying an integer tag
     * to be sent with this message.  Default tag is 0.
     */

    static void send(const int* buf,
                     const int length,
                     const int receiving_proc_number,
                     const bool send_length = true,
                     int tag = -1,
                     IBTK_MPI::comm communicator = MPI_COMM_NULL);

    /*!
     * @brief This function sends an MPI message with an array of bytes
     * (MPI_BYTES) to receiving_proc_number.
     *
     * This call must be paired with a matching call to IBTK_MPI::recvBytes.
     *
     * @param buf Void pointer to an array of number_bytes bytes to send.
     * @param number_bytes Integer number of bytes to send.
     * @param receiving_proc_number Receiving processor number.
     */
    static void sendBytes(const void* buf,
                          const int number_bytes,
                          const int receiving_proc_number,
                          IBTK_MPI::comm communicator = MPI_COMM_NULL);

    /*!
     * @brief This function receives an MPI message with an array of
     * max size number_bytes (MPI_BYTES) from any processer.
     *
     * This call must be paired with a matching call to IBTK_MPI::sendBytes.
     *
     * This function returns the processor number of the sender.
     *
     * @param buf Void pointer to a buffer of size number_bytes bytes.
     * @param number_bytes Integer number specifing size of buf in bytes.
     */
    static int recvBytes(void* buf, int number_bytes, IBTK_MPI::comm communicator = MPI_COMM_NULL);

    /*!
     * @brief This function receives an MPI message with an integer
     * array from another processer.
     *
     * If this processor knows in advance the length of the array,
     * use "get_length = false;" otherwise, the sending processor
     * will first send the length of the array, then send the data.
     * This call must be paired with a matching call to IBTK_MPI::send.
     *
     * @param buf Pointer to integer array buffer with capacity of
     * length integers.
     * @param length Maximum number of integers that can be stored in
     * buf.
     * @param sending_proc_number Processor number of sender.
     * @param get_length Optional boolean argument specifiying if
     * we first need to send a message to determine the array size.
     * Default value is true.
     * @param tag Optional integer argument specifying a tag which
     * must be matched by the tag of the incoming message. Default
     * tag is 0.
     */
    static void recv(int* buf,
                     int& length,
                     const int sending_proc_number,
                     const bool get_length = true,
                     int tag = -1,
                     IBTK_MPI::comm communicator = MPI_COMM_NULL);
    //@{
    /**
     * Each processor sends an array of integers or doubles to all other
     * processors; each processor's array may differ in length.
     * The x_out array must be pre-allocated to the correct length
     * (this is a bit cumbersome, but is necessary to avoid th allGather
     * function from allocating memory that is freed elsewhere).
     * To properly preallocate memory, before calling this method, call
     *
     *   size_out = IBTK_MPI::sumReduction(size_in)
     *
     * then allocate the x_out array.
     */
    static void
    allGather(const int* x_in, int size_in, int* x_out, int size_out, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static void allGather(const double* x_in,
                          int size_in,
                          double* x_out,
                          int size_out,
                          IBTK_MPI::comm communicator = MPI_COMM_NULL);
    //@}

    //@{
    /**
     * Each processor sends every other processor an integer or double.
     * The x_out array should be preallocated to a length equal
     * to the number of processors.
     */
    static void allGather(int x_in, int* x_out, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    static void allGather(double x_in, double* x_out, IBTK_MPI::comm communicator = MPI_COMM_NULL);

    //@}

private:
    /**
     * Private functions for MPI calls.
     */
    template <typename T>
    static void
    minReduction(T* x, const int n = 1, int* rank_of_min = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    template <typename T>
    static void
    maxReduction(T* x, const int n = 1, int* rank_of_max = NULL, IBTK_MPI::comm communicator = MPI_COMM_NULL);
    template <typename T>
    static void sumReduction(T* x, const int n = 1, IBTK_MPI::comm communicator = MPI_COMM_NULL);

    /**
     * Performs common functions needed by some of the allToAll methods.
     */
    static void
    allGatherSetup(int size_in, int size_out, int*& rcounts, int*& disps, IBTK_MPI::comm communicator = MPI_COMM_NULL);

    static IBTK_MPI::comm d_communicator;
};

} // namespace IBTK

#endif
