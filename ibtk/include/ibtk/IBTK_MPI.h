// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
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

#ifndef included_IBTK_MPI
#define included_IBTK_MPI

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <mpi.h>

#include <type_traits>
#include <utility>
#include <vector>

namespace IBTK
{
inline MPI_Datatype
mpi_type_id(const double)
{
    return MPI_DOUBLE;
}

inline MPI_Datatype
mpi_type_id(const std::pair<double, int>)
{
    return MPI_DOUBLE_INT;
}

inline MPI_Datatype
mpi_type_id(const int)
{
    return MPI_INT;
}

inline MPI_Datatype
mpi_type_id(const std::pair<int, int>)
{
    return MPI_2INT;
}

inline MPI_Datatype
mpi_type_id(const float)
{
    return MPI_FLOAT;
}

inline MPI_Datatype
mpi_type_id(const std::pair<float, int>)
{
    return MPI_FLOAT_INT;
}

inline MPI_Datatype
mpi_type_id(const char)
{
    return MPI_CHAR;
}

inline MPI_Datatype
mpi_type_id(const unsigned int)
{
    return MPI_UNSIGNED;
}

template <typename T>
inline MPI_Datatype
mpi_type_id(const T&)
{
    static_assert(!std::is_same<T, T>::value,
                  "The given type does not have a corresponding MPI_Datatype value. At this time only char, int, "
                  "unsigned int, float, double, std::pair<int, int>, std::pair<int, double>, and "
                  "std::pair<int, float> are supported by IBTK_MPI.");
    return MPI_CHAR;
}

/**
 * @brief Provides C++ wrapper around MPI routines.
 *
 * The IBTK_MPI struct provides simple interfaces to common MPI routines. All function calls allow a user-provided
 * communicator to be used.
 *
 * Note that this class is a utility class to group function calls in one
 * name space (all calls are to static functions).  Thus, you should never
 * attempt to instantiate a class of type IBTK_MPI; simply call the functions
 * as static functions using the IBTK_MPI::function(...) syntax.
 */

struct IBTK_MPI
{
    /**
     * Set the communicator that is used for the MPI communication routines.
     * The default communicator is MPI_COMM_WORLD.
     */
    static void setCommunicator(MPI_Comm communicator);

    /**
     * Get the current MPI communicator.  The default communicator is
     * MPI_COMM_WORLD.
     */
    static MPI_Comm getCommunicator();

    /**
     * Get SAMRAI World communicator.
     */
    static MPI_Comm getSAMRAIWorld();

    /**
     * Return the processor rank (identifier) from 0 through the number of
     * processors minus one.
     */
    static int getRank();

    /**
     * Return the number of processors (nodes).
     */
    static int getNodes();

    /**
     * Perform a global barrier across all processors.
     */
    static void barrier();

    //@{
    /**
     * Perform a min reduction on a data structure of type double, int, or float. Each processor
     * contributes an array of values and element-wise min is returned in the same array. If the rank_of_min is not
     * null, the rank of which processor the min is located is stored in the array.
     */
    template <typename T>
    static T minReduction(T x, int* rank_of_min = nullptr);
    template <typename T>
    static void minReduction(T* x, const int n = 1, int* rank_of_min = nullptr);

    //@}

    //@{
    /**
     * Perform a max reduction on a data structure of type double, int, or float. Each processor
     * contributes an array of values and element-wise max is returned in the same array. If the rank_of_max is not
     * null, the rank of which processor the max is located is stored in the array.
     */
    template <typename T>
    static T maxReduction(T x, int* rank_of_min = nullptr);
    template <typename T>
    static void maxReduction(T* x, const int n = 1, int* rank_of_min = nullptr);
    //@}

    //@{
    /**
     * Perform a sum reduction on a data structure of type double, int, or float. Each processor
     * contributes an array of values and element-wise sum is returned in the same array.
     */
    template <typename T>
    static T sumReduction(T);
    template <typename T>
    static void sumReduction(T* x, const int n = 1);
    //@}

    /**
     * Perform an all-to-one sum reduction on an integer array.
     * The final result is only available on the root processor.
     */
    static void allToOneSumReduction(int* x, const int n, const int root = 0);

    //@{
    /**
     * Broadcast integer array from specified root processor to all other
     * processors.  For the root processor, "array" and "length"
     * are treated as const.
     */
    template <typename T>
    static T bcast(const T x, const int root);
    template <typename T>
    static void bcast(T* x, int& length, const int root);
    //@}

    /*!
     * @brief This function sends an MPI message with an array to another processer.
     *
     * If the receiving processor knows in advance the length
     * of the array, use "send_length = false;"  otherwise,
     * this processor will first send the length of the array,
     * then send the data.  This call must be paired  with a
     * matching call to IBTK_MPI::recv.
     *
     * @param buf Pointer to a valid type array buffer with length integers.
     * @param length Number of integers in buf that we want to send.
     * @param receiving_proc_number Receiving processor number.
     * @param send_length Optional boolean argument specifiying if
     * we first need to send a message with the array size.
     * Default value is true.
     * @param tag Optional integer argument specifying an integer tag
     * to be sent with this message.  Default tag is 0.
     */
    template <typename T>
    static void
    send(const T* buf, const int length, const int receiving_proc_number, const bool send_length = true, int tag = 0);

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
    static void sendBytes(const void* buf, const int number_bytes, const int receiving_proc_number);

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
    static int recvBytes(void* buf, int number_bytes);

    /*!
     * @brief This function receives an MPI message with an array from another processer.
     *
     * If this processor knows in advance the length of the array,
     * use "get_length = false;" otherwise, the sending processor
     * will first send the length of the array, then send the data.
     * This call must be paired with a matching call to IBTK_MPI::send.
     *
     * @param buf Pointer to a valid type array buffer with capacity of
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
    template <typename T>
    static void recv(T* buf, int& length, const int sending_proc_number, const bool get_length = true, int tag = -1);

    //@{
    /**
     * Each processor sends an array of integers or doubles to all other
     * processors; each processor's array may differ in length.
     * The x_out array must be pre-allocated to the correct length
     * (this is a bit cumbersome, but is necessary to avoid the allGather
     * function from allocating memory that is freed elsewhere).
     * To properly preallocate memory, before calling this method, call
     *
     *   size_out = IBTK_MPI::sumReduction(size_in)
     *
     * then allocate the x_out array.
     */
    template <typename T>
    static void allGather(const T* x_in, int size_in, T* x_out, int size_out);
    template <typename T>
    static void allGather(T x_in, T* x_out);

    //@}

private:
    /**
     * Performs common functions needed by some of the allToAll methods.
     */
    static void allGatherSetup(int size_in, int size_out, std::vector<int>& rcounts, std::vector<int>& disps);

    template <typename T>
    static void minMaxReduction(T* x, const int n, int* rank, MPI_Op op);

    static MPI_Comm s_communicator;
};

} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/private/IBTK_MPI-inl.h> // IWYU pragma: keep

#endif
