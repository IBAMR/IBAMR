// Filename: FixedSizedStream.h
// Created on 14 Jun 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#ifndef included_FixedSizedStream
#define included_FixedSizedStream

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <vector>

#include "tbox/AbstractStream.h"
#include "tbox/Complex.h"

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class FixedSizedStream provides a fixed-size message buffer used by
 * various communication routines.
 *
 * This class implements the SAMRAI::tbox::AbstractStream interface.  Class
 * FixedSizedStream can packs and unpacks message streams via straight-forward
 * byte copying.
 *
 * \warning This class stores and communicates values in the native binary
 * machine format.  Consequently, this implementation is not suitable for
 * heterogeneous networks, which require a machine-independent storage format
 * such as XDR.
 */
class FixedSizedStream : public SAMRAI::tbox::AbstractStream
{
public:
    /*!
     * Create a message stream of the specified size in bytes.
     */
    FixedSizedStream(int bytes);

    /*!
     * Create a message stream with the specified buffer.
     */
    FixedSizedStream(const void* buffer, int bytes);

    /*!
     * Destructor for a message stream.
     */
    ~FixedSizedStream();

    /*!
     * Return a pointer to the start of the message buffer.
     */
    void* getBufferStart();

    /*!
     * Return a const pointer to the start of the message buffer.
     */
    const void* getBufferStart() const;

    /*!
     * Return the current size of the buffer in bytes.
     */
    int getCurrentSize() const;

    /*!
     * Return the current index into the buffer.
     */
    int getCurrentIndex() const;

    /*!
     * Set the current index into the buffer.  Further packing/unpacking will
     * begin at this new location.
     */
    void setCurrentIndex(int index);

    /*!
     * Reset the index to the beginning of the buffer.  This is the same as
     * setting the buffer index to zero via setCurrentIndex().
     */
    void resetIndex();

    /*!
     * \name Boolean Stream Primitives
     *
     * Pack and unpack booleans into and out of the message stream.
     */
    //\{

    /*!
     * Pack a single bool into the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator<<(const bool& data);

    /*!
     * Remove a single bool from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(bool& data);

    /*!
     * Pack an array of bools into the message stream.
     */
    void pack(const bool* data, int n = 1);

    /*!
     * Remove an array of bools from the message stream.
     */
    void unpack(bool* data, int n = 1);

    //\}

    /*!
     * \name Character Stream Primitives
     *
     * Pack and unpack chars into and out of the message stream.
     */
    //\{

    /*!
     * Pack a single char into the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator<<(const char& data);

    /*!
     * Remove a single char from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(char& data);

    /*!
     * Pack an array of chars into the message stream.
     */
    void pack(const char* data, int n = 1);

    /*!
     * Remove an array of chars from the message stream.
     */
    void unpack(char* data, int n = 1);

    //\}

    /*!
     * \name Double Complex Stream Primitives
     *
     * Pack and unpack double complex into and out of the message stream.
     */
    //\{

    /*!
     * Pack a single double complex into the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator<<(const dcomplex& data);

    /*!
     * Remove a single double complex from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(dcomplex& data);

    /*!
     * Pack an array of double complex into the message stream.
     */
    void pack(const dcomplex* data, int n = 1);

    /*!
     * Remove an array of double complex from the message stream.
     */
    void unpack(dcomplex* data, int n = 1);

    //\}

    /*!
     * \name Double Stream Primitives
     *
     * Pack and unpack doubles into and out of the message stream.
     */
    //\{

    /*!
     * Pack a single double into the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator<<(const double& data);

    /*!
     * Remove a single double from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(double& data);

    /*!
     * Pack an array of doubles into the message stream.
     */
    void pack(const double* data, int n = 1);

    /*!
     * Remove an array of doubles from the message stream.
     */
    void unpack(double* data, int n = 1);

    //\}

    /*!
     * \name Float Stream Primitives
     *
     * Pack and unpack floats into and out of the message stream.
     */
    //\{

    /*!
     * Pack a single float into the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator<<(const float& data);

    /*!
     * Remove a single float from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(float& data);

    /*!
     * Pack an array of floats into the message stream.
     */
    void pack(const float* data, int n = 1);

    /*!
     * Remove an array of floats from the message stream.
     */
    void unpack(float* data, int n = 1);

    //\}

    /*!
     * \name Integer Stream Primitives
     *
     * Pack and unpack integers into and out of the message stream.
     */
    //\{

    /*!
     * Pack a single integer into the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator<<(const int& data);

    /*!
     * Remove a single integer from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(int& data);

    /*!
     * Pack an array of integers into the message stream.
     */
    void pack(const int* data, int n = 1);

    /*!
     * Remove an array of integers from the message stream.
     */
    void unpack(int* data, int n = 1);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FixedSizedStream();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FixedSizedStream(const FixedSizedStream& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FixedSizedStream& operator=(const FixedSizedStream& that);

    /*!
     * \brief Return a pointer to buffer space and advance internal pointers to
     * reflect the allocated buffers space.
     */
    void* getPointerAndAdvanceCursor(int bytes);

    /*!
     * \brief Pack the specified data to the buffer.
     */
    template <typename T>
    void __pack(const T* m_data, unsigned int m_bytes);

    /*!
     * \brief Unpack the specified data from the buffer.
     */
    template <typename T>
    void __unpack(T* m_data, unsigned int m_bytes);

    /*
     * The size of the buffer.
     */
    const int d_buffer_size;

    /*
     * The current size of the buffer, i.e., the number of bytes in the buffer
     * which are currently in use.
     */
    int d_current_size;

    /*
     * The index of the first free element in the buffer.
     */
    int d_buffer_index;

    /*
     * The buffer.
     */
    std::vector<char> d_buffer;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/FixedSizedStream-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_FixedSizedStream
