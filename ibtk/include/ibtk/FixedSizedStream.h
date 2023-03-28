// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_FixedSizedStream
#define included_IBTK_FixedSizedStream

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "tbox/AbstractStream.h"
#include "tbox/Complex.h"

#include <vector>

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
    ~FixedSizedStream() = default;

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
    SAMRAI::tbox::AbstractStream& operator<<(const bool& data) override;

    /*!
     * Remove a single bool from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(bool& data) override;

    /*!
     * Pack an array of bools into the message stream.
     */
    void pack(const bool* data, int n = 1) override;

    /*!
     * Remove an array of bools from the message stream.
     */
    void unpack(bool* data, int n = 1) override;

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
    SAMRAI::tbox::AbstractStream& operator<<(const char& data) override;

    /*!
     * Remove a single char from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(char& data) override;

    /*!
     * Pack an array of chars into the message stream.
     */
    void pack(const char* data, int n = 1) override;

    /*!
     * Remove an array of chars from the message stream.
     */
    void unpack(char* data, int n = 1) override;

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
    SAMRAI::tbox::AbstractStream& operator<<(const dcomplex& data) override;

    /*!
     * Remove a single double complex from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(dcomplex& data) override;

    /*!
     * Pack an array of double complex into the message stream.
     */
    void pack(const dcomplex* data, int n = 1) override;

    /*!
     * Remove an array of double complex from the message stream.
     */
    void unpack(dcomplex* data, int n = 1) override;

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
    SAMRAI::tbox::AbstractStream& operator<<(const double& data) override;

    /*!
     * Remove a single double from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(double& data) override;

    /*!
     * Pack an array of doubles into the message stream.
     */
    void pack(const double* data, int n = 1) override;

    /*!
     * Remove an array of doubles from the message stream.
     */
    void unpack(double* data, int n = 1) override;

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
    SAMRAI::tbox::AbstractStream& operator<<(const float& data) override;

    /*!
     * Remove a single float from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(float& data) override;

    /*!
     * Pack an array of floats into the message stream.
     */
    void pack(const float* data, int n = 1) override;

    /*!
     * Remove an array of floats from the message stream.
     */
    void unpack(float* data, int n = 1) override;

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
    SAMRAI::tbox::AbstractStream& operator<<(const int& data) override;

    /*!
     * Remove a single integer from the message stream.
     */
    SAMRAI::tbox::AbstractStream& operator>>(int& data) override;

    /*!
     * Pack an array of integers into the message stream.
     */
    void pack(const int* data, int n = 1) override;

    /*!
     * Remove an array of integers from the message stream.
     */
    void unpack(int* data, int n = 1) override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FixedSizedStream() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FixedSizedStream(const FixedSizedStream& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FixedSizedStream& operator=(const FixedSizedStream& that) = delete;

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
    int d_current_size = 0;

    /*
     * The index of the first free element in the buffer.
     */
    int d_buffer_index = 0;

    /*
     * The buffer.
     */
    std::vector<char> d_buffer;
};
} // namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "ibtk/private/FixedSizedStream-inl.h" // IWYU pragma: keep

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_FixedSizedStream
