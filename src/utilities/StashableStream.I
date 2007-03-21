// Filename: StashableStream.I
// Created on 14 Jun 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <20.Mar.2007 22:21:33 griffith@box221.cims.nyu.edu>

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
void*
StashableStream::getBufferStart()
{
    return static_cast<void*>(d_buffer);
}// getBufferStart

inline
int
StashableStream::getCurrentSize() const
{
    return d_current_size;
}// getCurrentSize

inline
int
StashableStream::getCurrentIndex() const
{
    return d_buffer_index;
}// getCurrentIndex

inline
void
StashableStream::setCurrentIndex(
    const int index)
{
    d_buffer_index = index;
    return;
}// setCurrentIndex

inline
void
StashableStream::resetIndex()
{
    setCurrentIndex(0);
    return;
}// resetIndex

/*
*************************************************************************
*									*
* Packing and unpacking member functions for booleans.  Note that since	*
* the boolean representation is non-standard, boolean arrays are copied	*
* by converting into character arrays.	                                *
*									*
*************************************************************************
*/

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator<<(
    const bool& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator>>(
    bool& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

inline
void
StashableStream::pack(
    const bool* data,
    const int n)
{
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofBool(n);
    void* ptr = getPointerAndAdvanceCursor(bytes);
    char* c_ptr = static_cast<char*>(ptr);
    for (int i = 0; i < n; i++)
    {
        c_ptr[i] = (data[i] ? 1 : 0);
    }
    return;
}// pack

inline
void
StashableStream::unpack(
    bool* data,
    const int n)
{
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofBool(n);
    void* ptr = getPointerAndAdvanceCursor(bytes);
    const char* c_ptr = static_cast<const char*>(ptr);
    for (int i = 0; i < n; i++)
    {
        data[i] = (c_ptr[i] == 1);
    }
    return;
}// unpack

/*
*************************************************************************
*									*
* Packing and unpacking member functions for characters.		*
*									*
*************************************************************************
*/

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator<<(
    const char& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator>>(
    char& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

inline
void
StashableStream::pack(
    const char* data,
    const int n)
{
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofChar(n);
    __pack(data, bytes);
    return;
}// pack

inline
void
StashableStream::unpack(
    char* data,
    const int n)
{
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofChar(n);
    __unpack(data, bytes);
    return;
}// unpack

/*
*************************************************************************
*									*
* Packing and unpacking member functions for double complex.    	*
*									*
*************************************************************************
*/

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator<<(
    const dcomplex& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator>>(
    dcomplex& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

inline
void
StashableStream::pack(
    const dcomplex* data,
    const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(sizeof(dcomplex) == 2*sizeof(double));
#endif
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofDoubleComplex(n);
    __pack(data, bytes);
    return;
}// pack

inline
void
StashableStream::unpack(
    dcomplex* data,
    const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(sizeof(dcomplex) == 2*sizeof(double));
#endif
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofDoubleComplex(n);
    __unpack(data, bytes);
    return;
}// unpack

/*
*************************************************************************
*									*
* Packing and unpacking member functions for doubles.			*
*									*
*************************************************************************
*/

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator<<(
    const double& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator>>(
    double& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

inline
void
StashableStream::pack(
    const double* data,
    const int n)
{
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofDouble(n);
    __pack(data, bytes);
    return;
}// pack

inline
void
StashableStream::unpack(
    double* data,
    const int n)
{
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofDouble(n);
    __unpack(data, bytes);
    return;
}// unpack

/*
*************************************************************************
*									*
* Packing and unpacking member functions for floats.			*
*									*
*************************************************************************
*/

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator<<(
    const float& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator>>(
    float& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

inline
void
StashableStream::pack(
    const float* data,
    const int n)
{
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofFloat(n);
    __pack(data, bytes);
    return;
}// pack

inline
void
StashableStream::unpack(
    float* data,
    const int n)
{
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofFloat(n);
    __unpack(data, bytes);
    return;
}// unpack

/*
*************************************************************************
*									*
* Packing and unpacking member functions for integers.			*
*									*
*************************************************************************
*/

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator<<(
    const int& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

inline
SAMRAI::tbox::AbstractStream&
StashableStream::operator>>(
    int& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

inline
void
StashableStream::pack(
    const int* data,
    const int n)
{
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofInt(n);
    __pack(data, bytes);
    return;
}// pack

inline
void
StashableStream::unpack(
    int* data,
    const int n)
{
    const int bytes = SAMRAI::tbox::AbstractStream::sizeofInt(n);
    __unpack(data, bytes);
    return;
}// unpack

/////////////////////////////// PRIVATE //////////////////////////////////////

/*
*************************************************************************
*									*
* Packing/unpacking helper functions.  The member function              *
* getPointerAndAdvanceCursor() returns a pointer to buffer space and    *
* advances internal pointers to reflect the allocated buffers space.    *
* The two functions given below simplify packing and unpacking for      *
* the * numerous member functions below.                                *
*									*
*************************************************************************
*/

inline
void*
StashableStream::getPointerAndAdvanceCursor(
    const int bytes)
{
    void* ptr = &d_buffer[d_buffer_index];
    d_buffer_index += bytes;
    if (d_buffer_index > d_current_size)
    {
        d_current_size = d_buffer_index;
        if (d_buffer_index > d_buffer_size)
        {
            TBOX_ERROR("StashableStream::getPointerAndAdvanceCursor():\n"
                       << "  buffer overrun." << endl);
        }
    }
    return ptr;
}// getPointerAndAdvanceCursor

template <typename T>
inline
void
StashableStream::__pack(
    const T* const m_data,
    unsigned m_bytes)
{
    void* const ptr = getPointerAndAdvanceCursor(m_bytes);
    memcpy(ptr, static_cast<const void*>(m_data), m_bytes);
    return;
}// _pack

template <typename T>
inline
void
StashableStream::__unpack(
    T* const m_data,
    unsigned m_bytes)
{
    const void* const ptr = getPointerAndAdvanceCursor(m_bytes);
    memcpy(static_cast<void*>(m_data), ptr, m_bytes);
    return;
}// _unpack

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
