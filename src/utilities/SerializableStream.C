// Filename: StashableStream.C
// Created on 14 Jun 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <02.Oct.2006 11:34:25 boyce@boyce-griffiths-powerbook-g4-15.local>

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "StashableStream.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

bool StashableStream::s_use_xdr_translation = false;

/////////////////////////////// PUBLIC ///////////////////////////////////////

StashableStream::StashableStream(
    const int bytes,
    const StreamMode mode)
{
    d_buffer_size  = bytes;
    d_current_size = 0;
    d_buffer_index = 0;
    d_use_xdr      = s_use_xdr_translation;
    d_buffer       = new char[d_buffer_size];

#ifdef HAVE_XDR
    if (d_use_xdr)
    {
        xdr_op xop = ((mode==StashableStream::Read) ? XDR_DECODE : XDR_ENCODE);
        xdrmem_create(&d_xdr_stream, (caddr_t) d_buffer, d_buffer_size, xop);
        d_xdr_manager.setXDRStream(&d_xdr_stream);
    }
#endif
    return;
}// StashableStream

StashableStream::StashableStream(
    const int bytes,
    const StreamMode mode,
    const bool use_xdr)
{
    d_buffer_size  = bytes;
    d_current_size = 0;
    d_buffer_index = 0;
    d_use_xdr      = use_xdr;
    d_buffer       = new char[d_buffer_size];

#ifdef HAVE_XDR
    if (d_use_xdr)
    {
        xdr_op xop = ((mode==StashableStream::Read) ? XDR_DECODE : XDR_ENCODE);
        xdrmem_create(&d_xdr_stream, (caddr_t) d_buffer, d_buffer_size, xop);
        d_xdr_manager.setXDRStream(&d_xdr_stream);
    }
#endif
    return;
}// StashableStream

StashableStream::StashableStream(
    const void* const buffer,
    const int bytes,
    const StreamMode mode)
{
    d_buffer_size  = bytes;
    d_current_size = 0;
    d_buffer_index = 0;
    d_use_xdr      = s_use_xdr_translation;
    d_buffer       = new char[d_buffer_size];
    memcpy((void*)d_buffer, (void*)buffer, bytes);

#ifdef HAVE_XDR
    if (d_use_xdr)
    {
        xdr_op xop = ((mode==StashableStream::Read) ? XDR_DECODE : XDR_ENCODE);
        xdrmem_create(&d_xdr_stream, (caddr_t) d_buffer, d_buffer_size, xop);
        d_xdr_manager.setXDRStream(&d_xdr_stream);
    }
#endif
    return;
}// StashableStream

StashableStream::StashableStream(
    const void* const buffer,
    const int bytes,
    const StreamMode mode,
    const bool use_xdr)
{
    d_buffer_size  = bytes;
    d_current_size = 0;
    d_buffer_index = 0;
    d_use_xdr      = use_xdr;
    d_buffer       = new char[d_buffer_size];
    memcpy((void*)d_buffer, (void*)buffer, bytes);

#ifdef HAVE_XDR
    if (d_use_xdr)
    {
        xdr_op xop = ((mode==StashableStream::Read) ? XDR_DECODE : XDR_ENCODE);
        xdrmem_create(&d_xdr_stream, (caddr_t) d_buffer, d_buffer_size, xop);
        d_xdr_manager.setXDRStream(&d_xdr_stream);
    }
#endif
    return;
}// StashableStream

StashableStream::~StashableStream()
{
#ifdef HAVE_XDR
    if (d_use_xdr)
    {
#ifndef LACKS_PROPER_XDR_HEADER
        xdr_destroy(&d_xdr_stream);
#else
        if (d_xdr_stream.x_ops->x_destroy)
        {
            (*(void(*)(XDR*))(d_xdr_stream.x_ops->x_destroy))(&d_xdr_stream);
        }
#endif
    }
#endif
    delete[] d_buffer;
    return;
}// ~StashableStream

void
StashableStream::printClassData(
    ostream& os) const
{
    os << "Maximum buffer size = " << d_buffer_size << endl;
    os << "Current buffer size = " << d_current_size << endl;
    os << "Current buffer index = " << d_buffer_index << endl;
    os << "Pointer to buffer data = " << (void *) d_buffer << endl;
    os << "Using XDR translation = " << (d_use_xdr ? "true" : "false") << endl;
    return;
}// printClassData

/*
*************************************************************************
*									*
* Packing/unpacking helper functions and macros.  The member function	*
* getPointerAndAdvanceCursor() returns a pointer to buffer space and	*
* advances internal pointers to reflect the allocated buffers space.	*
* The two macros given below simplify packing and unpacking for the	*
* numerous member functions below.					*
*									*
*************************************************************************
*/

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
            TBOX_ERROR("MessageStream: Stream overrun of buffer...\n");
        }
    }
    return ptr;
}// getPointerAndAdvanceCursor

/*
*************************************************************************
*									*
* The following macros are used by all of the standard data types	*
* except bool for packing and unpacking the data buffer.		*
*									*
*************************************************************************
*/

#ifdef HAVE_XDR

#define PACK(m_data,m_size,m_bytes)                             \
    do                                                          \
    {                                                           \
        void *ptr = getPointerAndAdvanceCursor(m_bytes);        \
        if (d_use_xdr)                                          \
        {                                                       \
            d_xdr_manager.pack(m_data, m_size);                 \
        }                                                       \
        else                                                    \
        {                                                       \
            memcpy(ptr, (void *) m_data, m_bytes);              \
        }                                                       \
    }                                                           \
    while (0)

#define UNPACK(m_data,m_size,m_bytes)                           \
    do                                                          \
    {                                                           \
        void *ptr = getPointerAndAdvanceCursor(m_bytes);        \
        if (d_use_xdr)                                          \
        {                                                       \
            d_xdr_manager.unpack(m_data, m_size);               \
        }                                                       \
        else                                                    \
        {                                                       \
            memcpy((void *) m_data, ptr, m_bytes);              \
        }                                                       \
    }                                                           \
    while (0)

#else

#define PACK(m_data,m_size,m_bytes)                             \
    do                                                          \
    {                                                           \
        void *ptr = getPointerAndAdvanceCursor(m_bytes);        \
        memcpy(ptr, (void *) m_data, m_bytes);                  \
    }                                                           \
    while (0)

#define UNPACK(m_data,m_size,m_bytes)                           \
    do                                                          \
    {                                                           \
        void *ptr = getPointerAndAdvanceCursor(m_bytes);        \
        memcpy((void *) m_data, ptr, m_bytes);                  \
    }                                                           \
    while (0)

#endif


/*
*************************************************************************
*									*
* Packing and unpacking member functions for booleans.  Note that since	*
* the boolean representation is non-standard, boolean arrays are copied	*
* either using XDR or by converting into character arrays.		*
*									*
*************************************************************************
*/

AbstractStream&
StashableStream::operator<<(
    const bool& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

AbstractStream&
StashableStream::operator>>(
    bool& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

void
StashableStream::pack(
    const bool* data,
    const int n)
{
    const int bytes = AbstractStream::sizeofBool(n);
    void *ptr = getPointerAndAdvanceCursor(bytes);
    if (d_use_xdr)
    {
#ifdef HAVE_XDR
        d_xdr_manager.pack(data, n);
#endif
    }
    else
    {
        char *c_ptr = (char *) ptr;
        for (int i = 0; i < n; i++)
        {
            c_ptr[i] = (data[i] ? 1 : 0);
        }
    }
    return;
}// pack

void
StashableStream::unpack(
    bool* data,
    const int n)
{
    const int bytes = AbstractStream::sizeofBool(n);
    void *ptr = getPointerAndAdvanceCursor(bytes);
    if (d_use_xdr)
    {
#ifdef HAVE_XDR
        d_xdr_manager.unpack(data, 1);
#endif
    }
    else
    {
        const char *c_ptr = (const char *) ptr;
        for (int i = 0; i < n; i++)
        {
            data[i] = (c_ptr[i] ? true : false);
        }
    }
    return;
}// unpack

/*
*************************************************************************
*									*
* Packing and unpacking member functions for characters			*
*									*
*************************************************************************
*/

AbstractStream&
StashableStream::operator<<(
    const char& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

AbstractStream&
StashableStream::operator>>(
    char& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

void
StashableStream::pack(
    const char* data,
    const int n)
{
    const int bytes = AbstractStream::sizeofChar(n);
    PACK(data, n, bytes);
    return;
}// pack

void
StashableStream::unpack(
    char* data,
    const int n)
{
    const int bytes = AbstractStream::sizeofChar(n);
    UNPACK(data, n, bytes);
    return;
}// unpack

/*
*************************************************************************
*									*
* Packing and unpacking member functions for double complex		*
*									*
*************************************************************************
*/

AbstractStream&
StashableStream::operator<<(
    const dcomplex& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

AbstractStream&
StashableStream::operator>>(
    dcomplex& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

void
StashableStream::pack(
    const dcomplex* data,
    const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(sizeof(dcomplex) == 2*sizeof(double));
#endif
    const int bytes = AbstractStream::sizeofDoubleComplex(n);
    PACK(data, n, bytes);
    return;
}// pack

void
StashableStream::unpack(
    dcomplex* data,
    const int n)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(sizeof(dcomplex) == 2*sizeof(double));
#endif
    const int bytes = AbstractStream::sizeofDoubleComplex(n);
    UNPACK(data, n, bytes);
    return;
}// unpack

/*
*************************************************************************
*									*
* Packing and unpacking member functions for doubles			*
*									*
*************************************************************************
*/

AbstractStream&
StashableStream::operator<<(
    const double& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

AbstractStream&
StashableStream::operator>>(
    double& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

void
StashableStream::pack(
    const double* data,
    const int n)
{
    const int bytes = AbstractStream::sizeofDouble(n);
    PACK(data, n, bytes);
    return;
}// pack

void
StashableStream::unpack(
    double* data,
    const int n)
{
    const int bytes = AbstractStream::sizeofDouble(n);
    UNPACK(data, n, bytes);
    return;
}// unpack

/*
*************************************************************************
*									*
* Packing and unpacking member functions for floats			*
*									*
*************************************************************************
*/

AbstractStream&
StashableStream::operator<<(
    const float& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

AbstractStream&
StashableStream::operator>>(
    float& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

void
StashableStream::pack(
    const float* data,
    const int n)
{
    const int bytes = AbstractStream::sizeofFloat(n);
    PACK(data, n, bytes);
    return;
}// pack

void
StashableStream::unpack(
    float* data,
    const int n)
{
    const int bytes = AbstractStream::sizeofFloat(n);
    UNPACK(data, n, bytes);
    return;
}// unpack

/*
*************************************************************************
*									*
* Packing and unpacking member functions for integers			*
*									*
*************************************************************************
*/

AbstractStream&
StashableStream::operator<<(
    const int& data)
{
    pack(&data, 1);
    return *this;
}// operator<<

AbstractStream&
StashableStream::operator>>(
    int& data)
{
    unpack(&data, 1);
    return *this;
}// operator>>

void
StashableStream::pack(
    const int* data,
    const int n)
{
    const int bytes = AbstractStream::sizeofInt(n);
    PACK(data, n, bytes);
    return;
}// pack

void
StashableStream::unpack(
    int* data,
    const int n)
{
    const int bytes = AbstractStream::sizeofInt(n);
    UNPACK(data, n, bytes);
    return;
}// unpack

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
