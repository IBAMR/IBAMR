// Filename: StashableStream.I
// Created on 14 Jun 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)
// Last modified: <07.Oct.2006 22:03:00 boyce@bigboy.nyconnect.com>

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/Utilities.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
void
StashableStream::useXDR(
    const bool flag)
{
#ifdef HAVE_XDR
    s_use_xdr_translation = flag;
#else
    if (flag == true)
        TBOX_ERROR("StashableStream: Error in XDR was not used in this compile...\n");
#endif
}// useXDF

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
StashableStream::setCurrentIndex(const int index)
{
    d_buffer_index = index;

    // Is this broken?
    // If xdr is not used how is d_buffer_index moved?

#ifdef HAVE_XDR
    if (d_use_xdr)
    {
#ifndef LACKS_PROPER_XDR_HEADER
        xdr_setpos(&d_xdr_stream, d_buffer_index);
#else
        (*(void(*)(XDR*,unsigned int))(d_xdr_stream.x_ops->x_setpostn))(
            &d_xdr_stream, d_buffer_index);
#endif
    }
#endif
    return;
}// setCurrentIndex

inline
void
StashableStream::resetIndex()
{
    setCurrentIndex(0);
    return;
}// resetIndex

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
