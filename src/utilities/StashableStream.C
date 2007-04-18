// Filename: StashableStream.C
// Last modified: <17.Apr.2007 21:27:10 griffith@box221.cims.nyu.edu>
// Created on 14 Jun 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "StashableStream.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StashableStream::StashableStream(
    const int bytes,
    const StreamMode mode)
    : d_buffer_size(bytes),
      d_current_size(0),
      d_buffer_index(0),
      d_buffer(new char[d_buffer_size])
{
    // intentionally blank
    return;
}// StashableStream

StashableStream::StashableStream(
    const void* const buffer,
    const int bytes,
    const StreamMode mode)
    : d_buffer_size(bytes),
      d_current_size(0),
      d_buffer_index(0),
      d_buffer(new char[d_buffer_size])
{
    memcpy(static_cast<void*>(d_buffer), buffer, bytes);
    return;
}// StashableStream

StashableStream::~StashableStream()
{
    delete[] d_buffer;
    return;
}// ~StashableStream

/*
*************************************************************************
*									*
* Output class data.                            			*
*									*
*************************************************************************
*/

void
StashableStream::printClassData(
    std::ostream& os) const
{
    os << "\nStashableStream::printClassData..." << endl;
    os << "this = " << const_cast<StashableStream*>(this) << endl;
    os << "d_buffer_size = " << d_buffer_size << endl;
    os << "d_current_size = " << d_current_size << endl;
    os << "d_buffer_index = " << d_buffer_index << endl;
    os << "d_buffer = " << static_cast<void*>(d_buffer) << endl;
    return;
}// printClassData

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

//////////////////////////////////////////////////////////////////////////////
