#ifndef included_get_version
#define included_get_version

// Filename: get_version.h
// Last modified: <24.Oct.2006 14:35:55 boyce@bigboy.nyconnect.com>
// Created on 23 Sep 2006 by Boyce Griffith (boyce@boyce-griffiths-powerbook-g4-15.local)

/////////////////////////////// INCLUDES /////////////////////////////////////

// C++ STDLIB INCLUDES
#include <string>

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBAMR
{
/*!
 * @brief Return a string indicating the version of the IBAMR
 * library.
 */
std::string get_version();
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <stools/get_version.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_get_version
