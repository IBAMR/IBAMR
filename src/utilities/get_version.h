#ifndef included_get_version
#define included_get_version

// Filename: get_version.h
// Last modified: <16.Apr.2007 01:51:14 boyce@trasnaform2.local>
// Created on 23 Sep 2006 by Boyce Griffith (boyce@trasnaform2.local)

/////////////////////////////// INCLUDES /////////////////////////////////////

// C++ STDLIB INCLUDES
#include <string>

/////////////////////////////// FUNCTION DEFINITIONS /////////////////////////

namespace IBAMR
{
/*!
 * \brief Return a string indicating the version of the IBAMR library.
 */
std::string
get_version();
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/get_version.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_get_version
