// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
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

#ifndef included_IBTK_app_namespaces
#define included_IBTK_app_namespaces

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/namespaces.h>

//////////////////////////////////////////////////////////////////////////////

#ifndef DOXYGEN
/*!
 * Defines "using" declarations for all namespaces used in IBTK.  This header
 * file may be included in application codes, but it MUST NOT be included in any
 * other header (.h) or inline (-inl.h) file in the library.
 */

namespace Eigen
{
}
using namespace Eigen;

namespace std
{
}
using namespace std;

#endif

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_app_namespaces
