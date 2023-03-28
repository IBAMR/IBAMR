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

#ifndef included_IBTK_CHKERRQ
#define included_IBTK_CHKERRQ

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/compiler_hints.h"

#include "tbox/Utilities.h"

#include <ostream>

/////////////////////////////// MACRO DEFINITIONS ////////////////////////////

namespace IBTK
{
/*!
 * \brief Throw an error exception from within any C++ source code.
 *
 * This is is similar to the PETSc CHKERRQ(ierr) macro and is designed to be
 * invoked after a call to a PETSc library function.
 */
#define IBTK_CHKERRQ(ierr)                                                                                             \
    if (UNLIKELY(ierr))                                                                                                \
    {                                                                                                                  \
        std::ostringstream tboxos;                                                                                     \
        CHKERRCONTINUE(ierr);                                                                                          \
        SAMRAI::tbox::Utilities::abort(tboxos.str().c_str(), __FILE__, __LINE__);                                      \
    }
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_CHKERRQ
