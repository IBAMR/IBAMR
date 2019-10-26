// Filename: IBTK_CHKERRQ.h
// Created on 04 Nov 2004 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#ifndef included_IBTK_CHKERRQ
#define included_IBTK_CHKERRQ

/////////////////////////////// INCLUDES /////////////////////////////////////

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
