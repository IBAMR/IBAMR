// Filename: IBAnchorPointSpec.cpp
// Created on 18 Aug 2008 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBAnchorPointSpec.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep

#include "ibtk/StreamableFactory.h"
#include "ibtk/StreamableManager.h"

#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

#include <ostream>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBAnchorPointSpec::STREAMABLE_CLASS_ID = StreamableManager::getUnregisteredID();

void
IBAnchorPointSpec::registerWithStreamableManager()
{
    // We place MPI barriers here to ensure that all MPI processes actually
    // register the factory class with the StreamableManager, and to ensure that
    // all processes employ the same class ID for the IBAnchorPointSpec object.
    SAMRAI_MPI::barrier();
    if (!getIsRegisteredWithStreamableManager())
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(STREAMABLE_CLASS_ID == StreamableManager::getUnregisteredID());
#endif
        STREAMABLE_CLASS_ID = StreamableManager::getManager()->registerFactory(new IBAnchorPointSpecFactory());
    }
    SAMRAI_MPI::barrier();
    return;
} // registerWithStreamableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
