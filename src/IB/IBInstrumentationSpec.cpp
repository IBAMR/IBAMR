// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBInstrumentationSpec.h"

#include "ibtk/IBTK_MPI.h"
#include "ibtk/StreamableFactory.h"
#include "ibtk/StreamableManager.h"

#include <algorithm>
#include <string>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBInstrumentationSpec::STREAMABLE_CLASS_ID = StreamableManager::getUnregisteredID();

std::vector<std::string> IBInstrumentationSpec::s_instrument_names;

void
IBInstrumentationSpec::registerWithStreamableManager()
{
    // We place MPI barriers here to ensure that all MPI processes actually
    // register the factory class with the StreamableManager, and to ensure that
    // all processes employ the same class ID for the IBInstrumentationSpec
    // object.
    IBTK_MPI::barrier();
    if (!getIsRegisteredWithStreamableManager())
    {
#if !defined(NDEBUG)
        TBOX_ASSERT(STREAMABLE_CLASS_ID == StreamableManager::getUnregisteredID());
#endif
        STREAMABLE_CLASS_ID = StreamableManager::getManager()->registerFactory(new IBInstrumentationSpecFactory());
    }
    IBTK_MPI::barrier();
    return;
} // registerWithStreamableManager

void
IBInstrumentationSpec::setInstrumentNames(const std::vector<std::string>& names)
{
    s_instrument_names = names;
    return;
} // getInstrumentNames

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
