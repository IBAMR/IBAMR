// Filename: IBRodForceSpec.C
// Last modified: <27.Jun.2010 15:30:18 griffith@griffith-macbook-pro.local>
// Created on 23 Jun 2010 by Boyce Griffith (griffith@boyce-griffiths-mac-pro.local)

#include "IBRodForceSpec.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// IBAMR INCLUDES
#include <ibamr/IBRodForceSpecFactory.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/StashableManager.h>

// SAMRAI INCLUDES
#include <tbox/SAMRAI_MPI.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

bool IBRodForceSpec::s_registered_factory = false;
int  IBRodForceSpec::s_stashable_id = -1;

void
IBRodForceSpec::registerWithStashableManager()
{
    // We place an MPI barrier here to ensure that all MPI processes actually
    // register the stashable factory with the stashable manager, and to ensure
    // that all processes employ the same stashable id for the IBRodForceSpec
    // object.
    SAMRAI_MPI::barrier();
    if (!s_registered_factory)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(s_stashable_id == -1);
#endif
        s_stashable_id = StashableManager::getManager()->registerFactory(
            new IBRodForceSpecFactory());
        s_registered_factory = true;
    }
    return;
}// registerWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBRodForceSpec>;

//////////////////////////////////////////////////////////////////////////////
