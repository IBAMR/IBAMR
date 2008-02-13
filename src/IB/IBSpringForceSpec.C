// Filename: IBSpringForceSpec.C
// Last modified: <12.Feb.2008 21:43:57 griffith@box221.cims.nyu.edu>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

#include "IBSpringForceSpec.h"

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
#include <ibamr/IBSpringForceSpecFactory.h>
#include <ibamr/StashableManager.h>

// SAMRAI INCLUDES
#include <tbox/SAMRAI_MPI.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

bool IBSpringForceSpec::s_registered_factory = false;
int  IBSpringForceSpec::s_stashable_id = -1;

void
IBSpringForceSpec::registerWithStashableManager()
{
    // We place an MPI barrier here to ensure that all MPI processes actually
    // register the stashable factory with the stashable manager, and to ensure
    // that all processes employ the same stashable id for the IBSpringForceSpec
    // object.
    SAMRAI::tbox::SAMRAI_MPI::barrier();
    if (!s_registered_factory)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(s_stashable_id == -1);
#endif
        s_stashable_id = StashableManager::getManager()->registerFactory(
            new IBSpringForceSpecFactory());
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
template class SAMRAI::tbox::Pointer<IBAMR::IBSpringForceSpec>;

//////////////////////////////////////////////////////////////////////////////
