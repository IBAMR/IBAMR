// Filename: IBMovingTargetPointForceSpec.C
// Last modified: <15.Jun.2010 15:46:38 griffith@boyce-griffiths-mac-pro.local>
// Created on 14 Aug 2008 by Boyce Griffith (boyce@dm-linux.maths.gla.ac.uk)

#include "IBMovingTargetPointForceSpec.h"

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
#include <ibamr/IBMovingTargetPointForceSpecFactory.h>

// IBTK INCLUDES
#include <ibtk/StashableManager.h>

// SAMRAI INCLUDES
#include <tbox/SAMRAI_MPI.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

bool IBMovingTargetPointForceSpec::s_registered_factory = false;
int  IBMovingTargetPointForceSpec::s_stashable_id = -1;

void
IBMovingTargetPointForceSpec::registerWithStashableManager()
{
    // We place an MPI barrier here to ensure that all MPI processes actually
    // register the stashable factory with the stashable manager, and to ensure
    // that all processes employ the same stashable id for the
    // IBMovingTargetPointForceSpec object.
    SAMRAI::tbox::SAMRAI_MPI::barrier();
    if (!s_registered_factory)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(s_stashable_id == -1);
#endif
        s_stashable_id = IBTK::StashableManager::getManager()->registerFactory(
            new IBMovingTargetPointForceSpecFactory());
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
template class SAMRAI::tbox::Pointer<IBAMR::IBMovingTargetPointForceSpec>;

//////////////////////////////////////////////////////////////////////////////
