// Filename: IBInstrumentationSpec.C
// Last modified: <12.Mar.2008 23:00:45 griffith@box221.cims.nyu.edu>
// Created on 11 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBInstrumentationSpec.h"

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
#include <ibamr/IBInstrumentationSpecFactory.h>

// IBTK INCLUDES
#include <ibtk/StashableManager.h>

// SAMRAI INCLUDES
#include <tbox/SAMRAI_MPI.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

bool IBInstrumentationSpec::s_registered_factory = false;
int  IBInstrumentationSpec::s_stashable_id = -1;

std::vector<std::string> IBInstrumentationSpec::s_instrument_names;

void
IBInstrumentationSpec::registerWithStashableManager()
{
    // We place an MPI barrier here to ensure that all MPI processes actually
    // register the stashable factory with the stashable manager, and to ensure
    // that all processes employ the same stashable id for the IBInstrumentationSpec
    // object.
    SAMRAI::tbox::SAMRAI_MPI::barrier();
    if (!s_registered_factory)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        TBOX_ASSERT(s_stashable_id == -1);
#endif
        s_stashable_id = IBTK::StashableManager::getManager()->registerFactory(
            new IBInstrumentationSpecFactory());
        s_registered_factory = true;
    }
    return;
}// registerWithStashableManager

void
IBInstrumentationSpec::setInstrumentNames(
    const std::vector<std::string>& names)
{
    s_instrument_names = names;
    return;
}// getInstrumentNames

/////////////////////////////// PUBLIC ///////////////////////////////////////

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBInstrumentationSpec>;

//////////////////////////////////////////////////////////////////////////////
