// Filename: IBTargetPointForceSpecFactory.C
// Last modified: <21.Mar.2007 23:30:26 griffith@box221.cims.nyu.edu>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

#include "IBTargetPointForceSpecFactory.h"

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
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/StashableManager.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBTargetPointForceSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBTargetPointForceSpecFactory::IBTargetPointForceSpecFactory()
{
    setStashableID(StashableManager::getUnregisteredID());
    return;
}// IBTargetPointForceSpecFactory

IBTargetPointForceSpecFactory::~IBTargetPointForceSpecFactory()
{
    // intentionally blank
    return;
}// ~IBTargetPointForceSpecFactory

int
IBTargetPointForceSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
IBTargetPointForceSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

SAMRAI::tbox::Pointer<Stashable>
IBTargetPointForceSpecFactory::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    double kappa_target;
    stream.unpack(&kappa_target,1);
    std::vector<double> X_target(NDIM);
    stream.unpack(&X_target[0],NDIM);
    return new IBTargetPointForceSpec(kappa_target,X_target);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBTargetPointForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
