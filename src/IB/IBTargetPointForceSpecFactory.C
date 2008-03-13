// Filename: IBTargetPointForceSpecFactory.C
// Last modified: <12.Mar.2008 23:08:21 griffith@box221.cims.nyu.edu>
// Created on 21 Mar 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

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

// IBTK INCLUDES
#include <ibtk/StashableManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBTargetPointForceSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBTargetPointForceSpecFactory::IBTargetPointForceSpecFactory()
{
    setStashableID(IBTK::StashableManager::getUnregisteredID());
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

SAMRAI::tbox::Pointer<IBTK::Stashable>
IBTargetPointForceSpecFactory::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int mastr_idx;
    stream.unpack(&mastr_idx,1);
    double kappa_target;
    stream.unpack(&kappa_target,1);
    std::vector<double> X_target(NDIM);
    stream.unpack(&X_target[0],NDIM);
    return new IBTargetPointForceSpec(mastr_idx,kappa_target,X_target);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBTargetPointForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
