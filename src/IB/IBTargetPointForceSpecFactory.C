// Filename: IBTargetPointForceSpecFactory.C
// Last modified: <27.Jun.2010 15:31:25 griffith@griffith-macbook-pro.local>
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
#include <ibamr/namespaces.h>

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

Pointer<Stashable>
IBTargetPointForceSpecFactory::unpackStream(
    AbstractStream& stream,
    const IntVector<NDIM>& offset)
{
    int mastr_idx;
    stream.unpack(&mastr_idx,1);
    double kappa_target;
    stream.unpack(&kappa_target,1);
    double eta_target;
    stream.unpack(&eta_target,1);
    std::vector<double> X_target(NDIM);
    stream.unpack(&X_target[0],NDIM);
    return new IBTargetPointForceSpec(mastr_idx,kappa_target,eta_target,X_target);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBTargetPointForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
