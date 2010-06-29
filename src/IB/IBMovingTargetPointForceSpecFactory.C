// Filename: IBMovingTargetPointForceSpecFactory.C
// Last modified: <27.Jun.2010 15:30:12 griffith@griffith-macbook-pro.local>
// Created on 14 Aug 2008 by Boyce Griffith (boyce@dm-linux.maths.gla.ac.uk)

#include "IBMovingTargetPointForceSpecFactory.h"

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
#include <ibamr/IBMovingTargetPointForceSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/StashableManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBMovingTargetPointForceSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBMovingTargetPointForceSpecFactory::IBMovingTargetPointForceSpecFactory()
{
    setStashableID(StashableManager::getUnregisteredID());
    return;
}// IBMovingTargetPointForceSpecFactory

IBMovingTargetPointForceSpecFactory::~IBMovingTargetPointForceSpecFactory()
{
    // intentionally blank
    return;
}// ~IBMovingTargetPointForceSpecFactory

int
IBMovingTargetPointForceSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
IBMovingTargetPointForceSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

Pointer<Stashable>
IBMovingTargetPointForceSpecFactory::unpackStream(
    AbstractStream& stream,
    const IntVector<NDIM>& offset)
{
    int mastr_idx;
    stream.unpack(&mastr_idx,1);
    double kappa_target;
    stream.unpack(&kappa_target,1);
    double eta_target;
    stream.unpack(&eta_target,1);
    int spec_fcn_idx;
    stream.unpack(&spec_fcn_idx,1);
    std::vector<double> periodic_shift(NDIM);
    stream.unpack(&periodic_shift[0],NDIM);
    return new IBMovingTargetPointForceSpec(mastr_idx,kappa_target,eta_target,spec_fcn_idx,periodic_shift);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBMovingTargetPointForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
