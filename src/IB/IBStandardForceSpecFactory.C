// Filename: IBStandardForceSpecFactory.C
// Last modified: <19.Mar.2007 21:17:51 griffith@box221.cims.nyu.edu>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

#include "IBStandardForceSpecFactory.h"

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
#include <ibamr/IBStandardForceSpec.h>
#include <ibamr/StashableManager.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBStandardForceSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStandardForceSpecFactory::IBStandardForceSpecFactory()
{
    setStashableID(StashableManager::getUnregisteredID());
    return;
}// IBStandardForceSpecFactory

IBStandardForceSpecFactory::~IBStandardForceSpecFactory()
{
    // intentionally blank
    return;
}// ~IBStandardForceSpecFactory

int
IBStandardForceSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
IBStandardForceSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

SAMRAI::tbox::Pointer<Stashable>
IBStandardForceSpecFactory::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int num_links;
    stream.unpack(&num_links,1);
    vector<int> dst_idxs(num_links), force_fcn_idxs(num_links);
    vector<double> stiffnesses(num_links), rest_lengths(num_links);
    stream.unpack(&dst_idxs[0],num_links);
    stream.unpack(&force_fcn_idxs[0],num_links);
    stream.unpack(&stiffnesses[0],num_links);
    stream.unpack(&rest_lengths[0],num_links);
    vector<double> X_target(NDIM);
    double kappa_target;
    stream.unpack(&X_target[0],NDIM);
    stream.unpack(&kappa_target,1);
    return new IBStandardForceSpec(dst_idxs,force_fcn_idxs,stiffnesses,rest_lengths,X_target,kappa_target);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBStandardForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
