// Filename: IBSpringForceSpecFactory.C
// Last modified: <21.Mar.2007 21:26:07 griffith@box221.cims.nyu.edu>
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

#include "IBSpringForceSpecFactory.h"

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
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/StashableManager.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBSpringForceSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBSpringForceSpecFactory::IBSpringForceSpecFactory()
{
    setStashableID(StashableManager::getUnregisteredID());
    return;
}// IBSpringForceSpecFactory

IBSpringForceSpecFactory::~IBSpringForceSpecFactory()
{
    // intentionally blank
    return;
}// ~IBSpringForceSpecFactory

int
IBSpringForceSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
IBSpringForceSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

SAMRAI::tbox::Pointer<Stashable>
IBSpringForceSpecFactory::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int master_idx, num_springs;
    stream.unpack(&master_idx,1);
    stream.unpack(&num_springs,1);
    vector<int> slave_idxs(num_springs), force_fcn_idxs(num_springs);
    vector<double> stiffnesses(num_springs), rest_lengths(num_springs);
    stream.unpack(&slave_idxs[0],num_springs);
    stream.unpack(&force_fcn_idxs[0],num_springs);
    stream.unpack(&stiffnesses[0],num_springs);
    stream.unpack(&rest_lengths[0],num_springs);
    return new IBSpringForceSpec(master_idx,slave_idxs,force_fcn_idxs,stiffnesses,rest_lengths);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBSpringForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
