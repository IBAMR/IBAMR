// Filename: SpringForceSpecFactory.C
// Created on 14 Jul 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)
// Last modified: <25.Oct.2006 18:30:27 boyce@bigboy.nyconnect.com>

#include "SpringForceSpecFactory.h"

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
#include <ibamr/SpringForceSpec.h>
#include <ibamr/StashableManager.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int SpringForceSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

SpringForceSpecFactory::SpringForceSpecFactory()
{
    setStashableID(StashableManager::getUnregisteredID());
    return;
}// SpringForceSpecFactory

SpringForceSpecFactory::~SpringForceSpecFactory()
{
    // intentionally blank
    return;
}// ~SpringForceSpecFactory

int
SpringForceSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
SpringForceSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

SAMRAI::tbox::Pointer<Stashable>
SpringForceSpecFactory::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int num_links;
    stream.unpack(&num_links,1);
    vector<int> dst_idxs(num_links);
    vector<double> stiffnesses(num_links), rest_lengths(num_links);
    stream.unpack(&dst_idxs[0],num_links);
    stream.unpack(&stiffnesses[0],num_links);
    stream.unpack(&rest_lengths[0],num_links);

    return new SpringForceSpec(dst_idxs,stiffnesses,rest_lengths);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::SpringForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
