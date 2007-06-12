// Filename: IBInstrumentationSpecFactory.C
// Last modified: <11.Jun.2007 18:58:41 griffith@box221.cims.nyu.edu>
// Created on 11 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBInstrumentationSpecFactory.h"

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
#include <ibamr/IBInstrumentationSpec.h>
#include <ibamr/StashableManager.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBInstrumentationSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBInstrumentationSpecFactory::IBInstrumentationSpecFactory()
{
    setStashableID(StashableManager::getUnregisteredID());
    return;
}// IBInstrumentationSpecFactory

IBInstrumentationSpecFactory::~IBInstrumentationSpecFactory()
{
    // intentionally blank
    return;
}// ~IBInstrumentationSpecFactory

int
IBInstrumentationSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
IBInstrumentationSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

SAMRAI::tbox::Pointer<Stashable>
IBInstrumentationSpecFactory::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int master_idx;
    stream.unpack(&master_idx,1);
    int meter_idx;
    stream.unpack(&meter_idx,1);
    int node_idx;
    stream.unpack(&node_idx,1);
    return new IBInstrumentationSpec(master_idx,meter_idx,node_idx);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBInstrumentationSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
