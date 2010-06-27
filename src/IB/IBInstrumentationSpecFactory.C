// Filename: IBInstrumentationSpecFactory.C
// Last modified: <27.Jun.2010 15:29:35 griffith@griffith-macbook-pro.local>
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
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/StashableManager.h>

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

Pointer<Stashable>
IBInstrumentationSpecFactory::unpackStream(
    AbstractStream& stream,
    const IntVector<NDIM>& offset)
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
template class Pointer<IBAMR::IBInstrumentationSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
