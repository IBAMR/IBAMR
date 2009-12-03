// Filename: IBAnchorPointSpecFactory.C
// Last modified: <18.Aug.2008 13:35:44 boyce@dm-linux.maths.gla.ac.uk>
// Created on 18 Aug 2008 by Boyce Griffith (boyce@dm-linux.maths.gla.ac.uk)

#include "IBAnchorPointSpecFactory.h"

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
#include <ibamr/IBAnchorPointSpec.h>

// IBTK INCLUDES
#include <ibtk/StashableManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBAnchorPointSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBAnchorPointSpecFactory::IBAnchorPointSpecFactory()
{
    setStashableID(IBTK::StashableManager::getUnregisteredID());
    return;
}// IBAnchorPointSpecFactory

IBAnchorPointSpecFactory::~IBAnchorPointSpecFactory()
{
    // intentionally blank
    return;
}// ~IBAnchorPointSpecFactory

int
IBAnchorPointSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
IBAnchorPointSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

SAMRAI::tbox::Pointer<IBTK::Stashable>
IBAnchorPointSpecFactory::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int node_idx;
    stream.unpack(&node_idx,1);
    return new IBAnchorPointSpec(node_idx);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBAnchorPointSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
