// Filename: IBRodForceSpecFactory.C
// Last modified: <23.Jun.2010 16:42:20 griffith@boyce-griffiths-mac-pro.local>
// Created on 23 Jun 2010 by Boyce Griffith (griffith@boyce-griffiths-mac-pro.local)

#include "IBRodForceSpecFactory.h"

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
#include <ibamr/IBRodForceSpec.h>

// IBTK INCLUDES
#include <ibtk/StashableManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBRodForceSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBRodForceSpecFactory::IBRodForceSpecFactory()
{
    setStashableID(IBTK::StashableManager::getUnregisteredID());
    return;
}// IBRodForceSpecFactory

IBRodForceSpecFactory::~IBRodForceSpecFactory()
{
    // intentionally blank
    return;
}// ~IBRodForceSpecFactory

int
IBRodForceSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
IBRodForceSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

SAMRAI::tbox::Pointer<IBTK::Stashable>
IBRodForceSpecFactory::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int num_rods;
    stream.unpack(&num_rods,1);
    int master_idx;
    stream.unpack(&master_idx,1);
    std::vector<int> next_idxs(num_rods);
    stream.unpack(&next_idxs[0],num_rods);
    std::vector<std::vector<double> > material_params(num_rods,std::vector<double>(10));
    for (int n = 0; n < num_rods; ++n)
    {
        stream.unpack(&material_params[n][0],10);
    }
    return new IBRodForceSpec(master_idx,next_idxs,material_params);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBRodForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
