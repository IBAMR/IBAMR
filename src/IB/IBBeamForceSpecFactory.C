// Filename: IBBeamForceSpecFactory.C
// Last modified: <01.Sep.2008 13:46:48 boyce@dm-linux.maths.gla.ac.uk>
// Created on 22 Mar 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBBeamForceSpecFactory.h"

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
#include <ibamr/IBBeamForceSpec.h>

// IBTK INCLUDES
#include <ibtk/StashableManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBBeamForceSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBBeamForceSpecFactory::IBBeamForceSpecFactory()
{
    setStashableID(IBTK::StashableManager::getUnregisteredID());
    return;
}// IBBeamForceSpecFactory

IBBeamForceSpecFactory::~IBBeamForceSpecFactory()
{
    // intentionally blank
    return;
}// ~IBBeamForceSpecFactory

int
IBBeamForceSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
IBBeamForceSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

SAMRAI::tbox::Pointer<IBTK::Stashable>
IBBeamForceSpecFactory::unpackStream(
    SAMRAI::tbox::AbstractStream& stream,
    const SAMRAI::hier::IntVector<NDIM>& offset)
{
    int num_beams;
    stream.unpack(&num_beams,1);
    int master_idx;
    stream.unpack(&master_idx,1);
    std::vector<int> tmp_neighbor_idxs(2*num_beams);
    stream.unpack(&tmp_neighbor_idxs[0],2*num_beams);
    std::vector<std::pair<int,int> > neighbor_idxs(num_beams);
    for (int k = 0; k < num_beams; ++k)
    {
        neighbor_idxs[k].first  = tmp_neighbor_idxs[2*k  ];
        neighbor_idxs[k].second = tmp_neighbor_idxs[2*k+1];
    }
    std::vector<double> bend_rigidities(num_beams);
    stream.unpack(&bend_rigidities[0],num_beams);
    std::vector<double> mesh_dependent_curvatures(num_beams);
    stream.unpack(&mesh_dependent_curvatures[0],num_beams);
    return new IBBeamForceSpec(master_idx,neighbor_idxs,bend_rigidities,mesh_dependent_curvatures);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBBeamForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
