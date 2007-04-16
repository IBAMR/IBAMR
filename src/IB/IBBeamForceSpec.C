// Filename: IBBeamForceSpec.C
// Last modified: <16.Apr.2007 05:47:48 boyce@bigboy.nyconnect.com>
// Created on 22 Mar 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

#include "IBBeamForceSpec.h"

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
#include <ibamr/IBBeamForceSpecFactory.h>
#include <ibamr/StashableManager.h>

// SAMRAI INCLUDES
#include <tbox/MPI.h>
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

bool IBBeamForceSpec::s_registered_factory = false;
int  IBBeamForceSpec::s_stashable_id = -1;

void
IBBeamForceSpec::registerWithStashableManager()
{
    // We place an MPI barrier here to ensure that all MPI processes actually
    // register the stashable factory with the stashable manager, and to ensure
    // that all processes employ the same stashable id for the IBBeamForceSpec
    // object.
    SAMRAI::tbox::MPI::barrier();
    if (!s_registered_factory)
    {
#ifdef DEBUG_CHECK_ASSERTIONS
        assert(s_stashable_id == -1);
#endif
        s_stashable_id = StashableManager::getManager()->registerFactory(
            new IBBeamForceSpecFactory());
        s_registered_factory = true;
    }
    return;
}// registerWithStashableManager

bool
IBBeamForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBBeamForceSpec::IBBeamForceSpec(
    const int master_idx,
    const std::vector<pair_int_int>& neighbor_idxs,
    const std::vector<double>& bend_rigidities)
    : d_num_beams(static_cast<int>(neighbor_idxs.size())),
      d_master_idx(master_idx),
      d_neighbor_idxs(neighbor_idxs),
      d_bend_rigidities(bend_rigidities)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_beams == static_cast<int>(d_neighbor_idxs.size()));
    assert(d_num_beams == static_cast<int>(d_bend_rigidities.size()));
#endif
    if (!s_registered_factory)
    {
        TBOX_ERROR("IBBeamForceSpec::IBBeamForceSpec():\n"
                   << "  must call IBBeamForceSpec::registerWithStashableManager() before\n"
                   << "  creating any IBBeamForceSpec objects.\n");
    }
    return;
}// IBBeamForceSpec

IBBeamForceSpec::~IBBeamForceSpec()
{
    // intentionally blank
    return;
}// ~IBBeamForceSpec

unsigned
IBBeamForceSpec::getNumberOfBeams() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_beams == static_cast<int>(d_neighbor_idxs.size()));
    assert(d_num_beams == static_cast<int>(d_bend_rigidities.size()));
#endif
    return d_num_beams;
}// getNumberOfBeams

const int&
IBBeamForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

int&
IBBeamForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

const std::vector<IBBeamForceSpec::pair_int_int>&
IBBeamForceSpec::getNeighborNodeIndices() const
{
    return d_neighbor_idxs;
}// getNeighborNodeIndices

std::vector<IBBeamForceSpec::pair_int_int>&
IBBeamForceSpec::getNeighborNodeIndices()
{
    return d_neighbor_idxs;
}// getNeighborNodeIndices

const std::vector<double>&
IBBeamForceSpec::getBendingRigidities() const
{
    return d_bend_rigidities;
}// getBendingRigidities

std::vector<double>&
IBBeamForceSpec::getBendingRigidities()
{
    return d_bend_rigidities;
}// getBendingRigidities

int
IBBeamForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

size_t
IBBeamForceSpec::getDataStreamSize() const
{
    return ((2+2*d_num_beams)*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (  1*d_num_beams)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

void
IBBeamForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_beams == static_cast<int>(d_neighbor_idxs.size()));
    assert(d_num_beams == static_cast<int>(d_bend_rigidities.size()));
#endif
    std::vector<int> tmp_neighbor_idxs(2*d_num_beams);
    for (int k = 0; k < d_num_beams; ++k)
    {
        tmp_neighbor_idxs[2*k  ] = d_neighbor_idxs[k].first;
        tmp_neighbor_idxs[2*k+1] = d_neighbor_idxs[k].second;
    }
    stream.pack(&d_num_beams,1);
    stream.pack(&d_master_idx,1);
    stream.pack(&tmp_neighbor_idxs[0],2*d_num_beams);
    stream.pack(&d_bend_rigidities[0],1*d_num_beams);
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::IBBeamForceSpec>;

//////////////////////////////////////////////////////////////////////////////
