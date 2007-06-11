// Filename: IBBeamForceSpec.I
// Last modified: <11.Jun.2007 17:35:26 griffith@box221.cims.nyu.edu>
// Created on 11 Jun 2007 by Boyce Griffith (griffith@box221.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBBeamForceSpec::getIsRegisteredWithStashableManager()
{
    return s_registered_factory;
}// getIsRegisteredWithStashableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline
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

inline
IBBeamForceSpec::~IBBeamForceSpec()
{
    // intentionally blank
    return;
}// ~IBBeamForceSpec

inline unsigned
IBBeamForceSpec::getNumberOfBeams() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_num_beams == static_cast<int>(d_neighbor_idxs.size()));
    assert(d_num_beams == static_cast<int>(d_bend_rigidities.size()));
#endif
    return d_num_beams;
}// getNumberOfBeams

inline const int&
IBBeamForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
}// getMasterNodeIndex

inline int&
IBBeamForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
}// getMasterNodeIndex

inline const std::vector<IBBeamForceSpec::pair_int_int>&
IBBeamForceSpec::getNeighborNodeIndices() const
{
    return d_neighbor_idxs;
}// getNeighborNodeIndices

inline std::vector<IBBeamForceSpec::pair_int_int>&
IBBeamForceSpec::getNeighborNodeIndices()
{
    return d_neighbor_idxs;
}// getNeighborNodeIndices

inline const std::vector<double>&
IBBeamForceSpec::getBendingRigidities() const
{
    return d_bend_rigidities;
}// getBendingRigidities

inline std::vector<double>&
IBBeamForceSpec::getBendingRigidities()
{
    return d_bend_rigidities;
}// getBendingRigidities

inline int
IBBeamForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBBeamForceSpec::getDataStreamSize() const
{
    return ((2+2*d_num_beams)*SAMRAI::tbox::AbstractStream::sizeofInt() +
            (  1*d_num_beams)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

inline void
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

//////////////////////////////////////////////////////////////////////////////
