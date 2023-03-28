// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_IBBeamForceSpec_inl
#define included_IBAMR_IBBeamForceSpec_inl

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/IBBeamForceSpec.h"

#include "ibtk/StreamableManager.h"

#include "tbox/PIO.h"
#include "tbox/Utilities.h"

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

inline bool
IBBeamForceSpec::getIsRegisteredWithStreamableManager()
{
    return (STREAMABLE_CLASS_ID != IBTK::StreamableManager::getUnregisteredID());
} // getIsRegisteredWithStreamableManager

/////////////////////////////// PUBLIC ///////////////////////////////////////

inline IBBeamForceSpec::IBBeamForceSpec(const unsigned int num_beams)
    : d_neighbor_idxs(num_beams), d_bend_rigidities(num_beams), d_mesh_dependent_curvatures(num_beams)
{
#if !defined(NDEBUG)
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBBeamForceSpec::IBBeamForceSpec():\n"
                   << "  must call IBBeamForceSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBBeamForceSpec objects.\n");
    }
#endif
    return;
} // IBBeamForceSpec

inline IBBeamForceSpec::IBBeamForceSpec(const int master_idx,
                                        const std::vector<NeighborIdxs>& neighbor_idxs,
                                        const std::vector<double>& bend_rigidities,
                                        const std::vector<IBTK::Vector>& mesh_dependent_curvatures)
    : d_master_idx(master_idx),
      d_neighbor_idxs(neighbor_idxs),
      d_bend_rigidities(bend_rigidities),
      d_mesh_dependent_curvatures(mesh_dependent_curvatures)
{
#if !defined(NDEBUG)
    const size_t num_beams = d_neighbor_idxs.size();
    TBOX_ASSERT(num_beams == d_bend_rigidities.size());
    TBOX_ASSERT(num_beams == d_mesh_dependent_curvatures.size());
    if (!getIsRegisteredWithStreamableManager())
    {
        TBOX_ERROR("IBBeamForceSpec::IBBeamForceSpec():\n"
                   << "  must call IBBeamForceSpec::registerWithStreamableManager() before\n"
                   << "  creating any IBBeamForceSpec objects.\n");
    }
#endif
    return;
} // IBBeamForceSpec

inline IBBeamForceSpec::~IBBeamForceSpec()
{
    // intentionally blank
    return;
} // ~IBBeamForceSpec

inline unsigned int
IBBeamForceSpec::getNumberOfBeams() const
{
    const auto num_beams = static_cast<unsigned int>(d_neighbor_idxs.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(num_beams == d_bend_rigidities.size());
    TBOX_ASSERT(num_beams == d_mesh_dependent_curvatures.size());
#endif
    return num_beams;
} // getNumberOfBeams

inline const int&
IBBeamForceSpec::getMasterNodeIndex() const
{
    return d_master_idx;
} // getMasterNodeIndex

inline int&
IBBeamForceSpec::getMasterNodeIndex()
{
    return d_master_idx;
} // getMasterNodeIndex

inline const std::vector<IBBeamForceSpec::NeighborIdxs>&
IBBeamForceSpec::getNeighborNodeIndices() const
{
    return d_neighbor_idxs;
} // getNeighborNodeIndices

inline std::vector<IBBeamForceSpec::NeighborIdxs>&
IBBeamForceSpec::getNeighborNodeIndices()
{
    return d_neighbor_idxs;
} // getNeighborNodeIndices

inline const std::vector<double>&
IBBeamForceSpec::getBendingRigidities() const
{
    return d_bend_rigidities;
} // getBendingRigidities

inline std::vector<double>&
IBBeamForceSpec::getBendingRigidities()
{
    return d_bend_rigidities;
} // getBendingRigidities

inline const std::vector<IBTK::Vector>&
IBBeamForceSpec::getMeshDependentCurvatures() const
{
    return d_mesh_dependent_curvatures;
} // getMeshDependentCurvatures

inline std::vector<IBTK::Vector>&
IBBeamForceSpec::getMeshDependentCurvatures()
{
    return d_mesh_dependent_curvatures;
} // getMeshDependentCurvatures

inline int
IBBeamForceSpec::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
} // getStreamableClassID

inline size_t
IBBeamForceSpec::getDataStreamSize() const
{
    const size_t num_beams = d_neighbor_idxs.size();
    return ((2 + 2 * num_beams) * SAMRAI::tbox::AbstractStream::sizeofInt() +
            ((1 + NDIM) * num_beams) * SAMRAI::tbox::AbstractStream::sizeofDouble());
} // getDataStreamSize

inline void
IBBeamForceSpec::packStream(SAMRAI::tbox::AbstractStream& stream)
{
    const auto num_beams = static_cast<unsigned int>(d_neighbor_idxs.size());
#if !defined(NDEBUG)
    TBOX_ASSERT(num_beams == d_bend_rigidities.size());
    TBOX_ASSERT(num_beams == d_mesh_dependent_curvatures.size());
#endif
    std::vector<int> tmp_neighbor_idxs(2 * num_beams);
    for (unsigned int k = 0; k < num_beams; ++k)
    {
        tmp_neighbor_idxs[2 * k] = d_neighbor_idxs[k].first;
        tmp_neighbor_idxs[2 * k + 1] = d_neighbor_idxs[k].second;
    }
    stream << static_cast<int>(num_beams);
    stream.pack(&d_master_idx, 1);
    stream.pack(&tmp_neighbor_idxs[0], 2 * num_beams);
    stream.pack(&d_bend_rigidities[0], 1 * num_beams);
    for (unsigned k = 0; k < num_beams; ++k)
    {
        stream.pack(d_mesh_dependent_curvatures[k].data(), NDIM);
    }
    return;
} // packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_IBBeamForceSpec_inl
