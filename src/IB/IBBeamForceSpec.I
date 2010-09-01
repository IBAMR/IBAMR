// Filename: IBBeamForceSpec.I
// Created on 11 Jun 2007 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

// SAMRAI INCLUDES
#include <tbox/PIO.h>
#include <tbox/Utilities.h>

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
    const std::vector<NeighborIdxs>& neighbor_idxs,
    const std::vector<double>& bend_rigidities,
    const std::vector<std::vector<double> >& mesh_dependent_curvatures)
    : d_master_idx(master_idx),
      d_neighbor_idxs(neighbor_idxs),
      d_bend_rigidities(bend_rigidities),
      d_mesh_dependent_curvatures(mesh_dependent_curvatures)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    const int num_beams = d_neighbor_idxs.size();
    TBOX_ASSERT(num_beams == int(d_bend_rigidities.size()));
    TBOX_ASSERT(num_beams == int(d_mesh_dependent_curvatures.size()));
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
    const int num_beams = d_neighbor_idxs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_beams == int(d_bend_rigidities.size()));
    TBOX_ASSERT(num_beams == int(d_mesh_dependent_curvatures.size()));
#endif
    return num_beams;
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

inline const std::vector<IBBeamForceSpec::NeighborIdxs>&
IBBeamForceSpec::getNeighborNodeIndices() const
{
    return d_neighbor_idxs;
}// getNeighborNodeIndices

inline std::vector<IBBeamForceSpec::NeighborIdxs>&
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

inline const std::vector<std::vector<double> >&
IBBeamForceSpec::getMeshDependentCurvatures() const
{
    return d_mesh_dependent_curvatures;
}// getMeshDependentCurvatures

inline std::vector<std::vector<double> >&
IBBeamForceSpec::getMeshDependentCurvatures()
{
    return d_mesh_dependent_curvatures;
}// getMeshDependentCurvatures

inline int
IBBeamForceSpec::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

inline size_t
IBBeamForceSpec::getDataStreamSize() const
{
    const int num_beams = d_neighbor_idxs.size();
    return ((2+2*num_beams     )*SAMRAI::tbox::AbstractStream::sizeofInt() +
            ((1+NDIM)*num_beams)*SAMRAI::tbox::AbstractStream::sizeofDouble());
}// getDataStreamSize

inline void
IBBeamForceSpec::packStream(
    SAMRAI::tbox::AbstractStream& stream)
{
    const int num_beams = d_neighbor_idxs.size();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(num_beams == int(d_bend_rigidities.size()));
    TBOX_ASSERT(num_beams == int(d_mesh_dependent_curvatures.size()));
#endif
    std::vector<int> tmp_neighbor_idxs(2*num_beams);
    for (int k = 0; k < num_beams; ++k)
    {
        tmp_neighbor_idxs[2*k  ] = d_neighbor_idxs[k].first;
        tmp_neighbor_idxs[2*k+1] = d_neighbor_idxs[k].second;
    }
    stream.pack(&num_beams,1);
    stream.pack(&d_master_idx,1);
    stream.pack(&tmp_neighbor_idxs[0],2*num_beams);
    stream.pack(&d_bend_rigidities[0],1*num_beams);
    for (int k = 0; k < num_beams; ++k)
    {
        stream.pack(&d_mesh_dependent_curvatures[k][0],NDIM);
    }
    return;
}// packStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
