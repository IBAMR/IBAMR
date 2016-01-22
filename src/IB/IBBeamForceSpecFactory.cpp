// Filename: IBBeamForceSpecFactory.cpp
// Created on 22 Mar 2007 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#include <utility>
#include <vector>

#include "Eigen/Core"
#include "ibamr/IBBeamForceSpec.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/Streamable.h"
#include "ibtk/StreamableManager.h"
#include "tbox/AbstractStream.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class IntVector;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// PUBLIC ///////////////////////////////////////

IBBeamForceSpec::Factory::Factory()
{
    setStreamableClassID(StreamableManager::getUnregisteredID());
    return;
} // Factory

IBBeamForceSpec::Factory::~Factory()
{
    // intentionally blank
    return;
} // ~Factory

int
IBBeamForceSpec::Factory::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
} // getStreamableClassID

void
IBBeamForceSpec::Factory::setStreamableClassID(const int class_id)
{
    STREAMABLE_CLASS_ID = class_id;
    return;
} // setStreamableClassID

Pointer<Streamable>
IBBeamForceSpec::Factory::unpackStream(AbstractStream& stream, const IntVector<NDIM>& /*offset*/)
{
    int num_beams;
    stream.unpack(&num_beams, 1);
    Pointer<IBBeamForceSpec> ret_val = new IBBeamForceSpec(num_beams);
    stream.unpack(&ret_val->d_master_idx, 1);
    std::vector<int> tmp_neighbor_idxs(2 * num_beams);
    stream.unpack(&tmp_neighbor_idxs[0], 2 * num_beams);
    for (int k = 0; k < num_beams; ++k)
    {
        ret_val->d_neighbor_idxs[k].first = tmp_neighbor_idxs[2 * k];
        ret_val->d_neighbor_idxs[k].second = tmp_neighbor_idxs[2 * k + 1];
    }
    stream.unpack(&ret_val->d_bend_rigidities[0], num_beams);
    for (int k = 0; k < num_beams; ++k)
    {
        stream.unpack(ret_val->d_mesh_dependent_curvatures[k].data(), NDIM);
    }
    return ret_val;
} // unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
