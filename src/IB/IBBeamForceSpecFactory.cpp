// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBBeamForceSpec.h"

#include "ibtk/Streamable.h"
#include "ibtk/StreamableManager.h"

#include "tbox/AbstractStream.h"
#include "tbox/Pointer.h"

#include <memory>
#include <vector>

#include "ibamr/namespaces.h" // IWYU pragma: keep

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
    stream.unpack(tmp_neighbor_idxs.data(), 2 * num_beams);
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
