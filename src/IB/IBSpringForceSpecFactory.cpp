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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include "ibamr/IBSpringForceSpec.h"

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

IBSpringForceSpec::Factory::Factory()
{
    setStreamableClassID(StreamableManager::getUnregisteredID());
    return;
} // Factory

int
IBSpringForceSpec::Factory::getStreamableClassID() const
{
    return STREAMABLE_CLASS_ID;
} // getStreamableClassID

void
IBSpringForceSpec::Factory::setStreamableClassID(const int class_id)
{
    STREAMABLE_CLASS_ID = class_id;
    return;
} // setStreamableClassID

Pointer<Streamable>
IBSpringForceSpec::Factory::unpackStream(AbstractStream& stream, const IntVector<NDIM>& /*offset*/)
{
    int num_springs;
    stream.unpack(&num_springs, 1);
    Pointer<IBSpringForceSpec> ret_val = new IBSpringForceSpec(num_springs);
    stream.unpack(&ret_val->d_master_idx, 1);
    stream.unpack(&ret_val->d_slave_idxs[0], num_springs);
    stream.unpack(&ret_val->d_force_fcn_idxs[0], num_springs);
    for (int k = 0; k < num_springs; ++k)
    {
        int num_parameters;
        stream.unpack(&num_parameters);
        ret_val->d_parameters[k].resize(num_parameters);
        stream.unpack(&ret_val->d_parameters[k][0], num_parameters);
    }
    return ret_val;
} // unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
