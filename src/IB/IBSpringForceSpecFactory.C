// Filename: IBSpringForceSpecFactory.C
// Created on 14 Jul 2004 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "IBSpringForceSpecFactory.h"

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
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/StashableManager.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

int IBSpringForceSpecFactory::s_stashable_id = -1;

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBSpringForceSpecFactory::IBSpringForceSpecFactory()
{
    setStashableID(StashableManager::getUnregisteredID());
    return;
}// IBSpringForceSpecFactory

IBSpringForceSpecFactory::~IBSpringForceSpecFactory()
{
    // intentionally blank
    return;
}// ~IBSpringForceSpecFactory

int
IBSpringForceSpecFactory::getStashableID() const
{
    return s_stashable_id;
}// getStashableID

void
IBSpringForceSpecFactory::setStashableID(
    const int stashable_id)
{
    s_stashable_id = stashable_id;
    return;
}// setStashableID

Pointer<Stashable>
IBSpringForceSpecFactory::unpackStream(
    AbstractStream& stream,
    const IntVector<NDIM>& offset)
{
    int num_springs;
    stream.unpack(&num_springs,1);
    int master_idx;
    stream.unpack(&master_idx,1);
    std::vector<int> slave_idxs(num_springs);
    stream.unpack(&slave_idxs[0],num_springs);
    std::vector<int> force_fcn_idxs(num_springs);
    stream.unpack(&force_fcn_idxs[0],num_springs);
    std::vector<double> stiffnesses(num_springs);
    stream.unpack(&stiffnesses[0],num_springs);
    std::vector<double> rest_lengths(num_springs);
    stream.unpack(&rest_lengths[0],num_springs);
    return new IBSpringForceSpec(master_idx,slave_idxs,force_fcn_idxs,stiffnesses,rest_lengths);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBSpringForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
