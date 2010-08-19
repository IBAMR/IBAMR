// Filename: IBBeamForceSpecFactory.C
// Created on 22 Mar 2007 by Boyce Griffith
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
#include <ibamr/namespaces.h>

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
    setStashableID(StashableManager::getUnregisteredID());
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

Pointer<Stashable>
IBBeamForceSpecFactory::unpackStream(
    AbstractStream& stream,
    const IntVector<NDIM>& offset)
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
    std::vector<std::vector<double> > mesh_dependent_curvatures(num_beams,std::vector<double>(NDIM,0.0));
    for (int k = 0; k < num_beams; ++k)
    {
        stream.unpack(&mesh_dependent_curvatures[k][0],NDIM);
    }
    return new IBBeamForceSpec(master_idx,neighbor_idxs,bend_rigidities,mesh_dependent_curvatures);
}// unpackStream

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

/////////////////////////////// TEMPLATE INSTANTIATION ///////////////////////

#include <tbox/Pointer.C>
template class Pointer<IBAMR::IBBeamForceSpecFactory>;

//////////////////////////////////////////////////////////////////////////////
