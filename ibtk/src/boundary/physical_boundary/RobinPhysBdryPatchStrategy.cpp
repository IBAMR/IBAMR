// Filename: RobinPhysBdryPatchStrategy.cpp
// Created on 28 Apr 2012 by Boyce Griffith
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

#include <ostream>
#include <set>
#include <vector>

#include "Box.h"
#include "ComponentSelector.h"
#include "IntVector.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
} // namespace hier
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

RobinPhysBdryPatchStrategy::RobinPhysBdryPatchStrategy() : d_patch_data_indices(), d_bc_coefs(), d_homogeneous_bc(false)
{
    // intentionally blank
    return;
} // RobinPhysBdryPatchStrategy

RobinPhysBdryPatchStrategy::~RobinPhysBdryPatchStrategy()
{
    // intentionally blank
    return;
} // ~RobinPhysBdryPatchStrategy

void
RobinPhysBdryPatchStrategy::setPatchDataIndex(const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
} // setPatchDataIndex

void
RobinPhysBdryPatchStrategy::setPatchDataIndices(const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
} // setPatchDataIndices

void
RobinPhysBdryPatchStrategy::setPatchDataIndices(const ComponentSelector& patch_data_indices)
{
    std::set<int> patch_data_index_set;
    for (int l = 0; l < patch_data_indices.getSize(); ++l)
    {
        if (patch_data_indices.isSet(l))
        {
            const int patch_data_index = l;
            patch_data_index_set.insert(patch_data_index);
        }
    }
    setPatchDataIndices(patch_data_index_set);
    return;
} // setPatchDataIndices

void
RobinPhysBdryPatchStrategy::setPhysicalBcCoef(RobinBcCoefStrategy<NDIM>* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(1, bc_coef));
    return;
} // setPhysicalBcCoef

void
RobinPhysBdryPatchStrategy::setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if !defined(NDEBUG)
    for (unsigned int l = 0; l < bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(bc_coefs[l]);
    }
#endif
    d_bc_coefs = bc_coefs;
    return;
} // setPhysicalBcCoefs

void
RobinPhysBdryPatchStrategy::setHomogeneousBc(bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
} // setHomogeneousBc

bool
RobinPhysBdryPatchStrategy::getHomogeneousBc() const
{
    return d_homogeneous_bc;
} // getHomogeneousBc

void
RobinPhysBdryPatchStrategy::preprocessRefine(Patch<NDIM>& /*fine*/,
                                             const Patch<NDIM>& /*coarse*/,
                                             const Box<NDIM>& /*fine_box*/,
                                             const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
RobinPhysBdryPatchStrategy::postprocessRefine(Patch<NDIM>& /*fine*/,
                                              const Patch<NDIM>& /*coarse*/,
                                              const Box<NDIM>& /*fine_box*/,
                                              const IntVector<NDIM>& /*ratio*/)
{
    // intentionally blank
    return;
} // postprocessRefine

void
RobinPhysBdryPatchStrategy::accumulateFromPhysicalBoundaryData(Patch<NDIM>& /*patch*/,
                                                               double /*fill_time*/,
                                                               const IntVector<NDIM>& /*ghost_width_to_fill*/)
{
    TBOX_ERROR("RobinPhysBdryPatchStrategy::accumulateFromPhysicalBoundaryData(): unimplemented\n");
    return;
} // accumulateFromPhysicalBoundaryData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
