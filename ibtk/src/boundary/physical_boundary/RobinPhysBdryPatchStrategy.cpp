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

#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/samrai_compatibility_names.h"

#include "SAMRAIBox.h"
#include "SAMRAIComponentSelector.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIRobinBcCoefStrategy.h"

#include <set>
#include <vector>

#include "ibtk/app_namespaces.h" // IWYU pragma: keep

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
RobinPhysBdryPatchStrategy::setPatchDataIndices(const SAMRAIComponentSelector& patch_data_indices)
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
RobinPhysBdryPatchStrategy::setPhysicalBcCoef(SAMRAIRobinBcCoefStrategy* const bc_coef)
{
    setPhysicalBcCoefs(std::vector<SAMRAIRobinBcCoefStrategy*>(1, bc_coef));
    return;
} // setPhysicalBcCoef

void
RobinPhysBdryPatchStrategy::setPhysicalBcCoefs(const std::vector<SAMRAIRobinBcCoefStrategy*>& bc_coefs)
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
RobinPhysBdryPatchStrategy::preprocessRefine(SAMRAIPatch& /*fine*/,
                                             const SAMRAIPatch& /*coarse*/,
                                             const SAMRAIBox& /*fine_box*/,
                                             const SAMRAIIntVector& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
RobinPhysBdryPatchStrategy::postprocessRefine(SAMRAIPatch& /*fine*/,
                                              const SAMRAIPatch& /*coarse*/,
                                              const SAMRAIBox& /*fine_box*/,
                                              const SAMRAIIntVector& /*ratio*/)
{
    // intentionally blank
    return;
} // postprocessRefine

void
RobinPhysBdryPatchStrategy::accumulateFromPhysicalBoundaryData(SAMRAIPatch& /*patch*/,
                                                               double /*fill_time*/,
                                                               const SAMRAIIntVector& /*ghost_width_to_fill*/)
{
    TBOX_ERROR("RobinPhysBdryPatchStrategy::accumulateFromPhysicalBoundaryData(): unimplemented\n");
    return;
} // accumulateFromPhysicalBoundaryData

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
