// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
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

#include <ibtk/ibtk_utilities.h>
#include <ibtk/private/ArrayDataTransaction.h>
#include <ibtk/private/CachedGhostCellFill.h>
#include <ibtk/private/data_kind.h>

#include <tbox/Pointer.h>

#include <IntVector.h>
#include <PatchDescriptor.h>
#include <PatchHierarchy.h>
#include <VariableDatabase.h>

#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
namespace detail
{
/////////////////////////////// LOCAL ////////////////////////////////////////

namespace
{
// Whether two components differ in a member other than the patch data indices and the Robin coefficient objects, which
// resetTransactionComponents() accepts when the properties of the indices are the same. The fill patterns compare by
// shape, so a null pattern equals a default one, and the operator for the physical boundary values by pointer.
bool
differs_beyond_reset(const HierarchyGhostCellInterpolation::InterpolationTransactionComponent& a,
                     const HierarchyGhostCellInterpolation::InterpolationTransactionComponent& b)
{
    // Only the pattern matters, so use a ghost width that does not limit the stencil width of a NoCorners pattern.
    const IntVector<NDIM> wide(std::numeric_limits<int>::max());
    return a.d_refine_op_name != b.d_refine_op_name || a.d_use_cf_bdry_interpolation != b.d_use_cf_bdry_interpolation ||
           a.d_coarsen_op_name != b.d_coarsen_op_name || a.d_phys_bdry_extrap_type != b.d_phys_bdry_extrap_type ||
           a.d_consistent_type_2_bdry != b.d_consistent_type_2_bdry ||
           fill_pattern_shape(a.d_fill_pattern, wide) != fill_pattern_shape(b.d_fill_pattern, wide) ||
           a.d_phys_bdry_type != b.d_phys_bdry_type || a.d_sc_phys_bdry_op != b.d_sc_phys_bdry_op ||
           a.d_extend_coarse_data_linearly != b.d_extend_coarse_data_linearly;
} // differs_beyond_reset
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

HierarchyGhostCellInterpolation&
CachedGhostCellFill::get(const std::vector<HierarchyGhostCellInterpolation::InterpolationTransactionComponent>& comps,
                         Pointer<PatchHierarchy<NDIM>> hierarchy,
                         int coarsest_ln,
                         int finest_ln,
                         const bool homogeneous_bc)
{
    coarsest_ln = coarsest_ln == invalid_level_number ? 0 : coarsest_ln;
    finest_ln = finest_ln == invalid_level_number ? hierarchy->getFinestLevelNumber() : finest_ln;

    // resetTransactionComponents() can follow a change of the indices only if these properties stay equal.
    Pointer<PatchDescriptor<NDIM>> descriptor = VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor();
    std::vector<IndexProperties> index_properties(comps.size());
    for (std::size_t k = 0; k < comps.size(); ++k)
    {
        const auto& comp = comps[k];
        IndexProperties& properties = index_properties[k];
        properties.dst_kind = get_data_kind(comp.d_dst_data_idx);
        properties.src_kind = get_data_kind(comp.d_src_data_idx);
        properties.dst_ghosts = descriptor->getPatchDataFactory(comp.d_dst_data_idx)->getGhostCellWidth();
        properties.src_ghosts = descriptor->getPatchDataFactory(comp.d_src_data_idx)->getGhostCellWidth();
        properties.in_place = comp.d_dst_data_idx == comp.d_src_data_idx;
        for (const auto& coef : comp.d_robin_bc_coefs)
        {
            properties.null_robin_bc_coefs.push_back(!coef);
        }
    }

    bool rebuild = d_levels.changed(*hierarchy, coarsest_ln, finest_ln) || !d_fill || comps.size() != d_comps.size();
    bool reset = false;
    for (std::size_t k = 0; !rebuild && k < comps.size(); ++k)
    {
        const auto& a = comps[k];
        const auto& b = d_comps[k];
        if (differs_beyond_reset(a, b) || !(index_properties[k] == d_index_properties[k]))
        {
            rebuild = true;
        }
        else if (a.d_dst_data_idx != b.d_dst_data_idx || a.d_src_data_idx != b.d_src_data_idx ||
                 a.d_robin_bc_coefs != b.d_robin_bc_coefs)
        {
            reset = true;
        }
    }
    if (rebuild)
    {
        d_fill = new HierarchyGhostCellInterpolation();
        d_fill->initializeOperatorState(comps, hierarchy, coarsest_ln, finest_ln);
    }
    else if (reset)
    {
        d_fill->resetTransactionComponents(comps);
    }
    d_fill->setHomogeneousBc(homogeneous_bc);
    d_comps = comps;
    d_index_properties = std::move(index_properties);
    return *d_fill;
} // get

void
CachedGhostCellFill::clear()
{
    d_fill.setNull();
    d_comps.clear();
    d_index_properties.clear();
    d_levels = HierarchyLevelsIdentity();
    return;
} // clear

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
