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

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBAMR_private_StaggeredStokesPETScMatUtilities_inl_h
#define included_IBAMR_private_StaggeredStokesPETScMatUtilities_inl_h

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/StaggeredStokesPETScMatUtilities.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Narrow access surface for implementation files and tests that need to
 * inspect coupling-aware ASM internals.
 */
class StaggeredStokesPETScMatUtilitiesPrivateAccess
{
public:
    static const std::unordered_map<int, std::vector<int>>&
    getVelocityDofToAdjacentCellDofs(const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data)
    {
        return map_data.velocity_dof_to_adjacent_cell_dofs;
    }

    static const std::unordered_map<int, std::vector<int>>&
    getCellDofToClosureDofs(const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data)
    {
        return map_data.cell_dof_to_closure_dofs;
    }

    static const std::unordered_map<int, int>&
    getVelocityDofToComponentAxis(const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data)
    {
        return map_data.velocity_dof_to_component_axis;
    }

    static const std::unordered_map<int, std::vector<int>>& getVelocityDofToPairedSeedVelocityDofs(
        const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data)
    {
        return map_data.velocity_dof_to_paired_seed_velocity_dofs;
    }

    static SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>>
    getSourcePatchLevel(const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data)
    {
        return map_data.source_patch_level;
    }

    static int getSourceUDofIndex(const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data)
    {
        return map_data.source_u_dof_index_idx;
    }

    static int getSourcePDofIndex(const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data)
    {
        return map_data.source_p_dof_index_idx;
    }

    static bool velocityMapDataIsBuilt(const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data)
    {
        return map_data.velocity_maps_are_built;
    }

    static bool cellClosureMapIsBuilt(const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data)
    {
        return map_data.cell_closure_map_is_built;
    }

    static bool
    velocitySeedPairMapIsBuilt(const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data)
    {
        return map_data.velocity_seed_pair_map_is_built;
    }

    static bool isEmpty(const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data)
    {
        return !map_data.source_patch_level && map_data.source_u_dof_index_idx == IBTK::invalid_index &&
               map_data.source_p_dof_index_idx == IBTK::invalid_index && !map_data.velocity_maps_are_built &&
               !map_data.cell_closure_map_is_built && !map_data.velocity_seed_pair_map_is_built &&
               map_data.velocity_dof_to_adjacent_cell_dofs.empty() && map_data.cell_dof_to_closure_dofs.empty() &&
               map_data.velocity_dof_to_component_axis.empty() &&
               map_data.velocity_dof_to_paired_seed_velocity_dofs.empty();
    }

    static void computePatchLevelCouplingAwareASMSeedVelocityDofs(
        std::vector<int>& seed_velocity_dofs,
        int u_dof_index_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level,
        const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData& map_data,
        int seed_velocity_axis = 0,
        int seed_velocity_stride = 1,
#if (NDIM == 2)
        CouplingAwareASMSeedTraversalOrder seed_traversal_order = CouplingAwareASMSeedTraversalOrder::I_J)
#else
        CouplingAwareASMSeedTraversalOrder seed_traversal_order = CouplingAwareASMSeedTraversalOrder::I_J_K)
#endif
    {
        StaggeredStokesPETScMatUtilities::computePatchLevelCouplingAwareASMSeedVelocityDofs(seed_velocity_dofs,
                                                                                            u_dof_index_idx,
                                                                                            patch_level,
                                                                                            map_data,
                                                                                            seed_velocity_axis,
                                                                                            seed_velocity_stride,
                                                                                            seed_traversal_order);
    }

    static void
    findCoupledCellDofsFromA00(std::set<int>& involved_cell_dofs,
                               Mat A00_mat,
                               const std::set<int>& seed_velocity_dofs,
                               const std::unordered_map<int, std::vector<int>>& velocity_dof_to_adjacent_cell_dofs,
                               const std::unordered_map<int, std::vector<int>>& cell_dof_to_closure_dofs,
                               double relative_numerical_zero_tol = IBTK_RELATIVE_NUMERICAL_ZERO_TOL)
    {
        StaggeredStokesPETScMatUtilities::findCoupledCellDofsFromA00(involved_cell_dofs,
                                                                     A00_mat,
                                                                     seed_velocity_dofs,
                                                                     velocity_dof_to_adjacent_cell_dofs,
                                                                     cell_dof_to_closure_dofs,
                                                                     relative_numerical_zero_tol);
    }

    static void constructA00VelocitySubmatrix(Mat& A00_velocity_mat,
                                              Mat A00_mat,
                                              const std::vector<int>& num_dofs_per_proc,
                                              int u_dof_index_idx,
                                              int p_dof_index_idx,
                                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level)
    {
        StaggeredStokesPETScMatUtilities::constructA00VelocitySubmatrix(
            A00_velocity_mat, A00_mat, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, patch_level);
    }
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_private_StaggeredStokesPETScMatUtilities_inl_h
