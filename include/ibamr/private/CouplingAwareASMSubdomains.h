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

#ifndef included_IBAMR_CouplingAwareASMSubdomains
#define included_IBAMR_CouplingAwareASMSubdomains

#include <ibamr/config.h>

#include <ibamr/ibamr_enums.h>

#include <tbox/Pointer.h>

#include <petscmat.h>

#include <PatchLevel.h>

#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace IBAMR
{
/*! \brief Cached level geometry for coupling-aware patch construction.
 *
 * Recreate after changing the level or DOF data. The supplied patch data must
 * outlive this object. Construction semantics are defined by
 * StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains()
 * and StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches().
 */
class CouplingAwareASMSubdomains
{
public:
    /*! \brief Build velocity adjacency and cell closures from live DOF data. */
    CouplingAwareASMSubdomains(int u_idx, int p_idx, SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level);

    /*! \brief Construct ordered overlap sets and their first-owner partition. */
    void constructSubdomains(std::vector<std::set<int>>& overlap,
                             std::vector<std::set<int>>& nonoverlap,
                             const std::vector<int>& num_dofs_per_proc,
                             Mat matrix,
                             int seed_axis,
                             int seed_stride,
                             CouplingAwareASMSeedTraversalOrder order,
                             CouplingAwareASMClosurePolicy policy,
                             double relative_zero_tol);

    /*! \brief Construct pressure patches as defined by
     * StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches().
     */
    void constructPressureCellPatches(std::vector<std::set<int>>& patches,
                                      std::vector<int>& pressure_seeds,
                                      const std::vector<int>& num_dofs_per_proc,
                                      Mat elasticity,
                                      int seed_stride,
                                      CouplingAwareASMSeedTraversalOrder order,
                                      CouplingAwareASMClosurePolicy policy,
                                      double relative_zero_tol);

private:
    /*! \brief Add lower-face component pairing when STRICT first needs it. */
    void buildSeedPairs();

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> d_level;
    int d_u_idx, d_p_idx;
    std::unordered_map<int, std::set<int>> d_adjacent_cells, d_cell_closures, d_seed_pairs;
    std::unordered_set<int> d_velocity_dofs;
    bool d_pairs_built = false;
};
} // namespace IBAMR
#endif
