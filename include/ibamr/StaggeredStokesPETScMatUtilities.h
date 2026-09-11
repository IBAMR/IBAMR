// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#ifndef included_IBAMR_StaggeredStokesPETScMatUtilities
#define included_IBAMR_StaggeredStokesPETScMatUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/ibamr_enums.h>

#include <tbox/Pointer.h>

#include <petscao.h>
#include <petscmat.h>

#include <IntVector.h>
#include <PoissonSpecifications.h>

#include <set>
#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
template <int DIM>
class CoarseFineBoundary;
} // namespace hier
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StaggeredStokesPETScMatUtilities provides utility functions for
 * <A HREF="http://www.mcs.anl.gov/petsc">PETSc</A> Mat objects.
 */
class StaggeredStokesPETScMatUtilities
{
public:
    /*!
     * \name Methods acting on SAMRAI::hier::PatchLevel and
     * SAMRAI::hier::Variable objects.
     */
    //\{

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to a MAC
     * discretization of the time-dependent incompressible Stokes equations on a
     * single SAMRAI::hier::PatchLevel.
     */
    static void constructPatchLevelMACStokesOp(Mat& mat,
                                               const SAMRAI::solv::PoissonSpecifications& u_problem_coefs,
                                               const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                               double data_time,
                                               const std::vector<int>& num_dofs_per_proc,
                                               int u_dof_index_idx,
                                               int p_dof_index_idx,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level);

    /*!
     * \brief Partition the patch level into subdomains suitable to be used for
     * additive Schwarz method.
     */
    static void
    constructPatchLevelASMSubdomains(std::vector<std::set<int>>& is_overlap,
                                     std::vector<std::set<int>>& is_nonoverlap,
                                     const SAMRAI::hier::IntVector<NDIM>& box_size,
                                     const SAMRAI::hier::IntVector<NDIM>& overlap_size,
                                     const std::vector<int>& num_dofs_per_proc,
                                     int u_dof_index_idx,
                                     int p_dof_index_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM>> cf_boundary);

    /*! \brief Construct velocity-seeded coupling-aware ASM subdomains.
     *
     * The borrowed matrix uses the full coupled DOF numbering of the supplied
     * level and supports local row access. Only stored velocity columns enter
     * the numerical row graph: an entry is retained when its magnitude exceeds
     * max(n * epsilon, relative_zero_tol) times the velocity-row maximum, where
     * n includes stored velocity zeros. Pressure entries do not affect this test.
     * RELAXED closes all incident cells and retains expanded velocities. STRICT
     * adds lower-face seed components of the same cell and retains only cells
     * whose full velocity stencil is present after one row-neighbor expansion.
     *
     * Seeds of seed_axis are ordered by logical coordinates, de-duplicated and
     * then sampled by seed_stride. seed_axis must be in [0, NDIM), seed_stride
     * must be positive, and relative_zero_tol must be finite and nonnegative.
     * The traversal order must match NDIM. The outer vector preserves that order. Each
     * inner set uses global coupled IDs; nonoverlap assigns each locally owned
     * DOF to its first containing overlap. Incomplete local coverage is an error.
     * DOF data must have valid adjacent-cell ghosts and remain unchanged during
     * construction. This routine constructs subdomains; it does not apply them.
     */
    static void construct_patch_level_coupling_aware_asm_subdomains(
        std::vector<std::set<int>>& overlap,
        std::vector<std::set<int>>& nonoverlap,
        const std::vector<int>& num_dofs_per_proc,
        int u_dof_index_idx,
        int p_dof_index_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level,
        Mat matrix,
        int seed_axis = 0,
        int seed_stride = 1,
#if (NDIM == 2)
        CouplingAwareASMSeedTraversalOrder order = CouplingAwareASMSeedTraversalOrder::I_J,
#else
        CouplingAwareASMSeedTraversalOrder order = CouplingAwareASMSeedTraversalOrder::I_J_K,
#endif
        CouplingAwareASMClosurePolicy policy = CouplingAwareASMClosurePolicy::RELAXED,
        double relative_zero_tol = 1.0e-14);

    /*! \brief Construct ordered pressure-cell-seeded CAV patches.
     *
     * The borrowed elasticity matrix uses full global velocity-pressure IDs,
     * supports row access, and has numerically zero pressure rows and columns.
     * An entry is retained when its magnitude exceeds
     * max(n * epsilon, relative_zero_tol) times its row maximum, where n counts
     * all stored entries, including zeros. Each seed's incident MAC velocities
     * expand once through the retained row-or-column graph. Without expansion,
     * both policies return the standard Vanka patch. Otherwise RELAXED closes
     * all incident cells and retains expanded velocities; STRICT retains only
     * cells whose complete velocity stencil is supported.
     *
     * Pressure seeds are sorted in logical traversal order and de-duplicated
     * before applying the positive seed_stride. patches[k] corresponds to
     * pressure_seeds[k]; each set contains increasing unique global IDs.
     * The traversal order must match NDIM and relative_zero_tol must be finite
     * and nonnegative. DOF data need valid adjacent-cell ghosts and must remain
     * unchanged during construction. Only one MPI rank is supported.
     * The matrix is read anew on each call and is not retained. This operation
     * constructs patches, without partitioning ownership or applying corrections.
     */
    static void construct_patch_level_pressure_cell_seeded_cav_patches(
        std::vector<std::set<int>>& patches,
        std::vector<int>& pressure_seeds,
        const std::vector<int>& num_dofs_per_proc,
        int u_dof_index_idx,
        int p_dof_index_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level,
        Mat elasticity,
        int seed_stride = 1,
#if (NDIM == 2)
        CouplingAwareASMSeedTraversalOrder order = CouplingAwareASMSeedTraversalOrder::I_J,
#else
        CouplingAwareASMSeedTraversalOrder order = CouplingAwareASMSeedTraversalOrder::I_J_K,
#endif
        CouplingAwareASMClosurePolicy policy = CouplingAwareASMClosurePolicy::RELAXED,
        double relative_zero_tol = 1.0e-14);

    /*!
     * \brief Partition the patch level into subdomains suitable to be used for
     * PCFieldSplit preconditioner.
     */
    static void constructPatchLevelFields(std::vector<std::set<int>>& is_field,
                                          std::vector<std::string>& is_field_name,
                                          const std::vector<int>& num_dofs_per_proc,
                                          int u_dof_index_idx,
                                          int p_dof_index_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to data
     * prolongation from a coarser level to a finer level.
     */
    static void constructProlongationOp(Mat& mat,
                                        const std::string& u_op_type,
                                        const std::string& p_op_type,
                                        int u_dof_index_idx,
                                        int p_dof_index_idx,
                                        const std::vector<int>& num_fine_dofs_per_proc,
                                        const std::vector<int>& num_coarse_dofs_per_proc,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> fine_patch_level,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> coarse_patch_level,
                                        const AO& coarse_level_ao,
                                        const int u_coarse_ao_offset,
                                        const int p_coarse_ao_offset);

    //\}

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    StaggeredStokesPETScMatUtilities() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesPETScMatUtilities(const StaggeredStokesPETScMatUtilities& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesPETScMatUtilities& operator=(const StaggeredStokesPETScMatUtilities& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesPETScMatUtilities
