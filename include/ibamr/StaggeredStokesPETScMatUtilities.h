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
     *
     * The matrix is for homogeneous boundary conditions.  The momentum
     * equation for the normal velocity on a physical boundary at which the
     * normal velocity is not prescribed contains the pressure in the ghost cell
     * outside the domain.  If \a p_bc_coef is not null, that pressure is the
     * linear extrapolation defined by the Robin coefficients of \a p_bc_coef,
     * p_G = -((a*h - 2*b)/(a*h + 2*b))*p_I + 2*h*g/(a*h + 2*b), where p_I is
     * the pressure in the adjacent interior cell and h is the grid spacing
     * normal to the boundary, and the matrix includes the term with p_I.  If
     * \a p_bc_coef is null, the matrix does not include the ghost pressure.
     */
    static void constructPatchLevelMACStokesOp(Mat& mat,
                                               const SAMRAI::solv::PoissonSpecifications& u_problem_coefs,
                                               const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                               double data_time,
                                               const std::vector<int>& num_dofs_per_proc,
                                               int u_dof_index_idx,
                                               int p_dof_index_idx,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level,
                                               SAMRAI::solv::RobinBcCoefStrategy<NDIM>* p_bc_coef = nullptr);

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
