// Filename: PETScMatUtilities.h
// Created on 24 Aug 2010 by Boyce Griffith
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

#ifndef included_PETScMatUtilities
#define included_PETScMatUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <math.h>
#include <vector>

#include "PoissonSpecifications.h"
#include "petscao.h"
#include "petscmat.h"
#include "petscvec.h"
#include "tbox/Pointer.h"

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

namespace IBTK
{
/*!
 * \brief Class PETScMatUtilities provides utility functions for <A
 * HREF="http://www.mcs.anl.gov/petsc">PETSc</A> Mat objects.
 */
class PETScMatUtilities
{
public:
    /*!
     * \name Methods acting on SAMRAI::hier::PatchLevel and
     * SAMRAI::hier::Variable objects.
     */
    //\{

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to the
     * cell-centered Laplacian of a cell-centered variable restricted to a
     * single SAMRAI::hier::PatchLevel.
     */
    static void constructPatchLevelCCLaplaceOp(Mat& mat,
                                               const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                               SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
                                               double data_time,
                                               const std::vector<int>& num_dofs_per_proc,
                                               int dof_index_idx,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to the
     * cell-centered Laplacian of a cell-centered variable restricted to a
     * single SAMRAI::hier::PatchLevel.
     */
    static void constructPatchLevelCCLaplaceOp(Mat& mat,
                                               const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                               const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                               double data_time,
                                               const std::vector<int>& num_dofs_per_proc,
                                               int dof_index_idx,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to the
     * side-centered Laplacian of a side-centered variable restricted to a
     * single SAMRAI::hier::PatchLevel.
     */
    static void constructPatchLevelSCLaplaceOp(Mat& mat,
                                               const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                               const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                               double data_time,
                                               const std::vector<int>& num_dofs_per_proc,
                                               int dof_index_idx,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to the
     * side-centered IB interpolation operator for the provided kernel function.
     *
     * \warning This routine does not properly handle delta functions for which
     * interp_stencil is odd, nor does it properly handle physical boundary
     * conditions.
     */
    static void constructPatchLevelSCInterpOp(Mat& mat,
                                              void (*interp_fcn)(double r_lower, double* w),
                                              int interp_stencil,
                                              Vec& X_vec,
                                              const std::vector<int>& num_dofs_per_proc,
                                              int dof_index_idx,
                                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level);
    /*!
     * \brief Standard one-dimensional Peskin 4-pt delta function.
     *
     * \param r Normalized distance (by grid space h) between the IB point and
     * the lowermost stencil location.
     *
     * \param w Weights as a function of (normalized) distance between the IB
     * point and stencil locations. The first entry corresponds to the distance
     * \f$ r_0 = r \f$ between the IB point and lowermost stencil location. The
     * next normalized distance is taken as \f$ r_1 = r_0 + 1 \f$ and so on.
     */
    static void ib_4_interp_fcn(const double r, double* const w)
    {
        const double q = sqrt(-7.0 + 12.0 * r - 4.0 * r * r);
        w[0] = 0.125 * (5.0 - 2.0 * r - q);
        w[1] = 0.125 * (5.0 - 2.0 * r + q);
        w[2] = 0.125 * (-1.0 + 2.0 * r + q);
        w[3] = 0.125 * (-1.0 + 2.0 * r - q);
        return;
    } // ib_4_interp_fcn

    static const int ib_4_interp_stencil = 4;

    /*!
     * \brief Standard one-dimensional Piecewise linear interpolation function.
     */
    static void pwl_interp_fcn(const double r, double* const w)
    {
        w[0] = 1.0 - r;
        w[1] = r;
        return;
    } // pwl_interp_fcn

    static const int pwl_interp_stencil = 2;

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to data
     * prolongation from a coarser level to a finer level.
     */
    static void constructProlongationOp(Mat& mat,
                                        const std::string& op_type,
                                        int dof_index_idx,
                                        const std::vector<int>& num_fine_dofs_per_proc,
                                        const std::vector<int>& num_coarse_dofs_per_proc,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > fine_patch_level,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarse_patch_level,
                                        const AO& coarse_level_ao,
                                        const int coarse_ao_offset = 0);

    /*!
     * \brief Construct a diagonal matrix that scales the rows of adjoint
     * (transpose) of prolongation matrix to get a suitable restriction matrix.
     *
     * \NOTE We store the diagonal enteries into a Vec rather than a Mat.
     */
    static void constructRestrictionScalingOp(Mat& P, Vec& L);

    /*!
     * \brief Partition the patch level into subdomains suitable to be used for
     * additive Schwarz method.
     */
    static void
    constructPatchLevelASMSubdomains(std::vector<IS>& is_overlap,
                                     std::vector<IS>& is_nonoverlap,
                                     const SAMRAI::hier::IntVector<NDIM>& box_size,
                                     const SAMRAI::hier::IntVector<NDIM>& overlap_size,
                                     const std::vector<int>& num_dofs_per_proc,
                                     int dof_index_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM> > cf_boundary);

    //\}

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScMatUtilities();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScMatUtilities(const PETScMatUtilities& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScMatUtilities& operator=(const PETScMatUtilities& that);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to cc-data
     * and conservative prolongation from a coarser level to a finer level.
     */
    static void
    constructConservativeProlongationOp_cell(Mat& mat,
                                             int dof_index_idx,
                                             const std::vector<int>& num_fine_dofs_per_proc,
                                             const std::vector<int>& num_coarse_dofs_per_proc,
                                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > fine_patch_level,
                                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarse_patch_level,
                                             const AO& coarse_level_ao,
                                             const int coarse_ao_offset);
    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to sc-data
     * and RT0 prolongation from a coarser level to a finer level.
     */
    static void
    constructRT0ProlongationOp_side(Mat& mat,
                                    int dof_index_idx,
                                    const std::vector<int>& num_fine_dofs_per_proc,
                                    const std::vector<int>& num_coarse_dofs_per_proc,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > fine_patch_level,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarse_patch_level,
                                    const AO& coarse_level_ao,
                                    const int coarse_ao_offset);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to sc-data
     * and linear prolongation from a coarser level to a finer level.
     */
    static void
    constructLinearProlongationOp_side(Mat& mat,
                                       int dof_index_idx,
                                       const std::vector<int>& num_fine_dofs_per_proc,
                                       const std::vector<int>& num_coarse_dofs_per_proc,
                                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > fine_patch_level,
                                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarse_patch_level,
                                       const AO& coarse_level_ao,
                                       const int coarse_ao_offset);

    /*!
     * \brief Partition the patch level into subdomains suitable to be used for
     * additive Schwarz method for a cc-variable.
     */
    static void
    constructPatchLevelASMSubdomains_cell(std::vector<IS>& is_overlap,
                                          std::vector<IS>& is_nonoverlap,
                                          const SAMRAI::hier::IntVector<NDIM>& box_size,
                                          const SAMRAI::hier::IntVector<NDIM>& overlap_size,
                                          const std::vector<int>& num_dofs_per_proc,
                                          int dof_index_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM> > cf_boundary);
    /*!
     * \brief Partition the patch level into subdomains suitable to be used for
     * additive Schwarz method for a sc-variable.
     */
    static void
    constructPatchLevelASMSubdomains_side(std::vector<IS>& is_overlap,
                                          std::vector<IS>& is_nonoverlap,
                                          const SAMRAI::hier::IntVector<NDIM>& box_size,
                                          const SAMRAI::hier::IntVector<NDIM>& overlap_size,
                                          const std::vector<int>& num_dofs_per_proc,
                                          int dof_index_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM> > cf_boundary);
};
} // namespace IBTK

/////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PETScMatUtilities
