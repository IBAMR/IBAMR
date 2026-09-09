// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2026 by the IBAMR developers
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

#ifndef included_IBTK_PETScMatUtilities
#define included_IBTK_PETScMatUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/ibtk_enums.h>

#include <tbox/Pointer.h>

#include <petscao.h>
#include <petscmat.h>
#include <petscvec.h>

#include <Box.h>
#include <Index.h>
#include <PoissonSpecifications.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class PatchLevel;
template <int DIM>
class CoarseFineBoundary;
template <int DIM>
class IntVector;
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
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level);

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
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level);

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
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to the
     * side-centered viscous operator of a side-centered velocity variable
     * restricted to a single SAMRAI::hier::PatchLevel.
     *
     * \note The scaling factors of \f$ C \f$ and \f$ D \f$ variables in
     * the PoissonSpecification object are passed separately and are denoted
     * by \f$ \beta \f$ and \f$ \alpha \f$, respectively.
     */
    static void constructPatchLevelVCSCViscousOp(Mat& mat,
                                                 const SAMRAI::solv::PoissonSpecifications& poisson_spec,
                                                 double alpha,
                                                 double beta,
                                                 const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                                 double data_time,
                                                 const std::vector<int>& num_dofs_per_proc,
                                                 int dof_index_idx,
                                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level,
                                                 VCInterpType mu_interp_type = VC_HARMONIC_INTERP);

    /*!
     * \brief Construct a matrix mapping side-centered velocity on patch_level to IB points.
     *
     * X_vec contains consecutive NDIM coordinates for each IB point and
     * determines the matrix row ordering. The side-centered data at
     * dof_index_idx contain global column indices; num_dofs_per_proc gives
     * the column counts on each rank. The physical domain must be a single
     * box, and local index data must cover the stencils of local IB points.
     * Insufficient DOF ghost storage is a fatal error.
     *
     * For every velocity component 0 <= Axis < NDIM, Evaluator must provide:
     *
     * - static constexpr get_stencil_widths<Axis>(), returning
     *   std::array<int, NDIM> with strictly positive entries;
     * - evaluate<Axis>(const std::array<double, NDIM>&), callable on a const
     *   evaluator and returning exactly std::array<double, N>, where N is
     *   the product of those widths.
     *
     * Each r[d] is the displacement from the first stencil point to the IB
     * point, divided by the grid spacing. Result entries correspond to the
     * stencil points with coordinate zero varying fastest; grid-spacing
     * factors are not applied. The evaluator must supply weights consistent
     * with this ordering and the stencil placement below. These mathematical
     * requirements are the caller's responsibility, not compile-time checks.
     * IBKernelTensorProductEvaluator implements this interface.
     *
     * Odd widths use the nearest grid point, choosing the higher index at a
     * tie. Even widths bracket the point using the component's grid centering.
     * No kernel registration is required. The evaluator is borrowed for this call.
     * An existing mat is destroyed and replaced; the caller owns the new matrix.
     *
     * \warning Physical boundary conditions are not handled.
     */
    template <class Evaluator>
    static void constructPatchLevelSCInterpOp(Mat& mat,
                                              const Evaluator& evaluator,
                                              Vec& X_vec,
                                              const std::vector<int>& num_dofs_per_proc,
                                              int dof_index_idx,
                                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to data
     * prolongation from a coarser level to a finer level.
     */
    static void constructProlongationOp(Mat& mat,
                                        const std::string& op_type,
                                        int dof_index_idx,
                                        const std::vector<int>& num_fine_dofs_per_proc,
                                        const std::vector<int>& num_coarse_dofs_per_proc,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> fine_patch_level,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> coarse_patch_level,
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
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM>> cf_boundary);

    //\}

protected:
private:
    /*! \brief Interpolation stencil geometry and borrowed IB positions. */
    struct SCInterpOpData
    {
        /*! \brief Allocate the matrix and determine stencil boxes and local patches. */
        SCInterpOpData(Mat& mat,
                       Vec X,
                       const std::array<std::array<int, NDIM>, NDIM>& stencil_widths,
                       const std::vector<int>& num_dofs_per_proc,
                       int dof_index_idx,
                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level);
        /*! \brief Restore the borrowed position array. */
        ~SCInterpOpData();
        /*! \brief Disallow copying borrowed array access. */
        SCInterpOpData(const SCInterpOpData&) = delete;
        /*! \brief Disallow assigning borrowed array access. */
        SCInterpOpData& operator=(const SCInterpOpData&) = delete;
        /*! \brief Finish matrix assembly. */
        void assemble();

        //! Caller-owned matrix handle.
        Mat& d_mat;
        //! Borrowed vector; must remain alive through restoration of d_positions.
        Vec d_X;
        //! Array borrowed from VecGetArray until the matching VecRestoreArray.
        double* d_positions = nullptr;
        //! Grid spacings and physical domain origin.
        std::array<double, NDIM> d_dx, d_x_lower;
        //! Lower index of the physical domain.
        SAMRAI::hier::Index<NDIM> d_domain_lower;
        //! Number of local IB points and first local matrix row.
        int d_n_local_points = 0, d_row_lower = 0;
        //! Local patches and component stencil boxes for each IB point.
        std::vector<int> d_patch_numbers;
        std::vector<std::vector<SAMRAI::hier::Box<NDIM>>> d_stencil_boxes;
        //! Borrowed hierarchy data used to read global column indices.
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> d_level;
        int d_dof_index_idx;
    };

    /*! \brief Assemble matrix rows for one velocity component. */
    template <int Axis, class Evaluator>
    static void construct_sc_interp_op_axis(SCInterpOpData& data, const Evaluator& evaluator);

    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScMatUtilities() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScMatUtilities(const PETScMatUtilities& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScMatUtilities& operator=(const PETScMatUtilities& that) = delete;

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to cc-data
     * and conservative prolongation from a coarser level to a finer level.
     */
    static void
    constructConservativeProlongationOp_cell(Mat& mat,
                                             int dof_index_idx,
                                             const std::vector<int>& num_fine_dofs_per_proc,
                                             const std::vector<int>& num_coarse_dofs_per_proc,
                                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> fine_patch_level,
                                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> coarse_patch_level,
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
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> fine_patch_level,
                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> coarse_patch_level,
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
                                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> fine_patch_level,
                                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> coarse_patch_level,
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
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM>> cf_boundary);
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
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::CoarseFineBoundary<NDIM>> cf_boundary);
};
} // namespace IBTK

#include <ibtk/private/PETScMatUtilities-inl.h>

/////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_PETScMatUtilities
