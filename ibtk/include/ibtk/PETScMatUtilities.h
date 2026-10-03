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

#include <ibtk/IBKernelConcepts.h>
#include <ibtk/ibtk_enums.h>

#include <tbox/Pointer.h>

#include <petscao.h>
#include <petscmat.h>
#include <petscvec.h>

#include <Box.h>
#include <Index.h>
#include <PoissonSpecifications.h>

#include <array>
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
namespace pdat
{
template <int DIM, class TYPE>
class SideData;
} // namespace pdat
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
     * Each local point's cell must lie in a locally owned patch or its
     * one-cell neighborhood.
     *
     * Axis selects the side-normal coordinate. With data centering included in
     * the grid coordinate q, the first stencil index is floor(q+1/2)-(N-1)/2
     * for odd width N. For even N it is floor(q)-N/2+1 in the normal direction
     * and ceil(q)-N/2 in other directions. The evaluator receives r = q minus
     * the first stencil index, following the coordinate convention of
     * IBKernelEvaluatorScalar. The evaluator and X_vec are borrowed for this
     * call.
     *
     * An existing mat is destroyed and replaced; the caller owns the new matrix.
     *
     * \warning Physical boundary conditions are not handled: a stencil point outside the domain at a
     * non-periodic boundary has no DOF, so its weight is dropped from the row instead of being folded or
     * renormalized, and a row whose stencil crosses such a boundary sums to less than one. A one-time warning
     * is logged when this happens.
     */
    template <IBKernelEvaluatorCartesian<double, PetscScalar> Evaluator>
    static void constructPatchLevelSCInterpOp(Mat& mat,
                                              const Evaluator& evaluator,
                                              Vec X_vec,
                                              const std::vector<int>& num_dofs_per_proc,
                                              int dof_index_idx,
                                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level);

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to data
     * prolongation from a coarser level to a finer level.
     *
     * The data centering is that of the variable dof_index_idx. Cell-centered
     * data support the op_type "CONSERVATIVE" and "LINEAR", and side-centered
     * data support "RT0" and "LINEAR". For cell-centered data, "LINEAR" is the
     * matrix of IBTK::CartCellDoubleLinearRefine: it interpolates linearly in
     * each coordinate between the two nearest coarse cell centers, extends the
     * nearest coarse value as a constant between a physical boundary and the
     * first coarse cell center, and wraps at periodic boundaries. The restriction
     * built from this matrix with constructRestrictionScalingOp() agrees with
     * IBTK::CartCellDoubleLinearCoarsen on levels without a coarse-fine interface;
     * that class documents the difference at an interface.
     *
     * The rows are the fine-level DOFs and the columns are the coarse-level DOFs,
     * each numbered as in dof_index_idx; coarse_level_ao and coarse_ao_offset map
     * the coarse level's application ordering to the PETSc ordering. The AO is
     * borrowed and must remain valid during the call. The physical domain of the
     * coarse level must consist of a single box, which is checked only in debug
     * builds. Any existing mat is destroyed and replaced.
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
        //! Read-only array borrowed until the matching VecRestoreArrayRead.
        const double* d_positions = nullptr;
        //! Grid spacings and physical domain origin.
        std::array<double, NDIM> d_dx, d_x_lower;
        //! Lower index of the physical domain.
        SAMRAI::hier::Index<NDIM> d_domain_lower;
        //! Number of local IB points and first local matrix row.
        int d_n_local_points = 0, d_row_lower = 0;
        //! Component stencil boxes for each IB point.
        std::vector<std::array<SAMRAI::hier::Box<NDIM>, NDIM>> d_stencil_boxes;
        //! Shared DOF index data for each local IB point's patch, used to
        //! read global column indices without repeating the patch lookup.
        std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, int>>> d_dof_index_data;
    };

    /*! \brief Assemble matrix rows for one velocity component. */
    template <int Axis, IBKernelEvaluatorCartesian<double, PetscScalar> Evaluator>
    static void constructSCInterpOpAxis(SCInterpOpData& data, const Evaluator& evaluator);

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
     * \brief Construct the cell-centered "LINEAR" case of constructProlongationOp().
     *
     * The boundary behavior, numbering, AO lifetime, and replacement of mat are
     * those documented there.
     */
    static void
    constructLinearProlongationOp_cell(Mat& mat,
                                       int dof_index_idx,
                                       const std::vector<int>& num_fine_dofs_per_proc,
                                       const std::vector<int>& num_coarse_dofs_per_proc,
                                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> fine_patch_level,
                                       SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> coarse_patch_level,
                                       AO coarse_level_ao,
                                       int coarse_ao_offset);

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
