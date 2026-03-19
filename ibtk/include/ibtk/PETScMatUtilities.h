// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2023 by the IBAMR developers
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

#include <PoissonSpecifications.h>

#include <algorithm>
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
                                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level);
    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to the
     * side-centered IB interpolation operator for anisotropic kernels.
     *
     * \warning This routine does not properly handle odd stencil sizes, nor
     * does it properly handle physical boundary conditions.
     */
    static void constructPatchLevelSCInterpOp(Mat& mat,
                                              void (*component_interp_fcn)(double r_lower, double* w),
                                              int component_interp_stencil,
                                              void (*transverse_interp_fcn)(double r_lower, double* w),
                                              int transverse_interp_stencil,
                                              Vec& X_vec,
                                              const std::vector<int>& num_dofs_per_proc,
                                              int dof_index_idx,
                                              SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level);

    /*!
     * \brief Standard one-dimensional B-spline interpolation functions.
     */
    static void bspline_2_interp_fcn(const double r, double* const w)
    {
        w[0] = 1.0 - r;
        w[1] = r;
        return;
    } // bspline_2_interp_fcn

    static const int bspline_2_interp_stencil = 2;

    static void bspline_3_interp_fcn(const double r, double* const w)
    {
        const double rr[4] = { r, r - 1.0, r - 2.0, r - 3.0 };
        for (int i = 0; i < 4; ++i)
        {
            const double x = std::abs(rr[i]);
            const double rp = x + 1.5;
            if (x <= 0.5)
            {
                w[i] = 0.5 * (-2.0 * rp * rp + 6.0 * rp - 3.0);
            }
            else if (x <= 1.5)
            {
                w[i] = 0.5 * (rp * rp - 6.0 * rp + 9.0);
            }
            else
            {
                w[i] = 0.0;
            }
        }
        return;
    } // bspline_3_interp_fcn

    static const int bspline_3_interp_stencil = 4;

    static void bspline_4_interp_fcn(const double r, double* const w)
    {
        const double rr[4] = { r, r - 1.0, r - 2.0, r - 3.0 };
        for (int i = 0; i < 4; ++i)
        {
            const double x = std::abs(rr[i]);
            const double rp = x + 2.0;
            if (x <= 1.0)
            {
                w[i] = (1.0 / 6.0) * (3.0 * std::pow(rp, 3) - 24.0 * std::pow(rp, 2) + 60.0 * rp - 44.0);
            }
            else if (x <= 2.0)
            {
                w[i] = (1.0 / 6.0) * (-std::pow(rp, 3) + 12.0 * std::pow(rp, 2) - 48.0 * rp + 64.0);
            }
            else
            {
                w[i] = 0.0;
            }
        }
        return;
    } // bspline_4_interp_fcn

    static const int bspline_4_interp_stencil = 4;

    static void bspline_5_interp_fcn(const double r, double* const w)
    {
        const double rr[6] = { r, r - 1.0, r - 2.0, r - 3.0, r - 4.0, r - 5.0 };
        for (int i = 0; i < 6; ++i)
        {
            const double x = std::abs(rr[i]);
            const double rp = x + 2.5;
            if (x <= 0.5)
            {
                w[i] = (1.0 / 24.0) *
                       (6.0 * std::pow(rp, 4) - 60.0 * std::pow(rp, 3) + 210.0 * std::pow(rp, 2) - 300.0 * rp + 155.0);
            }
            else if (x <= 1.5)
            {
                w[i] = (1.0 / 24.0) *
                       (-4.0 * std::pow(rp, 4) + 60.0 * std::pow(rp, 3) - 330.0 * std::pow(rp, 2) + 780.0 * rp - 655.0);
            }
            else if (x <= 2.5)
            {
                w[i] = (1.0 / 24.0) *
                       (std::pow(rp, 4) - 20.0 * std::pow(rp, 3) + 150.0 * std::pow(rp, 2) - 500.0 * rp + 625.0);
            }
            else
            {
                w[i] = 0.0;
            }
        }
        return;
    } // bspline_5_interp_fcn

    static const int bspline_5_interp_stencil = 6;

    static void bspline_6_interp_fcn(const double r, double* const w)
    {
        const double rr[6] = { r, r - 1.0, r - 2.0, r - 3.0, r - 4.0, r - 5.0 };
        for (int i = 0; i < 6; ++i)
        {
            const double x = std::abs(rr[i]);
            const double rp = x + 3.0;
            if (x <= 1.0)
            {
                w[i] = (1.0 / 60.0) * (2193.0 - 3465.0 * rp + 2130.0 * std::pow(rp, 2) - 630.0 * std::pow(rp, 3) +
                                       90.0 * std::pow(rp, 4) - 5.0 * std::pow(rp, 5));
            }
            else if (x <= 2.0)
            {
                w[i] = (1.0 / 120.0) * (-10974.0 + 12270.0 * rp - 5340.0 * std::pow(rp, 2) + 1140.0 * std::pow(rp, 3) -
                                        120.0 * std::pow(rp, 4) + 5.0 * std::pow(rp, 5));
            }
            else if (x <= 3.0)
            {
                w[i] = (1.0 / 120.0) * (7776.0 - 6480.0 * rp + 2160.0 * std::pow(rp, 2) - 360.0 * std::pow(rp, 3) +
                                        30.0 * std::pow(rp, 4) - std::pow(rp, 5));
            }
            else
            {
                w[i] = 0.0;
            }
        }
        return;
    } // bspline_6_interp_fcn

    static const int bspline_6_interp_stencil = 6;

    /*!
     * \brief Standard one-dimensional IB interpolation functions.
     */
    static void ib_3_interp_fcn(const double r, double* const w)
    {
        const double rr[4] = { r, r - 1.0, r - 2.0, r - 3.0 };
        for (int i = 0; i < 4; ++i)
        {
            const double x = std::abs(rr[i]);
            if (x < 0.5)
            {
                w[i] = (1.0 / 3.0) * (1.0 + std::sqrt(1.0 - 3.0 * x * x));
            }
            else if (x < 1.5)
            {
                const double t = 1.0 - x;
                w[i] = (1.0 / 6.0) * (5.0 - 3.0 * x - std::sqrt(1.0 - 3.0 * t * t));
            }
            else
            {
                w[i] = 0.0;
            }
        }
        return;
    } // ib_3_interp_fcn

    static const int ib_3_interp_stencil = 4;
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
        const double q = std::sqrt(-7.0 + 12.0 * r - 4.0 * r * r);
        w[0] = 0.125 * (5.0 - 2.0 * r - q);
        w[1] = 0.125 * (5.0 - 2.0 * r + q);
        w[2] = 0.125 * (-1.0 + 2.0 * r + q);
        w[3] = 0.125 * (-1.0 + 2.0 * r - q);
        return;
    } // ib_4_interp_fcn

    static const int ib_4_interp_stencil = 4;

    static void ib_5_interp_fcn(const double r, double* const w)
    {
        static const double K = (38.0 - std::sqrt(69.0)) / 60.0;
        const double rr[6] = { r, r - 1.0, r - 2.0, r - 3.0, r - 4.0, r - 5.0 };
        for (int i = 0; i < 6; ++i)
        {
            const double x = std::abs(rr[i]);
            if (x > 2.5)
            {
                w[i] = 0.0;
                continue;
            }
            const double x2 = x * x;
            const double x4 = x2 * x2;
            const double x6 = x4 * x2;
            const double radicand =
                std::max(0.0,
                         2.0 * (3123.0 - 6840.0 * K + 3600.0 * K * K - 12440.0 * x2 + 25680.0 * K * x2 -
                                12600.0 * K * K * x2 + 8080.0 * x4 - 8400.0 * K * x4 - 1400.0 * x6));
            const double phi = (136.0 - 40.0 * K - 40.0 * x2 + std::sqrt(radicand)) / 280.0;
            if (x <= 0.5)
            {
                w[i] = phi;
            }
            else if (x <= 1.5)
            {
                const double xr = x - 1.0;
                w[i] = (4.0 - 4.0 * phi - K - 4.0 * xr + 3.0 * K * xr - xr * xr + xr * xr * xr) / 6.0;
            }
            else
            {
                const double xr = x - 2.0;
                w[i] = (-2.0 + 2.0 * phi + 2.0 * K + xr - 3.0 * K * xr + 2.0 * xr * xr - xr * xr * xr) / 12.0;
            }
        }
        return;
    } // ib_5_interp_fcn

    static const int ib_5_interp_stencil = 6;

    static void ib_6_interp_fcn(const double r, double* const w)
    {
        const double rl = r - 2.0;
        const double r2 = rl * rl;
        const double r3 = r2 * rl;
        const double r4 = r3 * rl;
        const double r5 = r4 * rl;
        static const double K = (59.0 / 60.0) * (1.0 - std::sqrt(1.0 - (3220.0 / 3481.0)));
        static const double K2 = K * K;
        static const double alpha = 28.0;
        const double beta = (9.0 / 4.0) - (3.0 / 2.0) * (K + r2) + ((22.0 / 3.0) - 7.0 * K) * rl - (7.0 / 3.0) * r3;
        const double gamma = (1.0 / 4.0) * (((161.0 / 36.0) - (59.0 / 6.0) * K + 5.0 * K2) * (1.0 / 2.0) * r2 +
                                            (-(109.0 / 24.0) + 5.0 * K) * (1.0 / 3.0) * r4 + (5.0 / 18.0) * r5 * rl);
        const double discr = std::max(0.0, beta * beta - 4.0 * alpha * gamma);
        w[0] = (-beta + std::copysign(1.0, (3.0 / 2.0) - K) * std::sqrt(discr)) / (2.0 * alpha);
        w[1] = -3.0 * w[0] - (1.0 / 16.0) + (1.0 / 8.0) * (K + r2) + (1.0 / 12.0) * (3.0 * K - 1.0) * rl +
               (1.0 / 12.0) * r3;
        w[2] = 2.0 * w[0] + (1.0 / 4.0) + (1.0 / 6.0) * (4.0 - 3.0 * K) * rl - (1.0 / 6.0) * r3;
        w[3] = 2.0 * w[0] + (5.0 / 8.0) - (1.0 / 4.0) * (K + r2);
        w[4] = -3.0 * w[0] + (1.0 / 4.0) - (1.0 / 6.0) * (4.0 - 3.0 * K) * rl + (1.0 / 6.0) * r3;
        w[5] = w[0] - (1.0 / 16.0) + (1.0 / 8.0) * (K + r2) - (1.0 / 12.0) * (3.0 * K - 1.0) * rl - (1.0 / 12.0) * r3;
        return;
    } // ib_6_interp_fcn

    static const int ib_6_interp_stencil = 6;

    /*!
     * \brief Standard one-dimensional Piecewise linear interpolation function.
     */
    static void pwl_interp_fcn(const double r, double* const w)
    {
        bspline_2_interp_fcn(r, w);
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
    static void
    constructLinearProlongationOp_cell(Mat& mat,
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

/////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_PETScMatUtilities
