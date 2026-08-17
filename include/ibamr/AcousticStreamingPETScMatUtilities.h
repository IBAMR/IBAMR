// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_AcousticStreamingPETScMatUtilities
#define included_IBAMR_AcousticStreamingPETScMatUtilities

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibtk/ibtk_enums.h"

#include "IntVector.h"
#include "tbox/Pointer.h"

#include "petscmat.h"

#include <array>
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
 * \brief Class AcousticStreamingPETScMatUtilities provides utility functions for
 * <A HREF="http://www.mcs.anl.gov/petsc">PETSc</A> Mat objects used in solving acoustic
 * streaming equations.
 *
 * see IBAMR::AcousticStreamingHierarchyIntegrator and IBAMR::FOAcousticStreamingPETScLevelSolver
 */
class AcousticStreamingPETScMatUtilities
{
public:
    /*!
     * \name Methods acting on SAMRAI::hier::PatchLevel and
     * SAMRAI::hier::Variable objects.
     */
    //\{

    /*!
     * \brief Construct a parallel PETSc Mat object corresponding to a MAC
     * discretization of the first order acoustic streaming equations on a
     * single SAMRAI::hier::PatchLevel.
     */
    static void constructPatchLevelFOAcousticStreamingOp(
        Mat& mat,
        double omega,
        double sound_speed,
        int rho_idx,
        int mu_idx,
        int lambda_idx,
        int chi_idx,
        const std::array<std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>, 2>& u_bc_coefs,
        double data_time,
        const std::vector<int>& num_dofs_per_proc,
        int u_dof_index_idx,
        int p_dof_index_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
        IBTK::VCInterpType mu_interp_type = IBTK::VC_HARMONIC_INTERP);

    /*!
     * \brief First-order acoustic/elastic/rigid matrix operator.
     *
     * The assembled real/imaginary block system is
     *
     *   [ +w rho I + L_e      L_a + P_r        0          G_a ] [v_r]
     *   [ L_a + P_r          -w rho I - L_e    G_a        0   ] [v_i]
     *   [ 0                   chi_a D_rho       w/c^2 I   0   ] [p_r] = 0
     *   [ chi_a D_rho         0                 0         -w/c^2 I] [p_i] = 0
     *
     * where
     *
     *   L_a(v) = -div[chi_a {mu (grad v + grad v^T) + lambda div(v) I}],
     *
     *   L_e(v) = (1/w) div[chi_e {Gamma (grad v + grad v^T) + zeta div(v) I}],
     *
     *   G_a(p) = grad(chi_a p),
     *
     *   P_r = chi_r/kappa.
     *
     *  Coefficient naming convention:
     *
     *   mu      : fluid shear viscosity
     *   lambda  : fluid bulk/dilatational viscosity
     *   Gamma   : elastic shear modulus
     *   zeta    : elastic Lame coefficient
     *
     *
     *  IMPORTANT:
     *
     * The constitutive coefficient fields are assumed to be ALREADY
     * multiplied by their respective indicators:
     *
     *   acoustic_mu_idx      = chi_a * mu
     *   acoustic_lambda_idx  = chi_a * lambda
     *   elastic_gamma_idx    = chi_e * Gamma
     *   elastic_zeta_idx     = chi_e * zeta
     *
     * rigid_penalty_idx stores chi_r/kappa.
     *
     * acoustic_indicator_idx stores chi_a = 1 - chi_e and is used for:
     *
     *   1. grad(chi_a p),
     *   2. chi_a D_rho(v) in the pressure/mass rows.
     */
    static void constructPatchLevelFOAcousticStreamingOp(
        Mat& mat,
        double omega,
        double sound_speed,
        int rho_idx,
        int acoustic_mu_idx,
        int acoustic_lambda_idx,
        int elastic_gamma_idx,
        int elastic_zeta_idx,
        int rigid_penalty_idx,
        int acoustic_indicator_idx,
        const std::array<std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>, 2>& u_bc_coefs,
        double data_time,
        const std::vector<int>& num_dofs_per_proc,
        int u_dof_index_idx,
        int p_dof_index_idx,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > patch_level,
        IBTK::VCInterpType coeff_interp_type = IBTK::VC_HARMONIC_INTERP);

    //\}

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AcousticStreamingPETScMatUtilities() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AcousticStreamingPETScMatUtilities(const AcousticStreamingPETScMatUtilities& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AcousticStreamingPETScMatUtilities& operator=(const AcousticStreamingPETScMatUtilities& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_AcousticStreamingPETScMatUtilities
