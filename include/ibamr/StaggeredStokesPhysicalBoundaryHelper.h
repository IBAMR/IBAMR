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

#ifndef included_IBAMR_StaggeredStokesPhysicalBoundaryHelper
#define included_IBAMR_StaggeredStokesPhysicalBoundaryHelper

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibtk/StaggeredPhysicalBoundaryHelper.h>
#include <ibtk/ibtk_enums.h>

#include <tbox/Pointer.h>

#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class SideData;
} // namespace pdat
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
 * \brief Class StaggeredStokesPhysicalBoundaryHelper provides helper functions
 * to enforce physical boundary conditions for a staggered grid discretization
 * of the incompressible (Navier-)Stokes equations.
 */
class StaggeredStokesPhysicalBoundaryHelper : public IBTK::StaggeredPhysicalBoundaryHelper
{
public:
    /*!
     * \brief Default constructor.
     */
    StaggeredStokesPhysicalBoundaryHelper() = default;

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesPhysicalBoundaryHelper() = default;

    /*!
     * \brief At Dirichlet boundaries, set values to enforce normal velocity
     * boundary conditions at the boundary.
     */
    void
    enforceNormalVelocityBoundaryConditions(int u_data_idx,
                                            int p_data_idx,
                                            const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                            double fill_time,
                                            bool homogeneous_bc,
                                            int coarsest_ln = IBTK::invalid_level_number,
                                            int finest_ln = IBTK::invalid_level_number) const;

    /*!
     * \brief Return the net outward flux through the physical boundary of the
     * normal velocity that \p u_bc_coefs prescribe at \p fill_time.
     *
     * \param u_data_idx Side-centered velocity data that the boundary condition
     * objects take as their target.
     * \param wgt_cc_idx Cell-centered composite-grid weights, which are zero on
     * covered cells.
     */
    double computeNetBoundaryFlux(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                  int u_data_idx,
                                  double fill_time,
                                  int wgt_cc_idx) const;

    /*!
     * \brief Add the viscous part of the TRACTION and PSEUDO_TRACTION conditions to \a f_data_idx at the boundary faces
     * where the normal velocity is not prescribed. Pass the same \a u_data_idx and \a viscous_coef as to
     * IBTK::HierarchyMathOps::laplace(). Pass the interpolation type of the pressure ghost values as
     * \a bdry_interp_type.
     *
     * At a boundary face with velocity \f$ u_B \f$, the function adds \f$ D t / h^2 \f$, in which \f$ D \f$ is
     * \a viscous_coef, \f$ h \f$ is the grid spacing normal to the boundary, \f$ u_I \f$ and \f$ u_{I2} \f$ are
     * the velocities at the next two faces inside the domain, and \f$ u_g^{\mathrm{div}} \f$ is the ghost value
     * that makes the divergence of the ghost cell zero:
     * - LINEAR: \f$ t = u_I - u_g^{\mathrm{div}} \f$, for TRACTION only;
     * - QUADRATIC: \f$ t = \frac{11}{3} u_I - 2 u_B - \frac{2}{3} u_{I2} - u_g^{\mathrm{div}} \f$, minus
     *   \f$ \frac{4}{3} (u_I - u_g^{\mathrm{div}}) \f$ for PSEUDO_TRACTION.
     *
     * For LINEAR, the PSEUDO_TRACTION term is the TRACTION term less \f$ D (u_I - u_g^{\mathrm{div}}) / h^2 \f$, which
     * is zero.
     */
    void addNormalTractionViscousTerm(int f_data_idx,
                                      int u_data_idx,
                                      double viscous_coef,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      IBTK::BdryInterpType bdry_interp_type,
                                      int coarsest_ln = IBTK::invalid_level_number,
                                      int finest_ln = IBTK::invalid_level_number) const;

    /*!
     * \brief Add the viscous part of the TRACTION and PSEUDO_TRACTION conditions to \a f_data on one patch. See the
     * version for a patch hierarchy.
     */
    void addNormalTractionViscousTerm(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> f_data,
                                      SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> u_data,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                      double viscous_coef,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      IBTK::BdryInterpType bdry_interp_type) const;

    /*!
     * \brief Add the viscous part of the TRACTION and PSEUDO_TRACTION conditions to \a f_data_idx at the boundary faces
     * where the normal velocity is not prescribed. Pass the same \a u_data_idx, \a viscous_coef_data_idx, and
     * \a viscous_coef_interp_type as to IBTK::HierarchyMathOps::vc_laplace(). Pass the interpolation type of the
     * pressure ghost values as \a bdry_interp_type.
     *
     * With the notation of the constant-coefficient version, in which \f$ u_G \f$ is the ghost value in
     * \a u_data_idx, and with \f$ D_I \f$, \f$ D_{I2} \f$, and \f$ D_G \f$ the viscous coefficients of the next
     * two cells inside the domain and of the ghost cell, the function adds \f$ t / h^2 \f$ with
     * - LINEAR: \f$ t = 2 [D_I (u_I - u_B) + D_G (u_B - u_G)] \f$;
     * - QUADRATIC: the LINEAR term plus \f$ 2 D_I (u_I - u_B) - \frac{2}{3} D_{I2} (u_{I2} - u_I) \f$.
     *
     * For PSEUDO_TRACTION, \f$ t \f$ also includes \f$ -c D_B (u_I - u_g^{\mathrm{div}}) \f$, in which
     * \f$ c = 1 \f$ for LINEAR and \f$ c = 4/3 \f$ for QUADRATIC, and \f$ D_B \f$ is the viscous coefficient at the
     * boundary face. This term differs from that of the constant-coefficient version because this operator is the
     * stress form, whose boundary row already includes the cross-derivative term of the viscous stress.
     */
    void addNormalTractionViscousTerm(int f_data_idx,
                                      int u_data_idx,
                                      int viscous_coef_data_idx,
                                      IBTK::VCInterpType viscous_coef_interp_type,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      IBTK::BdryInterpType bdry_interp_type,
                                      int coarsest_ln = IBTK::invalid_level_number,
                                      int finest_ln = IBTK::invalid_level_number) const;

    /*!
     * \brief Set \a r_data_idx to \a f_data_idx minus \a u_data_idx at the boundary faces at which the normal velocity
     * is prescribed, over the specified levels.
     */
    void setResidualAtDirichletBoundaries(int r_data_idx,
                                          int f_data_idx,
                                          int u_data_idx,
                                          int coarsest_ln = IBTK::invalid_level_number,
                                          int finest_ln = IBTK::invalid_level_number) const;

    /*!
     * \brief Setup physical boundary condition specification objects for
     * simultaneously filling velocity and pressure data.
     */
    static void setupBcCoefObjects(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                   SAMRAI::solv::RobinBcCoefStrategy<NDIM>* p_bc_coef,
                                   int u_target_data_idx,
                                   int p_target_data_idx,
                                   bool homogeneous_bc);

    /*!
     * \brief Reset physical boundary condition specification objects.
     */
    static void resetBcCoefObjects(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                   SAMRAI::solv::RobinBcCoefStrategy<NDIM>* p_bc_coef);

protected:
private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    StaggeredStokesPhysicalBoundaryHelper(const StaggeredStokesPhysicalBoundaryHelper& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    StaggeredStokesPhysicalBoundaryHelper& operator=(const StaggeredStokesPhysicalBoundaryHelper& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_StaggeredStokesPhysicalBoundaryHelper
