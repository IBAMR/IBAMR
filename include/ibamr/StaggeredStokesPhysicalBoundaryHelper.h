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
     * \brief Add the viscous part of the TRACTION conditions to \a f_data_idx at the boundary faces where the normal
     * velocity is not prescribed. Pass the same \a u_data_idx and \a viscous_coef as to
     * IBTK::HierarchyMathOps::laplace(). The term is exact only for \a bdry_interp_type LINEAR. See
     * checkNormalTractionBdryInterpType().
     */
    void addNormalTractionViscousTerm(int f_data_idx,
                                      int u_data_idx,
                                      double viscous_coef,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      IBTK::BdryInterpType bdry_interp_type,
                                      int coarsest_ln = IBTK::invalid_level_number,
                                      int finest_ln = IBTK::invalid_level_number) const;

    /*!
     * \brief Stop with an error if \a bdry_interp_type is QUADRATIC, \a viscous_coef_nonzero is true, and the normal
     * velocity of a TRACTION condition in \a u_bc_coefs (or of a PSEUDO_TRACTION condition, if
     * \a include_pseudo_traction is true) is not prescribed at some boundary face of \a hierarchy. The terms of
     * addNormalTractionViscousTerm() are exact only for LINEAR. Call this function when setting up the solver, not for
     * each term. Pass \a viscous_coef_nonzero true for a variable viscous coefficient.
     */
    void checkNormalTractionBdryInterpType(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                           IBTK::BdryInterpType bdry_interp_type,
                                           bool viscous_coef_nonzero,
                                           double fill_time,
                                           SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> hierarchy,
                                           bool include_pseudo_traction = false) const;

    /*!
     * \brief Add the viscous part of the TRACTION and PSEUDO_TRACTION conditions to \a f_data_idx at the boundary faces
     * where the normal velocity is not prescribed. Pass the same \a u_data_idx, \a viscous_coef_data_idx, and
     * \a viscous_coef_interp_type as to IBTK::HierarchyMathOps::vc_laplace(). The term is exact only for
     * \a bdry_interp_type LINEAR. See checkNormalTractionBdryInterpType().
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
     * \brief Add the viscous part of the TRACTION conditions to \a f_data on one patch. See the version for a patch
     * hierarchy.
     */
    void addNormalTractionViscousTerm(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> f_data,
                                      SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> u_data,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                      double viscous_coef,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs) const;

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
