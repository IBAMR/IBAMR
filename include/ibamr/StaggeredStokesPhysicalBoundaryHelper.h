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

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/StaggeredPhysicalBoundaryHelper.h>
#include <ibtk/ibtk_enums.h>

#include <tbox/Pointer.h>

#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
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
     * IBTK::HierarchyMathOps::laplace(). At a TRACTION boundary face where the normal velocity is not prescribed and
     * \a viscous_coef is nonzero, \a bdry_interp_type must be "LINEAR".
     */
    void addNormalTractionViscousTerm(int f_data_idx,
                                      int u_data_idx,
                                      double viscous_coef,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      const std::string& bdry_interp_type,
                                      int coarsest_ln = IBTK::invalid_level_number,
                                      int finest_ln = IBTK::invalid_level_number) const;

    /*!
     * \brief Add the viscous part of the TRACTION conditions to \a f_data on one patch. See the version for a patch
     * hierarchy.
     */
    void addNormalTractionViscousTerm(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> f_data,
                                      SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> u_data,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                      double viscous_coef,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      const std::string& bdry_interp_type) const;

    /*!
     * \brief Add the viscous part of the TRACTION and PSEUDO_TRACTION conditions to \a f_data_idx at the boundary faces
     * where the normal velocity is not prescribed. Pass the same \a u_data_idx, \a viscous_coef_data_idx, and
     * \a viscous_coef_interp_type as to IBTK::HierarchyMathOps::vc_laplace(). At a TRACTION or PSEUDO_TRACTION boundary
     * face where the normal velocity is not prescribed and the viscous coefficient is nonzero,
     * \a bdry_interp_type must be "LINEAR".
     */
    void addNormalTractionViscousTerm(int f_data_idx,
                                      int u_data_idx,
                                      int viscous_coef_data_idx,
                                      IBTK::VCInterpType viscous_coef_interp_type,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      const std::string& bdry_interp_type,
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
     * \brief Return the transaction component for a fill of the cell-centered pressure. Because the Robin data of
     * \a P_bc_coef are not a condition that the pressure satisfies, the coarse data are extended linearly beyond
     * physical boundaries.
     */
    static IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent makePressureTransactionComponent(
        int P_idx,
        const std::string& refine_type,
        bool use_cf_interpolation,
        const std::string& coarsen_type,
        const std::string& extrap_type,
        bool consistent_type_2_bdry,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef,
        SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern = nullptr,
        const std::string& bdry_interp_type = "LINEAR");

    /*!
     * \brief Return the transaction component that fills the ghost cells of the pressure data \a P_dst_idx from the
     * pressure data \a P_src_idx. The component is otherwise the same as for the version with one patch data index.
     */
    static IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent makePressureTransactionComponent(
        int P_dst_idx,
        int P_src_idx,
        const std::string& refine_type,
        bool use_cf_interpolation,
        const std::string& coarsen_type,
        const std::string& extrap_type,
        bool consistent_type_2_bdry,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef,
        SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM>> fill_pattern = nullptr,
        const std::string& bdry_interp_type = "LINEAR");

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
