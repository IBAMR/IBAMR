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

#include <tbox/Pointer.h>

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
     * \brief Boundary tags.
     */
    static const short int NORMAL_TRACTION_BDRY, NORMAL_VELOCITY_BDRY, ALL_BDRY;

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
     * \brief Set normal velocity ghost cell values to enforce discrete
     * divergence-free condition in the ghost cells abutting the physical boundary.
     *
     * \note The default behavior is to set these values only in cells adjacent to
     * boundary locations where normal traction conditions are imposed. Values can also
     * be set where normal velocity boundary conditions, or both.
     */
    void enforceDivergenceFreeConditionAtBoundary(int u_data_idx,
                                                  int coarsest_ln = IBTK::invalid_level_number,
                                                  int finest_ln = IBTK::invalid_level_number,
                                                  short int bdry_tag = NORMAL_TRACTION_BDRY) const;

    /*!
     * \brief Set normal velocity ghost cell values to enforce discrete
     * divergence-free condition in the ghost cells abutting the physical boundary.
     *
     * \note The default behavior is to set these values only in cells adjacent to
     * boundary locations where normal traction conditions are imposed. Values can also
     * be set where normal velocity boundary conditions, or both.
     */
    void enforceDivergenceFreeConditionAtBoundary(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> u_data,
                                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                  short int bdry_tag = NORMAL_TRACTION_BDRY) const;

    /*!
     * \brief Add the part of the viscous term that imposes TRACTION conditions at the boundary faces where the normal
     * velocity is not prescribed.
     *
     * At each such boundary face for which the entry of \a u_bc_coefs for the normal component is a
     * StokesBcCoefStrategy with TRACTION conditions, this function adds \f$ D (u_I - u_{div}) / h^2 \f$ to the normal
     * component of \a f_data_idx at the boundary face.  Here \f$ D \f$ is \a viscous_coef, \f$ u_I \f$ is the normal
     * velocity on the next face inside the domain, \f$ u_{div} \f$ is the normal velocity on the next face outside the
     * domain that makes the discrete divergence of \a u_data_idx vanish in the ghost cell, and \f$ h \f$ is the grid
     * spacing normal to the boundary.  Nothing is added at other boundary faces, and \a u_data_idx is not modified.
     *
     * The data \a f_data_idx must already hold the viscous term \f$ D \Delta u \f$ evaluated with the velocity ghost
     * values set by the velocity boundary conditions, which give the normal velocity ghost value that reflects the
     * normal velocity on the next face inside the domain.  The tangential velocity ghost values of \a u_data_idx must
     * have been set.  Together with the pressure boundary value \f$ p = -g \f$, the result imposes
     * \f$ -p + 2 \mu \partial u_n / \partial x_n = g \f$.
     *
     * The added term is exact only if the pressure ghost value at the boundary is the linear extrapolation
     * \f$ p_G = 2 p_b - p_I \f$, which is the case for the boundary interpolation type "LINEAR" and not for
     * "QUADRATIC".  Set \a linear_pressure_extrapolation to indicate whether the pressure ghost values are filled that
     * way.  It is an error for it to be false if the term is nonzero at some boundary face, that is, if \a viscous_coef
     * is nonzero and there is a boundary face with TRACTION conditions at which the normal velocity is not prescribed.
     */
    void addNormalTractionViscousTerm(int f_data_idx,
                                      int u_data_idx,
                                      double viscous_coef,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      bool linear_pressure_extrapolation,
                                      int coarsest_ln = IBTK::invalid_level_number,
                                      int finest_ln = IBTK::invalid_level_number) const;

    /*!
     * \brief Add the part of the viscous term that imposes TRACTION conditions on one patch; see the version of this
     * function for a patch hierarchy.
     */
    void addNormalTractionViscousTerm(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> f_data,
                                      SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> u_data,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                      double viscous_coef,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      bool linear_pressure_extrapolation) const;

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
