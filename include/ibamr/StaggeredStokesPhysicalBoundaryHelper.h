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
     * \brief Add the part of the variable-coefficient viscous term that imposes TRACTION and PSEUDO_TRACTION conditions
     * at the boundary faces where the normal velocity is not prescribed.
     *
     * The viscous term is \f$ \nabla \cdot D (\nabla u + \nabla u^T) \f$, as computed by
     * IBTK::HierarchyMathOps::vc_laplace(), in which the coefficient \f$ D \f$ is the node-centered (two dimensions) or
     * edge-centered (three dimensions) data \a viscous_coef_data_idx and \a viscous_coef_interp_type selects the
     * average that gives its cell-centered values.  At each boundary face at which the normal velocity is not
     * prescribed and for which the entry of \a u_bc_coefs for the normal component is a StokesBcCoefStrategy, this
     * function adds
     * \f[ 2 (D_I (u_I - u_B) + D_G (u_B - u_G)) / h^2 \f]
     * to the normal component of \a f_data_idx at the boundary face, and for PSEUDO_TRACTION conditions it also
     * subtracts \f$ D_B (u_I - u_{div}) / h^2 \f$.  Here \f$ u_B \f$, \f$ u_I \f$, and \f$ u_G \f$ are the normal
     * velocities on the boundary face and on the next faces inside and outside the domain, \f$ u_{div} \f$ is the
     * normal velocity on the next face outside the domain that makes the discrete divergence of \a u_data_idx vanish
     * in the ghost cell, \f$ D_I \f$ and \f$ D_G \f$ are the cell-centered coefficients in the cell abutting the
     * boundary and in the ghost cell, \f$ D_B \f$ is the average of the coefficient over the boundary face, and
     * \f$ h \f$ is the grid spacing normal to the boundary.  Nothing is added at other boundary faces, and
     * \a u_data_idx is not modified.
     *
     * The data \a f_data_idx must already hold the viscous term evaluated with the velocity ghost values set by the
     * velocity boundary conditions and with the same ghost values of \a viscous_coef_data_idx.  Together with the
     * pressure boundary value \f$ p = -g \f$, the result imposes
     * \f$ -p + 2 \mu \partial u_n / \partial x_n = g \f$ for TRACTION conditions and
     * \f$ -p + \mu \partial u_n / \partial x_n = g \f$ for PSEUDO_TRACTION conditions.
     *
     * The added term is exact only if the pressure ghost value at the boundary is the linear extrapolation
     * \f$ p_G = 2 p_b - p_I \f$.  Set \a linear_pressure_extrapolation to indicate whether the pressure ghost values
     * are filled that way.  It is an error for it to be false if there is a boundary face to which this function
     * applies.
     */
    void addNormalTractionViscousTerm(int f_data_idx,
                                      int u_data_idx,
                                      int viscous_coef_data_idx,
                                      IBTK::VCInterpType viscous_coef_interp_type,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      bool linear_pressure_extrapolation,
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
