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
     * \brief At open boundaries with TRACTION conditions, set the normal velocity ghost value that imposes the normal
     * traction condition through the viscous term.
     *
     * The value is set on the ghost face adjacent to the boundary, at each boundary face at which the normal velocity
     * is not prescribed and for which the entry of \a u_bc_coefs for the normal component is a StokesBcCoefStrategy
     * with TRACTION conditions.  It is the value for which the second-order viscous term at the boundary face,
     * together with the pressure boundary value \f$ p = -g \f$, imposes \f$ -p + 2 \mu \partial u_n / \partial x_n = g
     * \f$.  It is computed from the normal velocity on the boundary face and the first face inside the domain and the
     * tangential velocity ghost values, so this function must be called after those ghost values have been set.  No
     * data are modified at other boundary faces.
     */
    void setNormalTractionGhostValues(int u_data_idx,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                                      int coarsest_ln = IBTK::invalid_level_number,
                                      int finest_ln = IBTK::invalid_level_number) const;

    /*!
     * \brief Set the normal velocity ghost values that impose TRACTION conditions on one patch; see the version of
     * this function for a patch hierarchy.
     */
    void setNormalTractionGhostValues(SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double>> u_data,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                      const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs) const;

    /*!
     * \brief Return the normal velocity on the ghost face of the ghost cell \a i_g that makes the discrete
     * divergence of \a u_data vanish in that cell.
     *
     * The ghost face is the face of \a i_g normal to \a normal_axis that is farther from the physical boundary: the
     * lower face when \a is_lower is true and the upper face otherwise. The value of \a u_data on the ghost face is not
     * read; the other faces of \a i_g are.
     *
     * \param u_data      Velocity data.
     * \param i_g         Index of the ghost cell.
     * \param normal_axis Axis normal to the physical boundary.
     * \param is_lower    Whether the physical boundary is the lower boundary along \a normal_axis.
     * \param dx          Grid spacings.
     */
    static double divergenceFreeNormalGhostValue(const SAMRAI::pdat::SideData<NDIM, double>& u_data,
                                                 const SAMRAI::hier::Index<NDIM>& i_g,
                                                 unsigned int normal_axis,
                                                 bool is_lower,
                                                 const double* dx);

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
