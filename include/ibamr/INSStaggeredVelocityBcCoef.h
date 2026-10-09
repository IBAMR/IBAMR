// ---------------------------------------------------------------------
//
// Copyright (c) 2008 - 2026 by the IBAMR developers
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

#ifndef included_IBAMR_INSStaggeredVelocityBcCoef
#define included_IBAMR_INSStaggeredVelocityBcCoef

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include <ibamr/StokesBcCoefStrategy.h>
#include <ibamr/ibamr_enums.h>

#include <tbox/Pointer.h>

#include <BoxArray.h>
#include <IntVector.h>
#include <SideIndex.h>

#include <array>
#include <map>
#include <utility>
#include <vector>

namespace IBAMR
{
class INSStaggeredHierarchyIntegrator;
class StokesSpecifications;
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoundaryBox;
template <int DIM>
class Patch;
template <int DIM>
class Variable;
} // namespace hier
namespace pdat
{
template <int DIM, class TYPE>
class ArrayData;
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
 * \brief Class INSStaggeredVelocityBcCoef is a concrete StokesBcCoefStrategy
 * that is used to specify velocity boundary conditions for the staggered grid
 * incompressible Navier-Stokes solver.
 *
 * This class interprets pure Dirichlet boundary conditions on the velocity as
 * prescribed velocity boundary conditions, whereas pure Neumann boundary
 * conditions are interpreted as prescribed traction (stress) boundary
 * conditions.  These are translated into Dirichlet and generalized Neumann
 * boundary conditions, respectively, for the velocity.
 *
 * Dirichlet, true traction, and pseudo-traction boundary conditions are
 * all supported.
 *
 * With TRACTION boundary conditions, the value set for a tangential velocity
 * component includes the tangential derivative of the normal velocity along the
 * boundary. A corner is where the boundary meets an adjacent physical boundary
 * (an edge in three dimensions). The difference across a corner needs a normal
 * velocity value beyond it. If the adjacent boundary prescribes the normal
 * velocity, that value is 2*u_b - u, where u is the normal velocity on the
 * nearest face of the boundary and u_b is the prescribed value; otherwise it is
 * extrapolated linearly along the boundary. The difference is first-order
 * accurate at a corner, and so is the outermost tangential ghost column of the
 * patch data, where the difference is taken one cell further inside the data.
 */
class INSStaggeredVelocityBcCoef : public StokesBcCoefStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    INSStaggeredVelocityBcCoef(unsigned int comp_idx,
                               const INSStaggeredHierarchyIntegrator* fluid_solver,
                               const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                               TractionBcType traction_bc_type,
                               bool homogeneous_bc = false);

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredVelocityBcCoef() = default;

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \param bc_coefs  IBTK::Vector of boundary condition specification objects
     */
    void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    void setSolutionTime(double solution_time);

    /*!
     * \brief Set the current time interval.
     */
    void setTimeInterval(double current_time, double new_time);

    /*!
     * \name StokesBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Set the StokesSpecifications object used by this boundary condition
     * specification object.
     *
     * \param problem_coefs   Problem coefficients
     */
    void setStokesSpecifications(const StokesSpecifications* problem_coefs) override;

    /*!
     * \brief Set the target velocity data index to use when setting physical
     * boundary conditions and the time at which it is defined.
     */
    void setTargetVelocityPatchDataIndex(int u_target_data_idx) override;

    /*!
     * \brief Clear the target velocity data index used when setting physical
     * boundary conditions.
     */
    void clearTargetVelocityPatchDataIndex() override;

    /*!
     * \brief Set the target pressure data index to use when setting physical
     * boundary conditions and the time at which it is defined.
     */
    void setTargetPressurePatchDataIndex(int p_target_data_idx) override;

    /*!
     * \brief Clear the target pressure data index used when setting physical
     * boundary conditions.
     */
    void clearTargetPressurePatchDataIndex() override;

    //\}

    /*!
     * \name Extended SAMRAI::solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Set the target data index.
     */
    void setTargetPatchDataIndex(int target_idx) override;

    /*!
     * \brief Clear the target data index.
     */
    void clearTargetPatchDataIndex() override;

    /*!
     * \brief Set whether the class is filling homogeneous or inhomogeneous
     * boundary conditions.
     */
    void setHomogeneousBc(bool homogeneous_bc) override;

    /*!
     * \brief Accumulate the transpose of the dependence of \f$ \gamma \f$ on the normal velocity.
     *
     * With TRACTION boundary conditions, \f$ \gamma \f$ for a tangential component includes the tangential
     * derivative of the normal velocity along the boundary. This function adds the transpose of that term to the
     * normal velocity on the faces it depends on, in the patch interior or in its ghost cells.
     *
     * \see IBTK::ExtendedRobinBcCoefStrategy::accumulateFromBcCoefs()
     */
    void accumulateFromBcCoefs(const SAMRAI::pdat::ArrayData<NDIM, double>& gcoef_data,
                               const SAMRAI::hier::Patch<NDIM>& patch,
                               const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                               double fill_time) const override;

    //\}

    /*!
     * \brief Return the dependence of \f$ \gamma \f$ on the normal velocity at one location of a boundary.
     *
     * With TRACTION boundary conditions, \f$ \gamma \f$ for a tangential component includes the tangential
     * derivative of the normal velocity along the boundary. This function returns the normal velocity values that
     * \f$ \gamma \f$ depends on at the location \a i, each with the derivative of \f$ \gamma \f$ with respect to
     * it. It is the dependence whose transpose accumulateFromBcCoefs() applies. The result is empty if this object
     * does not impose TRACTION conditions or if its component is normal to the boundary.
     *
     * \param i  Location in the index space of the Robin coefficients of this component, at which the traction is
     *        prescribed.
     * \param patch  Patch with the geometry that setBcCoefs() requires.
     * \param bdry_box  Boundary box that contains the location.
     * \param ghost_box  Ghost box of the data that the returned side indices refer to. A location outside this box
     *        is extrapolated from the two nearest locations inside it.
     * \param fill_time  Time at which the boundary conditions of the normal component are evaluated.
     */
    std::vector<std::pair<SAMRAI::pdat::SideIndex<NDIM>, double>>
    getNormalVelocityStencil(const SAMRAI::hier::Index<NDIM>& i,
                             const SAMRAI::hier::Patch<NDIM>& patch,
                             const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                             const SAMRAI::hier::Box<NDIM>& ghost_box,
                             double fill_time) const;

    /*!
     * \name Implementation of SAMRAI::solv::RobinBcCoefStrategy interface.
     */
    //\{

    /*!
     * \brief Function to fill arrays of Robin boundary condition coefficients
     * at a patch boundary.
     *
     * \note In the original SAMRAI::solv::RobinBcCoefStrategy interface, it was
     * assumed that \f$ b = (1-a) \f$.  In the new interface, \f$a\f$ and
     * \f$b\f$ are independent.
     *
     * \see SAMRAI::solv::RobinBcCoefStrategy::setBcCoefs()
     *
     * \param acoef_data  Boundary coefficient data.
     *        The array will have been defined to include index range for
     *        corresponding to the boundary box \a bdry_box and appropriate for
     *        the alignment of the given variable.  If this is a null pointer,
     *        then the calling function is not interested in a, and you can
     *        disregard it.
     * \param bcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it is to be
     *        filled with the b coefficient.
     * \param gcoef_data  Boundary coefficient data.
     *        This array is exactly like \a acoef_data, except that it is to be
     *        filled with the g coefficient.
     * \param variable    Variable to set the coefficients for.
     *        If implemented for multiple variables, this parameter can be used
     *        to determine which variable's coefficients are being sought.
     * \param patch       Patch requiring bc coefficients.
     * \param bdry_box    Boundary box showing where on the boundary the coefficient data is
     *needed.
     * \param fill_time   Solution time corresponding to filling, for use when coefficients are
     *time-dependent.
     */
    void setBcCoefs(SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>& acoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>& bcoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double>>& gcoef_data,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>>& variable,
                    const SAMRAI::hier::Patch<NDIM>& patch,
                    const SAMRAI::hier::BoundaryBox<NDIM>& bdry_box,
                    double fill_time = 0.0) const override;

    /*
     * \brief Return how many cells past the edge or corner of the patch the
     * object can fill.
     *
     * The "extension" used here is the number of cells that a boundary box
     * extends past the patch in the direction parallel to the boundary.
     *
     * Note that the inability to fill the sufficient number of cells past the
     * edge or corner of the patch may preclude the child class from being used
     * in data refinement operations that require the extra data, such as linear
     * refinement.
     *
     * The boundary box that setBcCoefs() is required to fill should not extend
     * past the limits returned by this function.
     */
    SAMRAI::hier::IntVector<NDIM> numberOfExtensionsFillable() const override;

    //\}

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredVelocityBcCoef() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredVelocityBcCoef(const INSStaggeredVelocityBcCoef& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredVelocityBcCoef& operator=(const INSStaggeredVelocityBcCoef& that) = delete;

    /*
     * Component of the velocity which this boundary condition specification is
     * to operate on.
     */
    const unsigned int d_comp_idx;

    /*
     * The fluid solver.
     */
    const INSStaggeredHierarchyIntegrator* d_fluid_solver;

    /*
     * The boundary condition specification objects for the velocity.
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;

    /*
     * Return the physical domain of the level of patch, including its periodic
     * images, computed on first use for each refinement ratio.
     */
    const SAMRAI::hier::BoxArray<NDIM>& getPhysicalDomain(const SAMRAI::hier::Patch<NDIM>& patch) const;

    /*
     * The physical domain, with its periodic images, refined by each refinement
     * ratio for which it has been requested.
     */
    mutable std::map<std::array<int, NDIM>, SAMRAI::hier::BoxArray<NDIM>> d_physical_domain;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSStaggeredVelocityBcCoef
