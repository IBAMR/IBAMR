// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2020 by the IBAMR developers
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

#ifndef included_IBAMR_INSVCStaggeredPressureBcCoef
#define included_IBAMR_INSVCStaggeredPressureBcCoef

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/StokesBcCoefStrategy.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/ibtk_enums.h"

#include "IntVector.h"
#include "tbox/Pointer.h"

#include <vector>

namespace IBAMR
{
class INSVCStaggeredHierarchyIntegrator;
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
 * \brief Class INSVCStaggeredPressureBcCoef is a concrete StokesBcCoefStrategy
 * that is used to specify pressure boundary conditions for the staggered grid
 * incompressible Navier-Stokes solver with variable coefficients.
 *
 * This class interprets pure Dirichlet boundary conditions on the velocity as
 * prescribed velocity boundary conditions, whereas pure Neumann boundary
 * conditions are interpreted as prescribed traction (stress) boundary
 * conditions.  These are translated into Neumann and generalized Dirichlet
 * boundary conditions, respectively, for the pressure.
 */
class INSVCStaggeredPressureBcCoef : public StokesBcCoefStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    INSVCStaggeredPressureBcCoef(const INSVCStaggeredHierarchyIntegrator* fluid_solver,
                                 const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
                                 TractionBcType traction_bc_type,
                                 bool homogeneous_bc = false);

    /*!
     * \brief Destructor.
     */
    ~INSVCStaggeredPressureBcCoef() = default;

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

    //\}

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
    void setBcCoefs(SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& acoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& bcoef_data,
                    SAMRAI::tbox::Pointer<SAMRAI::pdat::ArrayData<NDIM, double> >& gcoef_data,
                    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& variable,
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

    /*!
     * \brief Set the interpolation type to bring cell centered viscosity
     * to side centers
     */
    void setViscosityInterpolationType(IBTK::VCInterpType mu_interp_type);

protected:
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSVCStaggeredPressureBcCoef() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSVCStaggeredPressureBcCoef(const INSVCStaggeredPressureBcCoef& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSVCStaggeredPressureBcCoef& operator=(const INSVCStaggeredPressureBcCoef& that) = delete;

    /*
     * The fluid solver.
     */
    const INSVCStaggeredHierarchyIntegrator* d_fluid_solver;

    /*
     * The boundary condition specification objects for the velocity.
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;

    /*
     * The type of interpolation to bring cell centered viscosity to side centers.
     */
    IBTK::VCInterpType d_mu_interp_type;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_INSVCStaggeredPressureBcCoef
