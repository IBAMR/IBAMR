// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_INSVCStaggeredConservativeMassMomentumRKIntegrator
#define included_IBAMR_INSVCStaggeredConservativeMassMomentumRKIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibamr/STSMassFluxIntegrator.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace hier
{
template <int DIM>
class BasePatchHierarchy;
} // namespace hier
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSVCStaggeredConservativeMassMomentumRKIntegrator integrates
 * the staggered density field
 *
 *  \f$ \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho u) = S(x,t) \f$
 *
 * and computes the conservative form of the convective operator
 * \f$ \nabla \cdot (\rho u u)\f$.
 *
 * This class implements the Forward Euler (RK-1) for single cycle and midpoint rule (RK-2) for
 * multiple cycles as the time-stepping schemes.
 *
 * Class INSVCStaggeredConservativeMassMomentumRKIntegrator computes the convective
 * derivative of a side-centered velocity field using various bounded-limiters
 * described by Patel and Natarajan.
 *
 * A side-centered density update is provided by this class, which is used in the
 * conservative discretization of the incompressible Navier-Stokes equation.
 * There is an optional mass density source term \f$ S(x,t) \f$ that can be set to check the order
 * of accuracy via manufactured solutions.
 *
 * References
 * Patel, JK. and Natarajan, G., <A HREF="https://www.sciencedirect.com/science/article/pii/S0045793014004009">
 * A generic framework for design of interface capturing schemes for multi-fluid flows</A>
 *
 * \note This class is specialized in that it computes a conservative discretization of the form
 * \f$N = \nabla \cdot (u \rho u)\f$, where the density \f$\rho\f$ can vary in space and time.
 * This operator is to be used in conjuction with the conservative form of the variable coefficient
 * Navier-Stokes equations, which will produce better results for high density ratio flows.
 *
 * \see INSVCStaggeredHierarchyIntegrator
 */
class INSVCStaggeredConservativeMassMomentumRKIntegrator : public STSMassFluxIntegrator
{
public:
    /*!
     * \brief Class constructor.
     */
    INSVCStaggeredConservativeMassMomentumRKIntegrator(std::string object_name,
                                                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual ~INSVCStaggeredConservativeMassMomentumRKIntegrator();

    /*!
     * \brief Integrate density and momentum field.
     */
    virtual void integrate(double dt) override;

    /*!
     * \name General operator functionality.
     */
    //\{

    /*!
     * \brief Compute hierarchy dependent data required for time integrating variables.
     */
    virtual void
    initializeSTSIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy) override;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSTSIntegrator().
     *
     * \note It is safe to call deallocateSTSIntegrator() when the time integrator
     * is already deallocated.
     *
     * \see initializeSTSIntegrator
     */
    virtual void deallocateSTSIntegrator() override;

    //\}
    /*
     * \brief Set the boundary condition object for the side-centered velocity.
     */
    void setSideCenteredVelocityBoundaryConditions(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_sc_bc_coefs);

    /*!
     * \brief Get the newly constructed side-centered density patch data index.
     *
     * \note This data is produced as a part of the apply() routine and should be used
     * in the linear operator for the INSVC solver.
     */
    int getUpdatedSideCenteredDensityPatchDataIndex();

    /*!
     * \brief Set an optional mass density source term.
     */
    void setMassDensitySourceTerm(const SAMRAI::tbox::Pointer<IBTK::CartGridFunction> S_fcn);

protected:
    /*!
     * \brief Compute the advection velocity using simple averages.
     */
    void
    computeAdvectionVelocity(std::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
                             const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > U_data,
                             const SAMRAI::hier::IntVector<NDIM>& patch_lower,
                             const SAMRAI::hier::IntVector<NDIM>& patch_upper,
                             const std::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes);

    /*!
     * \brief Compute the interpolation of a quantity Q onto Q_half, faces of the velocity DOF centered control volumes.
     */
    void interpolateSideQuantity(
        std::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> Q_half_data,
        const std::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > Q_data,
        const SAMRAI::hier::IntVector<NDIM>& patch_lower,
        const SAMRAI::hier::IntVector<NDIM>& patch_upper,
        const std::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes,
        const LimiterType& convective_limiter);

    /*!
     * \brief Compute div[rho_half*u_half*u_adv].
     */
    void computeConvectiveDerivative(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > N_data,
        std::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> P_half_data,
        const std::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const std::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> R_half_data,
        const std::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_half_data,
        const std::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes,
        const double* const dx,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch);

    /*!
     * \brief Compute the density update rho = a0*rho^0 + a1*rho^1 + a2*dt*(-div[u_adv*rho_half]) + a2*dt*S
     */
    void computeDensityUpdate(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > R_data,
        const double& a0,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > R0_data,
        const double& a1,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > R1_data,
        const double& a2,
        const std::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const std::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> R_half_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > S_data,
        const std::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes,
        const double& dt,
        const double* const dx);

    /*!
     * \brief Compute the error of the mass conservation equation using the integrated
     * density field pointwise.
     */
    void computeErrorOfMassConservationEquation(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > E_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > Rnew_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > Rold_data,
        const std::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const std::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> R_half_data,
        const std::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes,
        const double& dt,
        const double* const dx);
    /*!
     * \brief Enforce divergence free condition at the coarse-fine interface to ensure conservation of mass.
     */
    void enforceDivergenceFreeConditionAtCoarseFineInterface(const int U_idx);

    // Cached communications operators.
    std::string d_velocity_bdry_extrap_type = "CONSTANT";
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_v_transaction_comps;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_v_bdry_fill;

    // The limiter type for interpolation onto faces.
    LimiterType d_velocity_convective_limiter = CUI;

    // The required number of ghost cells for the chosen interpolation.
    int d_velocity_limiter_gcw = 1;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSVCStaggeredConservativeMassMomentumRKIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSVCStaggeredConservativeMassMomentumRKIntegrator(const INSVCStaggeredConservativeMassMomentumRKIntegrator& from) =
        delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSVCStaggeredConservativeMassMomentumRKIntegrator&
    operator=(const INSVCStaggeredConservativeMassMomentumRKIntegrator& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSVCStaggeredConservativeMassMomentumIntegrator
