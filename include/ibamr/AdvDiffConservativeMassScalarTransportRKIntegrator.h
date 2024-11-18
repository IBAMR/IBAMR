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

#ifndef included_IBAMR_AdvDiffConservativeMassScalarTransportRKIntegrator
#define included_IBAMR_AdvDiffConservativeMassScalarTransportRKIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////
#include <ibamr/config.h>

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
 * \brief Class AdvDiffConservativeMassScalarTransportRKIntegrator is a concrete class that integrates
 * the collocated density field
 *
 *  \f$ \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho u) = S(x,t) \f$
 *
 * and computes the conservative form of the convective operator
 * \f$ \nabla \cdot (\rho \gamma u Q)\f$ where \f$\gamma\f$ is a material property
 * and \f$Q\f$ is a transport quantity that is being transported. If \f$\gamma\f$ is not
 * registered to this class, then this class computes \f$ \nabla \cdot (\rho u Q)\f$.
 *
 * This class implements the Forward Euler (RK-1) for single cycle and midpoint rule (RK-2) for
 * multiple cycles as the time-stepping schemes.
 *
 * This class computes the convective derivative of a cell-centered transport quantity using
 * various bounded-limiters described by Patel and Natarajan.
 *
 * References
 * Patel, JK. and Natarajan, G., <A HREF="https://www.sciencedirect.com/science/article/pii/S0045793014004009">
 * A generic framework for design of interface capturing schemes for multi-fluid flows</A>
 *
 * A cell-centered density update is provided by this class, which is used in the
 * conservative discretization of the energy equation.
 *
 * \note This class is specialized in that it computes a conservative discretization of the form
 * \f$N = \nabla \cdot (u \rho \gamma Q)\f$, where the density \f$\rho\f$ can vary in space and time.
 * This operator is to be used in conjuction with the conservative form of the variable coefficient
 * energy equations.
 *
 */
class AdvDiffConservativeMassScalarTransportRKIntegrator : public STSMassFluxIntegrator
{
public:
    /*!
     * \brief Class constructor.
     */
    AdvDiffConservativeMassScalarTransportRKIntegrator(std::string object_name,
                                                       SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~AdvDiffConservativeMassScalarTransportRKIntegrator();

    /*!
     * \brief Integrate density and momentum field.
     */
    void integrate(double dt) override;

    /*!
     * \name General operator functionality.
     */
    //\{

    /*!
     * \brief Compute hierarchy dependent data required for time integrating variables.
     */
    void
    initializeSTSIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy) override;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeTimeIntegrator().
     *
     * \note It is safe to call deallocateTimeIntegrator() when the time integrator
     * is already deallocated.
     *
     * \see initializeTimeIntegrator
     */
    void deallocateSTSIntegrator() override;

    //\}
    /*!
     * \brief Set the current cell-centered material property patch data index.
     */
    void setCellCenteredMaterialPropertyPatchDataIndex(int gamma_cc_idx);

    /*!
     * \brief Set the current cell-centered transport quantity patch data index.
     */
    void setCellCenteredTransportQuantityPatchDataIndex(int Q_cc_idx);

    /*!
     * \brief Set the boundary condition object for the cell-centered density.
     */
    void setCellCenteredDensityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& rho_cc_bc_coefs);

    /*!
     * \brief Set the boundary condition object for the cell-centered material property.
     */
    void setCellCenteredMaterialPropertyBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& gamma_cc_bc_coefs);

    /*!
     * \brief Set the boundary condition object for the cell-centered transport quantity.
     */
    void setCellCenteredTransportQuantityBoundaryConditions(SAMRAI::solv::RobinBcCoefStrategy<NDIM>*& Q_cc_bc_coefs);

    /*!
     * \brief Get the newly constructed cell-centered density patch data index.
     */
    int getUpdatedCellCenteredDensityPatchDataIndex();

    /*!
     * \brief Set the patch data indices corresponding to the material property
     * to be used when computing the convective derivative
     *
     * \note This material property will be used to compute an approximation to material property required for computing
     * convective derivative. gamma_current_idx = n, gamma_new_idx = n+1,k (after a cycle). If
     * gamma_new_idx is not set, then it will degenerate to gamma_current_idx automatically, for the
     * very first simulation time step and cases where an AdvDiff cycle has not been executed, respectively.
     */
    void setMaterialPropertyPatchDataIndices(int gamma_current_idx, int gamma_new_idx);

    /*!
     * \brief Set the patch data indices corresponding to the transport quantity
     * to be used when computing the convective derivative.
     *
     * \note This transport quantity will be used to compute an approximation to transport quantity required for
     * computing convective derivative. Q_current_idx = n, Q_new_idx = n+1,k (after a
     * cycle). If Q_new_idx is not set, then it will degenerate to Q_current_idx automatically, for the very
     * first simulation time step and cases where an AdvDiff cycle has not been executed, respectively.
     */
    void setTransportQuantityPatchDataIndices(int Q_current_idx, int Q_new_idx);

    /*!
     * \brief Set the material property variable.
     */
    void setMaterialPropertyVariable(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > gamma_var);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffConservativeMassScalarTransportRKIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffConservativeMassScalarTransportRKIntegrator(const AdvDiffConservativeMassScalarTransportRKIntegrator& from) =
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
    AdvDiffConservativeMassScalarTransportRKIntegrator&
    operator=(const AdvDiffConservativeMassScalarTransportRKIntegrator& that) = delete;

    /*!
     * \brief Compute the interpolation of a quantity Q onto faces of the cell centered control volumes
     */
    void interpolateCellQuantity(SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > Q_half_data,
                                 SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > U_adv_data,
                                 const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > Q_data,
                                 const SAMRAI::hier::IntVector<NDIM>& patch_lower,
                                 const SAMRAI::hier::IntVector<NDIM>& patch_upper,
                                 const SAMRAI::hier::Box<NDIM>& patch_box,
                                 const LimiterType& convective_limiter);

    /*!
     * \brief Compute div[rho_lim * gamma_lim * u_adv * Q_lim]. Here, rho_lim * u_adv is obtained from
     * integrating the mass balance equation.
     */
    void computeConvectiveDerivative(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > N_data,
                                     SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > P_half_data,
                                     const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > U_adv_data,
                                     const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > R_half_data,
                                     const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > Q_half_data,
                                     const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > G_half_data,
                                     const SAMRAI::hier::Box<NDIM>& patch_box,
                                     const double* const dx);

    /*!
     * \brief Compute the density update rho = a0*rho^0 + a1*rho^1 + a2*dt*(-div[u_adv*rho_half]) + a2*dt*S
     */
    void computeDensityUpdate(SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > R_data,
                              const double& a0,
                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > R0_data,
                              const double& a1,
                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > R1_data,
                              const double& a2,
                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > U_adv_data,
                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > R_half_data,
                              const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > S_data,
                              const SAMRAI::hier::Box<NDIM>& patch_box,
                              const double& dt,
                              const double* const dx);

    /*!
     * \brief Compute the error of the mass conservation equation using the integrated
     * density field pointwise.
     */
    void computeErrorOfMassConservationEquation(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > E_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > Rnew_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM, double> > Rold_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > U_adv_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> > R_half_data,
        const SAMRAI::hier::Box<NDIM>& patch_box,
        const double& dt,
        const double* const dx);

    // Book keeping
    std::string d_object_name;

    std::string d_transport_quantity_bdry_extrap_type = "CONSTANT", d_material_property_bdry_extrap_type = "CONSTANT";

    // Cached communications operators.
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_Q_transaction_comps,
        d_gamma_transaction_comps;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_Q_bdry_fill, d_hier_gamma_bdry_fill;

    /*
     * Boundary condition objects
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_rho_cc_bc_coefs;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_gamma_cc_bc_coefs;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_Q_cc_bc_coefs;

    /*!
     * Variables.
     */

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_rho_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_gamma_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Q_cc_var;

    /*
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_gamma_cc_current_idx = IBTK::invalid_index, d_gamma_cc_scratch_idx = IBTK::invalid_index,
        d_gamma_cc_new_idx = IBTK::invalid_index;
    int d_Q_cc_current_idx = IBTK::invalid_index, d_Q_cc_scratch_idx = IBTK::invalid_index,
        d_Q_cc_new_idx = IBTK::invalid_index;

    /*
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_gamma_cc_composite_idx = IBTK::invalid_index;
    int d_Q_cc_composite_idx = IBTK::invalid_index;

    // The limiter type for interpolation onto faces.
    LimiterType d_transport_quantity_convective_limiter = CUI;
    LimiterType d_material_property_convective_limiter = CUI;

    // The required number of ghost cells for the chosen interpolation.
    int d_transport_quantity_limiter_gcw = 1, d_material_property_limiter_gcw = 1;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_AdvDiffConservativeMassScalarTransportRKIntegrator
