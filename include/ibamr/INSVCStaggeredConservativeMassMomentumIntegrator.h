// Filename: INSVCStaggeredConservativeMassMomentumIntegrator.h
// Created on 01 April 2018 by Nishant Nangia and Amneet Bhalla
//
// Copyright (c) 2002-2018, Nishant Nangia and Amneet Bhalla
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBAMR_INSVCStaggeredConservativeMassMomentumIntegrator
#define included_IBAMR_INSVCStaggeredConservativeMassMomentumIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "CoarseFineBoundary.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "SideVariable.h"
#include "boost/array.hpp"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "tbox/Database.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

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
 * \brief Class INSVCStaggeredConservativeMassMomentumIntegrator integrates
 * the staggered density field
 *
 *  \f$ \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho u) = S(x,t) \f$
 *
 * and computes the conservative form of the convective operator
 * \f$ \nabla \cdot (\rho u u)\f$.
 *
 * Class INSVCStaggeredConservativeMassMomentumIntegrator computes the convective
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
class INSVCStaggeredConservativeMassMomentumIntegrator : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Class constructor.
     */
    INSVCStaggeredConservativeMassMomentumIntegrator(const std::string& object_name,
                                                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual ~INSVCStaggeredConservativeMassMomentumIntegrator();

    /*!
     * \brief Integrate density and momentum field.
     */
    void integrate(double dt);

    /*!
     * \name General operator functionality.
     */
    //\{

    /*!
     * \brief Compute hierarchy dependent data required for time integrating variables.
     */
    void initializeTimeIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeTimeIntegrator().
     *
     * \note It is safe to call deallocateTimeIntegrator() when the time integrator
     * is already deallocated.
     *
     * \see initializeTimeIntegrator
     */
    void deallocateTimeIntegrator();

    //\}

    /*!
     * \brief Set the current side-centered density patch data index.
     */
    void setSideCenteredDensityPatchDataIndex(int rho_sc_idx);

    /*!
     * \brief Set the new side-centered convective derivative patch data index.
     */
    void setSideCenteredConvectiveDerivativePatchDataIndex(int N_sc_idx);

    /*
     * \brief Set the boundary condition object for the side-centered velocity.
     */
    void setSideCenteredVelocityBoundaryConditions(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_sc_bc_coefs);

    /*
     * \brief Set the boundary condition object for the side-centered density.
     */
    void setSideCenteredDensityBoundaryConditions(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& rho_sc_bc_coefs);

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

    /*!
     * \brief Set the patch data indices corresponding to the velocity at the previous time step
     * to be used when computing the density update
     *
     * \note This velocities will be used to compute an approximation to velocities required for SSPRK updates
     * V_old_idx = n-1, V_current_idx = n, V_new_idx = n+1,k (after an INS cycle)
     * If V_old_idx or V_new_idx are not set, then they will degenerate to V_current automatically, for the very first
     * simulation time step and cases where an INS cycle has not been executed, respectively.
     */
    void setFluidVelocityPatchDataIndices(int V_old_idx, int V_current_idx, int V_new_idx);

    /*!
     * \brief Set the cycle number currently being executed by the INS integrator. This will determine the rho advection
     * velocity.
     */
    void setCycleNumber(int cycle_num);

    /*!
     * \brief Set solution time.
     */
    void setSolutionTime(double solution_time);

    /*!
     * \brief Set the current time interval.
     */
    void setTimeInterval(double current_time, double new_time);

    /*!
     * \brief Get the current time interval.
     */
    std::pair<double, double> getTimeInterval() const;

    /*!
     * \brief Get the current time step size.
     */
    double getDt() const;

    /*!
     * \brief Set the previous time step value between times n - 1 (old) and n (current). This is used during velocity
     * extrapolation.
     */
    void setPreviousTimeStepSize(double dt_prev);

    /*!
     * \brief Set the HierarchyMathOps object used by the operator.
     */
    void setHierarchyMathOps(SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops);

    /*!
     * \brief Get the HierarchyMathOps object used by the operator.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> getHierarchyMathOps() const;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSVCStaggeredConservativeMassMomentumIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSVCStaggeredConservativeMassMomentumIntegrator(const INSVCStaggeredConservativeMassMomentumIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSVCStaggeredConservativeMassMomentumIntegrator&
    operator=(const INSVCStaggeredConservativeMassMomentumIntegrator& that);

    /*!
     * \brief Compute the advection velocity using simple averages
     */
    void computeAdvectionVelocity(
        boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > U_data,
        const SAMRAI::hier::IntVector<NDIM>& patch_lower,
        const SAMRAI::hier::IntVector<NDIM>& patch_upper,
        const boost::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes);

    /*!
     * \brief Compute the interpolation of a quantity Q onto Q_half, faces of the velocity DOF centered control volumes
     */
    void interpolateSideQuantity(
        boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> Q_half_data,
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > Q_data,
        const SAMRAI::hier::IntVector<NDIM>& patch_lower,
        const SAMRAI::hier::IntVector<NDIM>& patch_upper,
        const boost::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes,
        const LimiterType& convective_limiter);

    /*!
     * \brief Compute div[rho_half*u_half*u_adv]
     */
    void computeConvectiveDerivative(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > N_data,
        boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> P_half_data,
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> R_half_data,
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_half_data,
        const boost::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes,
        const double* const dx);

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
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> U_adv_data,
        const boost::array<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceData<NDIM, double> >, NDIM> R_half_data,
        const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideData<NDIM, double> > S_data,
        const boost::array<SAMRAI::hier::Box<NDIM>, NDIM>& side_boxes,
        const double& dt,
        const double* const dx);

    /*!
     * \brief Enforce divergence free condition at the coarse-fine interface to ensure conservation of mass.
     */
    void enforceDivergenceFreeConditionAtCoarseFineInterface(const int U_idx);

    // Book keeping
    std::string d_object_name;

    // Boundary condition helper object.
    SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> d_bc_helper;

    // Cached communications operators.
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
    std::string d_velocity_bdry_extrap_type, d_density_bdry_extrap_type;
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_rho_transaction_comps;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_rho_bdry_fill;
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_v_transaction_comps;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_v_bdry_fill;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    // Number of RK steps to take.
    int d_num_steps;

    // Boundary condition object for side-centered velocity field.
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_u_sc_bc_coefs;

    // Boundary condition object for side-centered density field.
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_rho_sc_bc_coefs;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_V_var;
    int d_V_scratch_idx, d_V_old_idx, d_V_current_idx, d_V_new_idx, d_V_composite_idx, d_N_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_rho_sc_var;
    int d_rho_sc_current_idx, d_rho_sc_scratch_idx, d_rho_sc_new_idx;

    // Hierarchy operation objects.
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> > d_hier_sc_data_ops;

    // Mathematical operators.
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;
    bool d_hier_math_ops_external;

    // Boolean value to indicate whether the integrator is presently
    // initialized.
    bool d_is_initialized;

    // Logging configuration.
    bool d_enable_logging;

    // The limiter type for interpolation onto faces.
    LimiterType d_velocity_convective_limiter;
    LimiterType d_density_convective_limiter;

    // The required number of ghost cells for the chosen interpolation.
    int d_velocity_limiter_gcw, d_density_limiter_gcw;

    // Variable to indicate the density update time-stepping type.
    TimeSteppingType d_density_time_stepping_type;

    // Source term variable and function for the mass density update.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_S_var;
    int d_S_scratch_idx;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_S_fcn;

    // Variable to indicate the cycle number.
    int d_cycle_num;

    // Variable to indicate time and time step size.
    double d_current_time, d_new_time, d_solution_time, d_dt, d_dt_prev;

    // Coarse-fine boundary object
    std::vector<SAMRAI::hier::CoarseFineBoundary<NDIM>*> d_cf_boundary;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_INSVCStaggeredConservativeMassMomentumIntegrator
