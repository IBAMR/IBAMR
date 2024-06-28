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

#ifndef included_IBAMR_STSMassFluxIntegrator
#define included_IBAMR_STSMassFluxIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/ibtk_utilities.h"

#include "CoarseFineBoundary.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "SideVariable.h"
#include "tbox/Database.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

#include <array>
#include <string>
#include <vector>

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
 * \brief Class STSMassFluxIntegrator is an abstract class which integrates
 * the density field
 *
 *  \f$ \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{u}) = S(x,t) \f$
 *
 * and computes the conservative form of the convective operator
 * \f$ \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u})\f$ for the momentum and
 * \f$ \nabla \cdot (\rho \mathbf{u} Q)\f$ for the energy equation where \f$ Q = C_\textrm{p} T\f$
 * or \f$ Q = h\f$ based on the form of the energy equation.
 *
 * \note STSMassFluxIntegrator is an intermediate or a lightweight integrator that is used within one or a single time
 * step (STS) to get a consistent mass flux \f$ m_\rho = \rho \mathbf{u} \f$ for advecting either momentum \f$
 * \mathbf{u} \f$ through \f$ \nabla \cdot (\rho \mathbf{u} \otimes \mathbf{u})  = \nabla \cdot (m_\rho \otimes
 * \mathbf{u} ) \f$ or a scalar \f$ Q \f$ through  \f$ \nabla \cdot (\rho \mathbf{u} Q)  = \nabla \cdot (m_\rho  Q) \f$.
 * It does not maintain any hierarchy data, child or parent integrator, as done in implementations of
 * HierarchyIntegrator.
 *
 * The core concept behind STSMassFluxIntegrator is to obtain a consistent mass flux $\mathbf{m_\rho} = \rho \mathbf{u}$
 * for stabilizing high density ratio flows. It also ensures a consistency condition that:
 *
 * 1. When there are no viscous, pressure forces or any other body forces in the domain, the momentum equation becomes
 *    discretely equal to the mass equation.
 * 2. Under isothermal/isenthalpic conditions, with no additional heat source, the energy/enthalpy equation reverts to
 *    mass balance equation.
 *
 * To achieve these conditions, the time integrator scheme of STSMassFluxIntegrator is tightly coupled to how
 * INSVCStaggeredHierarchyIntegrator or PhaseChangeHierarchyIntegrator integrate their respective variables.
 * Note that when flow is incompressible, one does not need to solve mass balance equation.
 * Here, we are solving a redundant mass balance equation through STSMassFluxIntegrator class in order to achieve
 * consistency and stability.
 *
 * Reference
 * Nangia et. al, <A HREF="https://www.sciencedirect.com/science/article/pii/S0021999119302256">
 * A robust incompressible Navier-Stokes solver for high density ratio multiphase flows </A>
 *
 * Thirumalaisamy and Bhalla, <A HREF="https://arxiv.org/abs/2301.06256">
 * A low Mach enthalpy method to model non-isothermal gas-liquid-solid flows with melting and solidification </A>
 */
class STSMassFluxIntegrator : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Class constructor.
     */
    STSMassFluxIntegrator(std::string object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual ~STSMassFluxIntegrator() = default;

    /*!
     * \brief Integrate density and compute the convective operator of the momentum and/or energy.
     */
    virtual void integrate(double dt) = 0;

    /*!
     * \name General time stepping functionality.
     */
    //\{

    /*!
     * \brief Compute hierarchy dependent data required for time integrating variables.
     */
    virtual void
    initializeSTSIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy) = 0;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSTSIntegrator().
     *
     * \note It is safe to call deallocateSTSIntegrator() when the time integrator
     * is already deallocated.
     *
     * \see initializeSTSIntegrator
     */
    virtual void deallocateSTSIntegrator() = 0;

    //\}

    /*!
     * \brief Set the current cell-centered density patch data index.
     */
    void setDensityPatchDataIndex(int rho_idx);

    /*!
     * \brief Set the new convective derivative patch data index.
     */
    void setConvectiveDerivativePatchDataIndex(int N_idx);

    /*!
     * \brief Set the boundary condition object for the density.
     */
    void setDensityBoundaryConditions(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& rho_sc_bc_coefs);

    /*!
     * \brief Get the newly constructed density patch data index.
     */
    int getUpdatedDensityPatchDataIndex();

    /*!
     * \brief Set the patch data indices corresponding to the velocity at the previous time step
     * to be used when computing the density update
     *
     * \note This velocities will be used to compute an approximation to velocities required for time integrator.
     * For example, SSPRK updates V_old_idx = n-1, V_current_idx = n, V_new_idx = n+1,k (after an INS cycle).
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
    double getTimeStepSize() const;

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

protected:
    // Book keeping
    std::string d_object_name;

    // Boundary condition helper object.
    SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> d_bc_helper;

    // Cached communications operators.
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
    std::string d_density_bdry_extrap_type = "CONSTANT";
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_rho_transaction_comps;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_rho_bdry_fill;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    // Boundary condition object for velocity field.
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_u_bc_coefs;

    // Boundary condition object for side-centered density field.
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_rho_bc_coefs;

    // Hierarchy operation objects.
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM, double> > d_hier_fc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> > d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> > d_hier_cc_data_ops;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_V_var;
    int d_V_scratch_idx = IBTK::invalid_index, d_V_old_idx = IBTK::invalid_index, d_V_current_idx = IBTK::invalid_index,
        d_V_new_idx = IBTK::invalid_index, d_V_composite_idx = IBTK::invalid_index;

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_rho_var;
    int d_rho_current_idx = IBTK::invalid_index, d_rho_scratch_idx = IBTK::invalid_index,
        d_rho_new_idx = IBTK::invalid_index, d_rho_composite_idx = IBTK::invalid_index;
    int d_N_idx = IBTK::invalid_index;

    // Source term variable and function for the mass density update.
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_S_var;
    int d_S_scratch_idx = IBTK::invalid_index;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_S_fcn;

    // Variable and index to store the error of mass conservation equation.
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_E_var;
    int d_E_scratch_idx = IBTK::invalid_index;

    // Mathematical operators.
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;
    bool d_hier_math_ops_external = false;

    // Boolean value to indicate whether the integrator is presently
    // initialized.
    bool d_is_initialized = false;

    // Logging configuration.
    bool d_enable_logging = false;

    // The limiter type for interpolation onto faces.
    LimiterType d_density_convective_limiter = CUI;

    // The required number of ghost cells for the chosen interpolation.
    int d_density_limiter_gcw = 1;

    // Variable to indicate the density update time-stepping type.
    TimeSteppingType d_density_time_stepping_type = FORWARD_EULER;

    // Variable to indicate the cycle number.
    int d_cycle_num = -1;

    // Variable to indicate time and time step size.
    double d_current_time = std::numeric_limits<double>::quiet_NaN(),
           d_new_time = std::numeric_limits<double>::quiet_NaN(),
           d_solution_time = std::numeric_limits<double>::quiet_NaN(), d_dt = std::numeric_limits<double>::quiet_NaN(),
           d_dt_prev = -1.0;

    // Coarse-fine boundary objects.
    std::vector<SAMRAI::hier::CoarseFineBoundary<NDIM> > d_cf_boundary;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    STSMassFluxIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    STSMassFluxIntegrator(const STSMassFluxIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    STSMassFluxIntegrator& operator=(const STSMassFluxIntegrator& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_STSMassFluxIntegrator
