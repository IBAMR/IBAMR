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

#ifndef included_IBAMR_MassIntegrator
#define included_IBAMR_MassIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

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
 * \brief Class MassIntegrator is an abstract class which integrates
 * the density field
 *
 *  \f$ \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho u) = S(x,t) \f$
 *
 * and computes the conservative form of the convective operator
 * \f$ \nabla \cdot (\rho u u)\f$ for momentum and \f$ \nabla \cdot (\rho u Q)\f$ for
 * energy equation where \f$ Q = C_\textrm{p} T\f$ or \f$ Q = h\f$ based on the form of the
 * energy equation.
 */
class MassIntegrator : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Class constructor.
     */
    MassIntegrator(std::string object_name, SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    virtual ~MassIntegrator() = default;

    /*!
     * \brief Integrate density and compute the convective operator of the momentum and/or energy.
     */
    virtual void integrate(double dt) = 0;

    /*!
     * \name General operator functionality.
     */
    //\{

    /*!
     * \brief Compute hierarchy dependent data required for time integrating variables.
     */
    virtual void
    initializeTimeIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy) = 0;

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeTimeIntegrator().
     *
     * \note It is safe to call deallocateTimeIntegrator() when the time integrator
     * is already deallocated.
     *
     * \see initializeTimeIntegrator
     */
    virtual void deallocateTimeIntegrator() = 0;

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

    // Number of RK steps to take. Default is set for RK3.
    int d_num_steps = 3;

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
    MassIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    MassIntegrator(const MassIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    MassIntegrator& operator=(const MassIntegrator& that) = delete;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_MassIntegrator
