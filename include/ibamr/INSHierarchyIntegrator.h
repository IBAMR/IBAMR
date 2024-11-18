// ---------------------------------------------------------------------
//
// Copyright (c) 2006 - 2024 by the IBAMR developers
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

#ifndef included_IBAMR_INSHierarchyIntegrator
#define included_IBAMR_INSHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/PoissonSolver.h"

#include "FaceVariable.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "Variable.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Patch;
template <int DIM>
class PatchLevel;
} // namespace hier
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
 * \brief Class INSHierarchyIntegrator provides an abstract interface for a time
 * integrator for the incompressible Navier-Stokes equations on an AMR grid
 * hierarchy, along with basic data management for variables defined on that
 * hierarchy.
 */
class INSHierarchyIntegrator : public IBTK::HierarchyIntegrator
{
public:
    /*!
     * The destructor for class INSHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~INSHierarchyIntegrator();

    /*!
     * Set the type of viscous time integration scheme being employed by the
     * incompressible flow solver.
     *
     * Different implementations may support different time stepping schemes.
     */
    void setViscousTimeSteppingType(TimeSteppingType viscous_time_stepping_type);

    /*!
     * Get the type of viscous time integration scheme being employed by the
     * incompressible flow solver.
     *
     * Different implementations may support different time stepping schemes.
     */
    TimeSteppingType getViscousTimeSteppingType() const;

    /*!
     * Set the type of convective time integration scheme being employed by the
     * incompressible flow solver.
     *
     * Different implementations may support different time stepping schemes.
     */
    void setConvectiveTimeSteppingType(TimeSteppingType convective_time_stepping_type);

    /*!
     * Get the type of convective time integration scheme being employed by the
     * incompressible flow solver.
     *
     * Different implementations may support different time stepping schemes.
     */
    TimeSteppingType getConvectiveTimeSteppingType() const;

    /*!
     * Set the type of convective time integration scheme being employed by the
     * incompressible flow solver during the initial time step.
     *
     * Different implementations may support different time stepping schemes.
     *
     * \note This is used \em only when the basic convective time stepping
     * scheme uses a multi-step method such as Adams-Bashforth.
     */
    void setInitialConvectiveTimeSteppingType(TimeSteppingType init_convective_time_stepping_type) const;

    /*!
     * Get the type of convective time integration scheme being employed by the
     * incompressible flow solver during the initial time step.
     *
     * Different implementations may support different time stepping schemes.
     *
     * \note This is used \em only when the basic convective time stepping
     * scheme uses a multi-step method such as Adams-Bashforth.
     */
    TimeSteppingType getInitialConvectiveTimeSteppingType() const;

    /*!
     * Register an advection-diffusion solver with this incompressible
     * Navier-Stokes solver.
     */
    void registerAdvDiffHierarchyIntegrator(SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator);

    /*!
     * Set the problem coefficients used by the solver.
     */
    void setStokesSpecifications(StokesSpecifications problem_coefs);

    /*!
     * Get a const pointer to the problem coefficients object used by the
     * solver.
     */
    const StokesSpecifications* getStokesSpecifications() const;

    /*!
     * Supply a physical boundary conditions specificaion for the velocity
     * field. Boundary conditions take the form of \f$ a\mathbf{u} + b\tau\cdot\mathbf{n} = \mathbf{g}\f$ where \f$\tau
     * = -p\mathbf{I} + \frac{\mu}{2}\left(\nabla\mathbf{u} + \nabla\mathbf{u}^T\right)\f$ is the Newtonian fluid
     * stress and \f$\mu\f$ is the fluid viscosity.
     *
     * \note Periodic boundaries take presidence over physical boundaries.
     * The type of boundary can be adjusted through the `CartesianGeometry` database in the
     * input file, specifically with the flag `periodic_dimension`.
     *
     * \note Current implementations of physical boundary conditions require
     * that at any physical location the values of \f$ a\f$ and \f$ b \f$
     * stasify that either \f$ a = 1 \f$ or \f$b = 1\f$ as well as
     * \f$a + b = 1\f$. An error will occur if this does not happen.
     *
     * \see IBTK::muParserRobinBcCoefs
     */
    void registerPhysicalBoundaryConditions(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * Get a vector of pointers to the velocity boundary condition specification
     * objects.
     *
     * \note Implementations may return a vector of nullptr pointers.
     */
    virtual const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& getVelocityBoundaryConditions() const;

    /*!
     * Get a pointer to the pressure boundary condition specification object.
     *
     * \note Implementations may return a nullptr pointer.
     */
    virtual SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getPressureBoundaryConditions() const;

    /*!
     * Supply initial conditions for the velocity field.
     */
    void registerVelocityInitialConditions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> U_init);

    /*!
     * Supply initial conditions for the pressure.
     *
     * \note Initial conditions are not required for the pressure, but when
     * available, they can speed the convergence of the solver during the
     * initial time step.
     */
    void registerPressureInitialConditions(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> P_init);

    /*!
     * Supply a body force.
     */
    void registerBodyForceFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * Supply a fluid source/sink distribution.
     *
     * @deprecated Use registerVelocityDivergenceFunction() instead.
     */
    void registerFluidSourceFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> Q_fcn);

    /*!
     * Supply a CartGridFunction that specifies \f$ \nabla \cdot \mathbf{u} \f$.
     *
     * @note If \param Q_fcn is not specified, then \f$ \nabla \cdot \mathbf{u}  = 0 \f$ is imposed.
     */
    void registerVelocityDivergenceFunction(SAMRAI::tbox::Pointer<IBTK::CartGridFunction> Q_fcn);

    /*!
     * Return a pointer to the fluid velocity variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getVelocityVariable() const;

    /*!
     * Return a pointer to the fluid pressure state variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getPressureVariable() const;

    /*!
     * Return a pointer to the body force variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getBodyForceVariable() const;

    /*!
     * Return a pointer to the source strength variable.
     *
     * @deprecated Use getVelocityDivergenceVariable() instead.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getFluidSourceVariable() const;

    /*!
     * Return a pointer to the variable that specifies the divergence of the velocity.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > getVelocityDivergenceVariable() const;

    /*!
     * Return a pointer to a fluid velocity variable that can be used to advect
     * quantities registered with an advection-diffusion solver.
     *
     * Data are allocated for this variable using the current context.  Patch
     * data for this variable are allocated only when an advection-diffusion
     * solver is registered with the Navier-Stokes solver.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > getAdvectionVelocityVariable() const;

    /*!
     * Get a vector of pointers to the intermediate velocity boundary condition
     * specification objects.
     */
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> getIntermediateVelocityBoundaryConditions() const;

    /*!
     * Get a pointer to the projection Poisson problem boundary condition
     * specification object.
     */
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* getProjectionBoundaryConditions() const;

    /*!
     * \brief Set the convective operator type to be used by the solver.
     */
    void setConvectiveOperatorType(const std::string& op_type);

    /*!
     * \brief Get the convective operator type used by the solver.
     */
    const std::string& getConvectiveOperatorType() const;

    /*!
     * \brief Set the convective differencing form to be used by the
     * solver.
     */
    void setConvectiveDifferencingType(ConvectiveDifferencingType difference_form);

    /*!
     * \brief Get the convective differencing form used by the solver.
     */
    ConvectiveDifferencingType getConvectiveDifferencingType() const;

    /*!
     * \brief Set whether the integrator solves the time-dependent (creeping)
     * Stokes equations.  Otherwise, the integrator solves the time-dependent
     * Navier-Stokes equations.
     */
    void setCreepingFlow(bool creeping_flow);

    /*!
     * \brief Get whether the integrator solves the time-dependent (creeping)
     * Stokes equations.  Otherwise, the integrator solves the time-dependent
     * Navier-Stokes equations.
     */
    bool getCreepingFlow() const;

    /*!
     * Register an operator to compute the convective acceleration term u*grad
     * u.
     *
     * If the supplied operator is nullptr, then the integrator will solve the
     * time-dependent (creeping) Stokes equations instead of the Navier-Stokes
     * equations.
     */
    void setConvectiveOperator(SAMRAI::tbox::Pointer<ConvectiveOperator> convective_op);

    /*!
     * Get the convective operator being used by this solver class.
     *
     * If the time integrator is configured to solve the time-dependent
     * (creeping) Stokes equations, then the returned pointer will be nullptr.
     *
     * If the convective operator has not already been constructed, then this
     * function will initialize a default convective operator.
     */
    virtual SAMRAI::tbox::Pointer<ConvectiveOperator> getConvectiveOperator() = 0;

    /*!
     * Indicate that the convective operator should be (re-)initialized before
     * the next time step.
     */
    void setConvectiveOperatorNeedsInit();

    /*!
     * Register a solver for the velocity subsystem.
     */
    void setVelocitySubdomainSolver(SAMRAI::tbox::Pointer<IBTK::PoissonSolver> velocity_solver);

    /*!
     * Get the subdomain solver for the velocity subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     *
     * If the velocity subdomain solver has not already been constructed, then
     * this function will initialize a default solver.
     */
    virtual SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getVelocitySubdomainSolver() = 0;

    /*!
     * Indicate that the velocity subdomain solver should be (re-)initialized
     * before the next time step.
     */
    void setVelocitySubdomainSolverNeedsInit();

    /*!
     * Register a solver for the pressure subsystem.
     */
    void setPressureSubdomainSolver(SAMRAI::tbox::Pointer<IBTK::PoissonSolver> pressure_solver);

    /*!
     * Get the subdomain solver for the pressure subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     *
     * If the pressure subdomain solver has not already been constructed, then
     * this function will initialize a default solver.
     */
    virtual SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getPressureSubdomainSolver() = 0;

    /*!
     * Indicate that the velocity subdomain solver should be (re-)initialized
     * before the next time step.
     */
    void setPressureSubdomainSolverNeedsInit();

    /*!
     * Returns the number of cycles to perform for the present time step.
     */
    int getNumberOfCycles() const override;

    /*!
     * Finish postprocessing the hierarchy by computing the current CFL number.
     */
    virtual void postprocessIntegrateHierarchy(double current_time,
                                               double new_time,
                                               bool skip_synchronize_new_state_data,
                                               int num_cycles = 1) override;

    /*!
     * Return the current CFL number (i.e., the CFL number computed from the
     * current velocity field). This is typically computed at the end of each
     * time step.
     */
    virtual double getCurrentCFLNumber() const;

protected:
    /*!
     * The constructor for class INSHierarchyIntegrator sets some default
     * values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     *
     * This constructor sets the default coarsen and refine operator types to:
     *
     * - "CONSERVATIVE_COARSEN" and "CONSERVATIVE_LINEAR_REFINE" for U
     * - "CONSERVATIVE_COARSEN" and "LINEAR_REFINE" for P
     * - "CONSERVATIVE_COARSEN" and "CONSERVATIVE_LINEAR_REFINE" for F
     * - "CONSERVATIVE_COARSEN" and "CONSTANT_REFINE" for Q
     *
     * The other constructor allows these default values to be overridden.
     */
    INSHierarchyIntegrator(std::string object_name,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > U_var,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > P_var,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > F_var,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > Q_var,
                           bool register_for_restart);

    /*!
     * The constructor for class INSHierarchyIntegrator sets some default
     * values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    INSHierarchyIntegrator(std::string object_name,
                           SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > U_var,
                           std::string U_default_coarsen_type,
                           std::string U_default_refine_type,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > P_var,
                           std::string P_default_coarsen_type,
                           std::string P_default_refine_type,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > F_var,
                           std::string F_default_coarsen_type,
                           std::string F_default_refine_type,
                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > Q_var,
                           std::string Q_default_coarsen_type,
                           std::string Q_default_refine_type,
                           bool register_for_restart);

    /*!
     * Pure virtual method to project the velocity field following a regridding operation.
     */
    virtual void regridProjection() = 0;

    /*!
     * Update the current CFL number (i.e., at the end of a timestep).
     *
     * @note this method can handle both cell-centered and side-centered
     * velocities.
     */
    virtual void updateCurrentCFLNumber(const int data_idx, const double dt);

    /*!
     * Compute the maximum vorticity magnitude at any given point.
     *
     * @note This function does not read ghost data from @p Omega_idx.
     */
    double getMaximumVorticityMagnitude(const int Omega_idx);

    /*!
     * Tag cells based on the vorticity magnitude.
     *
     * @note This function is typically called by applyGradientDetectorSpecialized() in inheriting classes.
     */
    void tagCellsByVorticityMagnitude(const int level_number, const int Omega_idx, const int tag_idx);

    /*!
     * Return the maximum stable time step size.
     */
    double getMaximumTimeStepSizeSpecialized() override;

    /*!
     * Determine the largest stable timestep on an individual patch level.
     */
    double getStableTimestep(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level) const;

    /*!
     * Determine the largest stable timestep on an individual patch.
     */
    virtual double getStableTimestep(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const = 0;

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*
     * Boolean value that indicates whether the integrator has been initialized.
     */
    bool d_integrator_is_initialized = false;

    /*!
     * Enum indicating the time integration employed for the implicit
     * discretization of the viscous terms.
     */
    TimeSteppingType d_viscous_time_stepping_type = TRAPEZOIDAL_RULE;

    /*!
     * Enum indicating the time integration employed for the explicit
     * discretization of the convective terms.
     */
    TimeSteppingType d_convective_time_stepping_type = ADAMS_BASHFORTH;

    /*!
     * Enum indicating the time integration employed for the explicit
     * discretization of the convective terms during the \em initial time step.
     */
    TimeSteppingType d_init_convective_time_stepping_type = MIDPOINT_RULE;

    /*!
     * Problem coeficients.
     */
    StokesSpecifications d_problem_coefs;

    /*
     * The AdvDiffHierarchyIntegrator is used to provide time integration
     * capability for quantities transported by the fluid velocity field.
     */
    std::vector<SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> > d_adv_diff_hier_integrators;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM, double> > d_U_adv_diff_var;

    /*!
     * Current CFL number.
     */
    double d_cfl_current = std::numeric_limits<double>::quiet_NaN();

    /*!
     * The maximum CFL number.
     */
    double d_cfl_max = 1.0;

    /*!
     * Cell tagging criteria based on the relative and absolute magnitudes of
     * the local vorticity.
     */
    bool d_using_vorticity_tagging = false;
    SAMRAI::tbox::Array<double> d_Omega_rel_thresh, d_Omega_abs_thresh;

    /*!
     * This boolean value determines whether the pressure is normalized to have
     * zero mean (i.e., discrete integral) at the end of each timestep.
     */
    bool d_normalize_pressure = false;

    /*!
     * This boolean value determines whether the velocity is normalized to have
     * zero mean (i.e., discrete integral) at the end of each timestep.
     *
     * This parameter only affects the case in which rho=0 (i.e. the steady
     * Stokes equations).
     */
    bool d_normalize_velocity = false;

    /*!
     * This boolean value determines whether the convective acceleration term is
     * included in the momentum equation.  (If it is not, this solver
     * effectively solves the so-called creeping Stokes equations.)
     */
    bool d_creeping_flow = false;

    /*!
     * Threshold that determines whether the velocity field needs to be
     * reprojected following adaptive regridding.
     */
    double d_regrid_max_div_growth_factor = 1.1;

    /*!
     * Double precision values are (optional) factors used to rescale the
     * velocity, pressure, and force for plotting.
     *
     * Boolean values indicates whether to output various quantities for
     * visualization.
     */
    double d_U_scale = 1.0, d_P_scale = 1.0, d_F_scale = 1.0, d_Q_scale = 1.0, d_Omega_scale = 1.0, d_Div_U_scale = 1.0,
           d_EE_scale = 1.0;
    bool d_output_U = true, d_output_P = true, d_output_F = false, d_output_Q = false, d_output_Omega = true,
         d_output_Div_U = true, d_output_EE = false;

    /*!
     * Fluid solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_U_var;
    std::string d_U_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_U_refine_type = "CONSERVATIVE_LINEAR_REFINE";

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_P_var;
    std::string d_P_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_P_refine_type = "LINEAR_REFINE";

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_F_var;
    std::string d_F_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_F_refine_type = "CONSERVATIVE_LINEAR_REFINE";

    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_Q_var;
    std::string d_Q_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_Q_refine_type = "CONSTANT_REFINE";

    /*!
     * Objects to set initial conditions, boundary conditions, body forces, and
     * fluid source/sink distributions.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_U_init, d_P_init;
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM> d_default_bc_coefs;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs, d_U_bc_coefs, d_U_star_bc_coefs;
    TractionBcType d_traction_bc_type = TRACTION;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef;
    std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM> > d_Phi_bc_coef;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_F_fcn, d_Q_fcn;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_U_bdry_bc_fill_op, d_P_bdry_bc_fill_op,
        d_Q_bdry_bc_fill_op, d_no_fill_op;

    bool d_use_div_sink_drag_term = false;

    /*!
     * Hierarchy operators and solvers and related configuration data.
     */
    int d_coarsest_reset_ln, d_finest_reset_ln;

    std::string d_convective_op_type = "DEFAULT";
    ConvectiveDifferencingType d_convective_difference_form = ADVECTIVE;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_convective_op_input_db;
    SAMRAI::tbox::Pointer<ConvectiveOperator> d_convective_op;
    bool d_convective_op_needs_init;

    std::string d_velocity_solver_type, d_velocity_precond_type, d_velocity_sub_precond_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_velocity_solver_db, d_velocity_precond_db,
        d_velocity_sub_precond_db;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_velocity_solver;
    bool d_velocity_solver_needs_init;

    std::string d_pressure_solver_type, d_pressure_precond_type, d_pressure_sub_precond_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_pressure_solver_db, d_pressure_precond_db,
        d_pressure_sub_precond_db;
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> d_pressure_solver;
    bool d_pressure_solver_needs_init;

    std::string d_regrid_projection_solver_type, d_regrid_projection_precond_type, d_regrid_projection_sub_precond_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_regrid_projection_solver_db, d_regrid_projection_precond_db,
        d_regrid_projection_sub_precond_db;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSHierarchyIntegrator(const INSHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSHierarchyIntegrator& operator=(const INSHierarchyIntegrator& that) = delete;

    /*!
     * Read input values from a given database.
     */
    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSHierarchyIntegrator
