// Filename: INSStaggeredHierarchyIntegrator.h
// Created on 20 Mar 2008 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
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
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
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

#ifndef included_INSStaggeredHierarchyIntegrator
#define included_INSStaggeredHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscsys.h>

// IBAMR INCLUDES
#include <ibamr/INSHierarchyIntegrator.h>
#include <ibamr/INSStaggeredStokesOperator.h>

// IBTK INCLUDES
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonFACOperator.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/SCLaplaceOperator.h>
#include <ibtk/SCPoissonFACOperator.h>
#include <ibtk/SideDataSynchronization.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredHierarchyIntegrator provides a staggered-grid solver
 * for the incompressible Navier-Stokes equations on an AMR grid hierarchy.
 */
class INSStaggeredHierarchyIntegrator
    : public INSHierarchyIntegrator
{
public:
    /*!
     * The constructor for class INSStaggeredHierarchyIntegrator sets some
     * default values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    INSStaggeredHierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        bool register_for_restart=true);

    /*!
     * The destructor for class INSStaggeredHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~INSStaggeredHierarchyIntegrator();

    /*!
     * Set the problem coefficients used by the solver.
     */
    void
    setINSProblemCoefs(
        INSProblemCoefs problem_coefs);

    /*!
     * Get the problem coefficients used by the solver.
     */
    INSProblemCoefs
    getINSProblemCoefs() const;

    /*!
     * Register an operator to compute the convective acceleration term u*grad
     * u.
     *
     * By default, this solver uses class INSStaggeredPPMConvectiveOperator to
     * compute the convective acceleration term.
     *
     * If a NULL pointer is passed as the argument to this function, the
     * convective acceleration term is dropped, so that the solver effectively
     * solves the time-dependent (creeping) Stokes equations.
     */
    void
    setConvectiveOperator(
        SAMRAI::tbox::Pointer<IBTK::GeneralOperator> convective_op,
        ConvectiveDifferencingType convective_difference_form=UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE);

    /*!
     * Get the convective operator being used by this solver class.
     *
     * \note This operator is not initialized until the patch hierarchy is
     * initialized.
     */
    SAMRAI::tbox::Pointer<IBTK::GeneralOperator>
    getConvectiveOperator() const;

    /*!
     * Get the type of convective differencing being used by this solver class.
     */
    ConvectiveDifferencingType
    getConvectiveDifferenceForm() const;

    /*!
     * Register a preconditioner for the incompressible Stokes equations.
     *
     * By default, this solver uses class INSStaggeredProjectionPreconditioner
     * to precondition the incompressible Stokes system.
     */
    void
    setStokesPreconditioner(
        SAMRAI::tbox::Pointer<IBTK::LinearSolver> stokes_pc);

    /*!
     * Get the preconditioner being used by this solver class.
     *
     * \note This preconditioner is not initialized until the patch hierarchy is
     * initialized.
     */
    SAMRAI::tbox::Pointer<IBTK::LinearSolver>
    getStokesPreconditioner();

    /*!
     * Return a subdomain solver for the velocity subsystem.  Such solvers can
     * be useful in constructing block preconditioners.
     *
     * \note Data management for the preconditioner object is managed by this
     * solver class.
     */
    SAMRAI::tbox::Pointer<IBTK::LinearSolver>
    buildVelocitySubdomainSolver();

    /*!
     * Return a subdomain solver for the pressure subsystem.  Such solvers can
     * be useful in constructing block preconditioners.
     *
     * \note Data management for the preconditioner object is managed by this
     * solver class.
     */
    SAMRAI::tbox::Pointer<IBTK::LinearSolver>
    buildPressureSubdomainSolver();

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void
    initializeHierarchyIntegrator(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void
    preprocessIntegrateHierarchy(
        const double current_time,
        const double new_time,
        const int num_cycles=1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void
    integrateHierarchy(
        const double current_time,
        const double new_time,
        const int cycle_num=0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void
    postprocessIntegrateHierarchy(
        const double current_time,
        const double new_time,
        const bool skip_synchronize_new_state_data,
        const int num_cycles=1);

    /*!
     * Regrid the patch hierarchy.
     */
    void
    regridHierarchy();

protected:
    /*!
     * Return the maximum stable time step size.
     */
    double
    getTimeStepSizeSpecialized();

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     */
    void
    initializeLevelDataSpecialized(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
        const bool allocate_data);

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    void
    applyGradientDetectorSpecialized(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index,
        const bool initial_time,
        const bool uses_richardson_extrapolation_too);

    /*!
     * Prepare variables for plotting.
     */
    void
    setupPlotDataSpecialized();

    /*!
     * Write out specialized object state to the given database.
     */
    void
    putToDatabaseSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredHierarchyIntegrator(
        const INSStaggeredHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredHierarchyIntegrator&
    operator=(
        const INSStaggeredHierarchyIntegrator& that);

    /*!
     * Compute the appropriate source term that must be added to the momentum
     * equation when the fluid contains internal sources and sinks.
     */
    void
    computeDivSourceTerm(
        const int F_idx,
        const int Q_idx,
        const int U_idx);

    /*!
     * Project the velocity field following a regridding operation.
     */
    void
    regridProjection();

    /*!
     * Initialize the operators and solvers used by the hierarchy integrator.
     */
    void
    initializeOperatorsAndSolvers();

    /*!
     * Reinitialize the operators and solvers used by the hierarchy integrator.
     */
    void
    reinitializeOperatorsAndSolvers(
        const double current_time,
        const double new_time);

    /*!
     * Determine the largest stable timestep on an individual patch level.
     */
    double
    getStableTimestep(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level) const;

    /*!
     * Determine the largest stable timestep on an individual patch.
     */
    double
    getStableTimestep(
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

    /*!
     * Read input values from a given database.
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        const bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void
    getFromRestart();

    /*
     * Boolean value that indicates whether the integrator has been initialized.
     */
    bool d_integrator_is_initialized;

    /*
     * The fluid density (rho), dynamic viscosity (mu), kinematic viscosity
     * (nu), and (optional) drag coefficient (lambda).
     *
     * \note rho_water = 1.00 g cm^-3
     *       mu_water  = 0.01 g cm^-1 s^-1
     *       nu_water  = 0.01 cm^2 s^-1
     */
    INSProblemCoefs d_problem_coefs;

    /*
     * The maximum CFL number.
     *
     * NOTE: Empirical tests suggest that the maximum stable CFL number is
     * determined by the number of fixed-point iterations:
     *
     *    num_cycles == 1 ===> CFL <= 0.25
     *    num_cycles == 2 ===> CFL <= 0.5
     *    num_cycles == 3 ===> CFL <= 1.0
     */
    double d_cfl_max;

    /*
     * Cell tagging criteria based on the relative and absolute magnitudes of
     * the local vorticity.
     */
    bool d_using_vorticity_tagging;
    SAMRAI::tbox::Array<double> d_Omega_rel_thresh, d_Omega_abs_thresh;
    double d_Omega_max;

    /*
     * This boolean value determines whether the pressure is normalized to have
     * zero mean (i.e., discrete integral) at the end of each timestep.
     */
    bool d_normalize_pressure;

    /*
     * This boolean value determines whether the convective acceleration term is
     * included in the momentum equation.  (If it is not, this solver
     * effectively solves the so-called creeping Stokes equations.)
     */
    bool d_creeping_flow;

    /*
     * Threshold that determines whether the velocity field needs to be
     * reprojected following adaptive regridding.
     */
    double d_regrid_max_div_growth_factor;

    /*
     * Double precision values are (optional) factors used to rescale the
     * velocity, pressure, and force for plotting.
     *
     * Boolean values indicates whether to output various quantities for
     * visualization.
     */
    double d_U_scale, d_P_scale, d_F_scale, d_Q_scale, d_Omega_scale, d_Div_U_scale;
    bool d_output_U, d_output_P, d_output_F, d_output_Q, d_output_Omega, d_output_Div_U;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > d_hier_sc_data_ops;

    /*
     * Boundary condition and data synchronization operators.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_U_bdry_bc_fill_op, d_Q_bdry_bc_fill_op, d_no_fill_op;
    SAMRAI::tbox::Pointer<INSStaggeredPhysicalBoundaryHelper> d_U_bc_helper;
    blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM> d_U_bc_coefs, d_U_star_bc_coefs;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_Phi_bc_coef;
    SAMRAI::tbox::Pointer<IBTK::SideDataSynchronization> d_side_synch_op;

    /*
     * Hierarchy operators and solvers.
     */
    bool d_vectors_need_init;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_U_scratch_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_U_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_U_half_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_N_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_P_scratch_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_P_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_sol_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_nul_vec;

    bool d_stokes_op_needs_init;
    SAMRAI::tbox::Pointer<INSStaggeredStokesOperator> d_stokes_op;

    bool d_convective_op_needs_init;
    SAMRAI::tbox::Pointer<IBTK::GeneralOperator> d_convective_op;
    ConvectiveDifferencingType d_convective_difference_form;

    bool d_helmholtz_solver_needs_init;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>          d_helmholtz_hypre_pc_db, d_helmholtz_fac_pc_db;
    SAMRAI::tbox::Pointer<IBTK::SCLaplaceOperator>         d_helmholtz_op;
    SAMRAI::solv::PoissonSpecifications*                   d_helmholtz_spec;
    SAMRAI::tbox::Pointer<IBTK::SCPoissonHypreLevelSolver> d_helmholtz_hypre_pc;
    SAMRAI::tbox::Pointer<IBTK::SCPoissonFACOperator>      d_helmholtz_fac_op;
    SAMRAI::tbox::Pointer<IBTK::FACPreconditioner>         d_helmholtz_fac_pc;
    SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver>   d_helmholtz_solver;

    bool d_poisson_solver_needs_init;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>          d_poisson_hypre_pc_db, d_poisson_fac_pc_db;
    SAMRAI::tbox::Pointer<IBTK::CCLaplaceOperator>         d_poisson_op;
    SAMRAI::solv::PoissonSpecifications*                   d_poisson_spec;
    SAMRAI::tbox::Pointer<IBTK::CCPoissonHypreLevelSolver> d_poisson_hypre_pc;
    SAMRAI::tbox::Pointer<IBTK::CCPoissonFACOperator>      d_poisson_fac_op;
    SAMRAI::tbox::Pointer<IBTK::FACPreconditioner>         d_poisson_fac_pc;
    SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver>   d_poisson_solver;

    bool d_stokes_solver_needs_init;
    SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver>   d_stokes_solver;
    SAMRAI::tbox::Pointer<IBTK::LinearSolver>              d_stokes_pc;

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>          d_regrid_projection_fac_pc_db;

    /*!
     * State and scratch variables not defined in INSHierarchyIntegrator base
     * class.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_U_sc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_U_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_P_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_F_sc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_F_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_Q_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_Omega_var;
#if (NDIM == 3)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_Omega_Norm_var;
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_Div_U_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_U_regrid_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_U_src_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_indicator_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_F_div_var;

    /*
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int          d_U_current_idx,          d_U_new_idx,          d_U_scratch_idx;
    int       d_U_cc_current_idx,       d_U_cc_new_idx,       d_U_cc_scratch_idx;
    int          d_P_current_idx,          d_P_new_idx,          d_P_scratch_idx;
    int          d_F_current_idx,          d_F_new_idx,          d_F_scratch_idx;
    int       d_F_cc_current_idx,       d_F_cc_new_idx,       d_F_cc_scratch_idx;
    int          d_Q_current_idx,          d_Q_new_idx,          d_Q_scratch_idx;
    int      d_Omega_current_idx,      d_Omega_new_idx,      d_Omega_scratch_idx;
#if (NDIM == 3)
    int d_Omega_Norm_current_idx, d_Omega_Norm_new_idx, d_Omega_Norm_scratch_idx;
#endif
    int      d_Div_U_current_idx,      d_Div_U_new_idx,      d_Div_U_scratch_idx;

    /*
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_U_regrid_idx, d_U_src_idx, d_indicator_idx, d_F_div_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSStaggeredHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredHierarchyIntegrator
