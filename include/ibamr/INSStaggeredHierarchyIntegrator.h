// ---------------------------------------------------------------------
//
// Copyright (c) 2008 - 2023 by the IBAMR developers
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

#ifndef included_IBAMR_INSStaggeredHierarchyIntegrator
#define included_IBAMR_INSStaggeredHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/config.h>

#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/StaggeredStokesSolverManager.h"
#include "ibamr/ibamr_enums.h"

#include "ibtk/SideDataSynchronization.h"

#include "CellVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <string>
#include <vector>

namespace IBAMR
{
class ConvectiveOperator;
} // namespace IBAMR
namespace IBTK
{
class PoissonSolver;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BasePatchLevel;
template <int DIM>
class Patch;
template <int DIM>
class PatchHierarchy;
template <int DIM>
class BasePatchHierarchy;
} // namespace hier
namespace mesh
{
template <int DIM>
class GriddingAlgorithm;
} // namespace mesh
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredHierarchyIntegrator provides a staggered-grid solver
 * for the incompressible Navier-Stokes equations on an AMR grid hierarchy.
 */
class INSStaggeredHierarchyIntegrator : public INSHierarchyIntegrator
{
public:
    /*!
     * The constructor for class INSStaggeredHierarchyIntegrator sets some
     * default values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    INSStaggeredHierarchyIntegrator(std::string object_name,
                                    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                    bool register_for_restart = true);

    /*!
     * The destructor for class INSStaggeredHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~INSStaggeredHierarchyIntegrator();

    /*!
     * Get the convective operator being used by this solver class.
     *
     * If the time integrator is configured to solve the time-dependent
     * (creeping) Stokes equations, then the returned pointer will be NULL.
     *
     * If the convective operator has not already been constructed, and if the
     * time integrator is not configured to solve the time-dependent (creeping)
     * Stokes equations, then this function will initialize the default type of
     * convective operator, which may be set in the class input database.
     */
    SAMRAI::tbox::Pointer<ConvectiveOperator> getConvectiveOperator() override;

    /*!
     * Get the subdomain solver for the velocity subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getVelocitySubdomainSolver() override;

    /*!
     * Get the subdomain solver for the pressure subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getPressureSubdomainSolver() override;

    /*!
     * Register a solver for the time-dependent incompressible Stokes equations.
     */
    void setStokesSolver(SAMRAI::tbox::Pointer<StaggeredStokesSolver> stokes_solver);

    /*!
     * Get the solver for the time-dependent incompressible Stokes equations
     * used by this solver class.
     */
    SAMRAI::tbox::Pointer<StaggeredStokesSolver> getStokesSolver();

    /*!
     * Indicate that the Stokes solver should be (re-)initialized before the
     * next time step.
     */
    void setStokesSolverNeedsInit();

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
    initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Initialize the AMR patch hierarchy and data defined on the hierarchy at
     * the start of a computation.  If the computation is begun from a restart
     * file, the patch hierarchy and patch data are read from the hierarchy
     * database.  Otherwise, the patch hierarchy and patch data are initialized
     * by the gridding algorithm associated with the integrator object.
     *
     * The implementation of this function assumes that the hierarchy exists
     * upon entry to the function, but that it contains no patch levels.  On
     * return from this function, the state of the integrator object will be
     * such that it is possible to step through time via the advanceHierarchy()
     * function.
     */
    void initializePatchHierarchy(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg) override;

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1) override;

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1) override;

    /*!
     * Setup solution and RHS vectors using state data maintained by the
     * integrator.
     */
    void setupSolverVectors(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol_vec,
                            const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& rhs_vec,
                            double current_time,
                            double new_time,
                            int cycle_num);

    /*!
     * Reset solution and RHS vectors using state data maintained by the
     * integrator, and copy the solution data into the state data maintained by
     * the integrator.
     */
    void resetSolverVectors(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol_vec,
                            const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& rhs_vec,
                            double current_time,
                            double new_time,
                            int cycle_num);

    /*!
     * Explicitly remove nullspace components from a solution vector.
     */
    void removeNullSpace(const SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> >& sol_vec);

protected:
    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchySpecialized(double current_time, double new_time, int cycle_num = 0) override;

    /*!
     * L1 norm of the discrete divergence of the fluid velocity before regridding.
     */
    double d_div_U_norm_1_pre = 0.0;

    /*!
     * L2 norm of the discrete divergence of the fluid velocity before regridding.
     */
    double d_div_U_norm_2_pre = 0.0;

    /*!
     * L-infinity norm of the discrete divergence of the fluid velocity before regridding.
     */
    double d_div_U_norm_oo_pre = 0.0;

    /*!
     * L1 norm of the discrete divergence of the fluid velocity after regridding.
     */
    double d_div_U_norm_1_post = 0.0;

    /*!
     * L2 norm of the discrete divergence of the fluid velocity after regridding.
     */
    double d_div_U_norm_2_post = 0.0;

    /*!
     * L-infinity norm of the discrete divergence of the fluid velocity after regridding.
     */
    double d_div_U_norm_oo_post = 0.0;

    /*!
     * Whether we need to perform a regrid projection when (re-)initializing composite hierarchy data.
     */
    bool d_do_regrid_projection = false;

    /*!
     * Determine the largest stable timestep on an individual patch.
     */
    double getStableTimestep(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const override;

    /*!
     * Write out specialized object state to the given database.
     */
    void putToDatabaseSpecialized(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db) override;

    /*!
     * Prepare the current hierarchy for regridding. Here we calculate the divergence.
     */
    void regridHierarchyBeginSpecialized() override;

    /*!
     * Update the current hierarchy data after regridding. Here we recalculate
     * the divergence and, if it has grown by a factor more than
     * d_regrid_max_div_growth_factor, we indicate that a regrid projection is needed
     * to (re)initialize the composite hierarchy data.
     */
    void regridHierarchyEndSpecialized() override;

    /*!
     * Perform data initialization after the entire hierarchy has been constructed.
     */
    void initializeCompositeHierarchyDataSpecialized(double init_data_time, bool initial_time) override;

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     */
    void initializeLevelDataSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                        int level_number,
                                        double init_data_time,
                                        bool can_be_refined,
                                        bool initial_time,
                                        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
                                        bool allocate_data) override;

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level) override;

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    void applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                          int level_number,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too) override;

    /*!
     * Prepare variables for plotting.
     */
    void setupPlotDataSpecialized() override;

    /*!
     * Project the velocity field following a regridding operation.
     */
    void regridProjection() override;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredHierarchyIntegrator() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredHierarchyIntegrator(const INSStaggeredHierarchyIntegrator& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredHierarchyIntegrator& operator=(const INSStaggeredHierarchyIntegrator& that) = delete;

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();

    /*!
     * Compute the appropriate source term that must be added to the momentum
     * equation when the fluid contains internal sources and sinks.
     */
    void computeDivSourceTerm(int F_idx, int Q_idx, int U_idx);

    /*!
     * Reinitialize the operators and solvers used by the hierarchy integrator.
     */
    void reinitializeOperatorsAndSolvers(double current_time, double new_time);

    /*!
     * Determine the convective time stepping type for the current time step and
     * cycle number.
     */
    TimeSteppingType getConvectiveTimeSteppingType(int cycle_num);

    /*!
     * Determine the time step size ratio.
     */
    double getTimeStepSizeRatio() const;

    /*!
     * Apply divergence-preserving (and curl-preserving) prolongation to the specified staggered-grid velocity field.
     */
    void applyDivergencePreservingProlongation(int U_idx,
                                               int level_number,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > old_level,
                                               SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                               double init_data_time) const;

    /*!
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM, double> > d_hier_fc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> > d_hier_sc_data_ops;

    /*
     * Boundary condition and data synchronization operators.
     */
    SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> d_bc_helper;
    SAMRAI::tbox::Pointer<IBTK::SideDataSynchronization> d_side_synch_op;

    /*
     * Hierarchy operators and solvers.
     */
    int d_coarsest_reset_ln, d_finest_reset_ln;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_U_scratch_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_U_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_U_adv_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_N_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_P_scratch_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_P_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_sol_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_rhs_vec;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_nul_vecs;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > > d_U_nul_vecs;
    bool d_vectors_need_init, d_explicitly_remove_nullspace;

    std::string d_stokes_solver_type = StaggeredStokesSolverManager::UNDEFINED,
                d_stokes_precond_type = StaggeredStokesSolverManager::UNDEFINED,
                d_stokes_sub_precond_type = StaggeredStokesSolverManager::UNDEFINED;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_stokes_solver_db, d_stokes_precond_db, d_stokes_sub_precond_db;
    SAMRAI::tbox::Pointer<StaggeredStokesSolver> d_stokes_solver;
    bool d_stokes_solver_needs_init;

    /*!
     * Fluid solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_P_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_F_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_F_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Q_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_N_old_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_old_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Omega_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_Omega_nc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Div_U_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_regrid_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_src_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_indicator_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_F_div_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_EE_var;

    std::string d_N_coarsen_type = "CONSERVATIVE_COARSEN";
    std::string d_N_refine_type = "CONSERVATIVE_LINEAR_REFINE";

    std::string d_U_P_bdry_interp_type = "LINEAR";

    /*!
     * Variables for graphical output.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_U_nc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_P_nc_var;

    /*
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_U_current_idx = IBTK::invalid_index, d_U_new_idx = IBTK::invalid_index, d_U_scratch_idx = IBTK::invalid_index;
    int d_P_current_idx = IBTK::invalid_index, d_P_new_idx = IBTK::invalid_index, d_P_scratch_idx = IBTK::invalid_index;
    int d_F_current_idx = IBTK::invalid_index, d_F_new_idx = IBTK::invalid_index, d_F_scratch_idx = IBTK::invalid_index;
    int d_Q_current_idx = IBTK::invalid_index, d_Q_new_idx = IBTK::invalid_index, d_Q_scratch_idx = IBTK::invalid_index;
    int d_N_old_current_idx = IBTK::invalid_index, d_N_old_new_idx = IBTK::invalid_index,
        d_N_old_scratch_idx = IBTK::invalid_index;
    int d_U_old_current_idx = IBTK::invalid_index, d_U_old_new_idx = IBTK::invalid_index,
        d_U_old_scratch_idx = IBTK::invalid_index;

    /*
     * Patch data descriptor for state variables which are only present in the current context.
     */
    int d_Div_U_idx = IBTK::invalid_index, d_Omega_idx = IBTK::invalid_index;

    /*
     * Patch data descriptor indices for all "plot" variables managed by the
     * integrator. These are only used for directly computing graphical output.
     * In particular, this list does *not* contain some scratch variables used
     * to generate graphical output. d_Div_U_idx is also directly plotted.
     *
     * Plot variables have one context: plot.
     */
    int d_U_nc_idx = IBTK::invalid_index, d_P_nc_idx = IBTK::invalid_index, d_F_cc_idx = IBTK::invalid_index,
        d_Omega_nc_idx = IBTK::invalid_index, d_EE_idx = IBTK::invalid_index;

    /**
     * Vector of all such plot-only data indices.
     */
    std::vector<int> d_plot_indices;

    /*
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_U_regrid_idx = IBTK::invalid_index, d_U_src_idx = IBTK::invalid_index, d_indicator_idx = IBTK::invalid_index,
        d_F_div_idx = IBTK::invalid_index;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBAMR_INSStaggeredHierarchyIntegrator
