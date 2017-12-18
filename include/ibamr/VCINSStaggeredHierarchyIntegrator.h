// Filename: VCINSStaggeredHierarchyIntegrator.h
// Created on 25 Sep 2017 by Nishant Nangia
//
// Copyright (c) 2002-2014, Boyce Griffith
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

#ifndef included_IBAMR_VCINSStaggeredHierarchyIntegrator
#define included_IBAMR_VCINSStaggeredHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "CellVariable.h"
#include "EdgeVariable.h"
#include "HierarchyCellDataOpsReal.h"
#include "HierarchyEdgeDataOpsReal.h"
#include "HierarchyFaceDataOpsReal.h"
#include "HierarchyNodeDataOpsReal.h"
#include "HierarchySideDataOpsReal.h"
#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "NodeVariable.h"
#include "SAMRAIVectorReal.h"
#include "SideVariable.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/ibamr_enums.h"
#include "ibtk/SideDataSynchronization.h"
#include "ibtk/ibtk_enums.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

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
 * \brief Class VCINSStaggeredHierarchyIntegrator provides a staggered-grid solver
 * for the incompressible Navier-Stokes equations on an AMR grid hierarchy, with variable
 * coefficients.
 *
 * This class always uses a non-conservative discretization of the form of the momentum equation
 * \f$\rho(\frac{\partial u}{\partial t} + N(u)) = -\nabla p + \nabla \cdot \mu (\nabla u) + (\nabla u)^T )\f$ 
 * where \f$ N(u) = u \cdot \nabla u \f$ for convective_difference_form = ADVECTIVE and
 * \f$ N(u) = \nabla \cdot (u u) \f$ for convective_difference_form = CONSERVATIVE.
 * 
 * In other words, this class will NEVER treat the left-hand side of the momentum equation in conservative form
 * i.e. \frac{\partial \rho u}{\partial t} + \nabla \cdot (\rho u u)
 */
class VCINSStaggeredHierarchyIntegrator : public INSHierarchyIntegrator
{
public:
    /*!
     * The constructor for class VCINSStaggeredHierarchyIntegrator sets some
     * default values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     */
    VCINSStaggeredHierarchyIntegrator(const std::string& object_name,
                                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                      bool register_for_restart = true);

    /*!
     * The destructor for class VCINSStaggeredHierarchyIntegrator unregisters the
     * integrator object with the restart manager when the object is so
     * registered.
     */
    ~VCINSStaggeredHierarchyIntegrator();

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
    SAMRAI::tbox::Pointer<ConvectiveOperator> getConvectiveOperator();

    /*!
     * Get the subdomain solver for the velocity subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getVelocitySubdomainSolver();

    /*!
     * Get the subdomain solver for the pressure subsystem.  Such solvers can be
     * useful in constructing block preconditioners.
     */
    SAMRAI::tbox::Pointer<IBTK::PoissonSolver> getPressureSubdomainSolver();

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
    void initializeHierarchyIntegrator(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
                                       SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

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
                                  SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void preprocessIntegrateHierarchy(double current_time, double new_time, int num_cycles = 1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void integrateHierarchy(double current_time, double new_time, int cycle_num = 0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void postprocessIntegrateHierarchy(double current_time,
                                       double new_time,
                                       bool skip_synchronize_new_state_data,
                                       int num_cycles = 1);

    /*!
     * Regrid the patch hierarchy.
     */
    void regridHierarchy();

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
     * Copy the solution data into the state data maintained by
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

    /*!
     * Set the interpolation type used for material properties rho and mu
     */
    void setVCInterpolationType(const IBTK::VCInterpType vc_interp_type);

    /*!
     * \brief Function to reset fluid density or viscosity if they are
     * maintained by this integrator.
     */
    typedef void (*ResetFluidPropertiesFcnPtr)(int property_idx,
                                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                               int cycle_num,
                                               double time,
                                               double new_time,
                                               void* ctx);

    /*!
     * \brief Register interface neighborhood locating functions.
     */
    void registerResetFluidDensityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Register interface neighborhood locating functions.
     */
    void registerResetFluidViscosityFcn(ResetFluidPropertiesFcnPtr callback, void* ctx);

    /*!
     * \brief Get whether or not density is constant
     */
    inline bool rhoIsConstant() const
    {
        return d_rho_is_const;
    }

    /*!
     * \brief Get whether or not viscosity is constant
     */
    inline bool muIsConstant() const
    {
        return d_mu_is_const;
    }

    /*!
     * \brief Get cell centered viscosity patch data index
     */
    inline int getMuPatchDataIndex() const
    {
        return d_mu_scratch_idx;
    }

    /*!
     * \brief Get interpolated node (2D) or edge (3D) centered viscosity patch data index
     */
    inline int getInterpolatedMuPatchDataIndex() const
    {
        return d_mu_interp_idx;
    }

    /*!
     * \brief Get the interpolated side centered density patch data index
     */
    inline int getInterpolatedRhoPatchDataIndex() const
    {
        return d_rho_interp_idx;
    }

    /*!
     * \brief Get the scaling factor used for A, p and u_rhs
     */
    inline double getPressureScalingFactor() const
    {
        return d_A_scale;
    }

protected:
    /*!
     * Determine the largest stable timestep on an individual patch.
     */
    double getStableTimestep(SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch) const;

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
                                        bool allocate_data);

    /*!
     * Reset cached hierarchy dependent data.
     */
    void
    resetHierarchyConfigurationSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                           int coarsest_level,
                                           int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    void applyGradientDetectorSpecialized(SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
                                          int level_number,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too);

    /*!
     * Prepare variables for plotting.
     */
    void setupPlotDataSpecialized();

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    VCINSStaggeredHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    VCINSStaggeredHierarchyIntegrator(const VCINSStaggeredHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    VCINSStaggeredHierarchyIntegrator& operator=(const VCINSStaggeredHierarchyIntegrator& that);

    /*!
     * Compute the appropriate source term that must be added to the momentum
     * equation when the fluid contains internal sources and sinks.
     */
    void computeDivSourceTerm(int F_idx, int Q_idx, int U_idx);

    /*!
     * Project the velocity field following a regridding operation.
     */
    void regridProjection();

    /*!
     * Determine the convective time stepping type for the current time step and
     * cycle number.
     */
    TimeSteppingType getConvectiveTimeSteppingType(int cycle_num);

    /*!
     * Preprocess the operators and solvers used by the hierarchy integrator.
     */
    void preprocessOperatorsAndSolvers(double current_time, double new_time);

    /*!
     * Update the operators and solvers to account for changes due to time-dependent coefficients
     */
    void updateOperatorsAndSolvers(double current_time, double new_time);

    /*!
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM, double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyFaceDataOpsReal<NDIM, double> > d_hier_fc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM, double> > d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyNodeDataOpsReal<NDIM, double> > d_hier_nc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyEdgeDataOpsReal<NDIM, double> > d_hier_ec_data_ops;

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

    std::string d_stokes_solver_type, d_stokes_precond_type;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_stokes_solver_db, d_stokes_precond_db;
    SAMRAI::tbox::Pointer<StaggeredStokesSolver> d_stokes_solver;
    bool d_stokes_solver_needs_init;

    /*!
     * Fluid solver variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_U_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_P_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_F_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_F_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Q_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_N_old_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Omega_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Div_U_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_Omega_Norm_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_regrid_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_U_src_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_indicator_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_F_div_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_EE_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_rho_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_mu_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_pressure_D_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_pressure_rhs_D_var;
#if (NDIM == 2)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_velocity_D_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_velocity_rhs_D_var;
#elif (NDIM == 3)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > d_velocity_D_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > d_velocity_rhs_D_var;
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_velocity_D_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_velocity_C_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_velocity_rhs_C_var;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_N_full_var;

    /*!
     * Interpolated material property variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_rho_interp_var;
#if (NDIM == 2)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::NodeVariable<NDIM, double> > d_mu_interp_var;
#elif (NDIM == 3)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::EdgeVariable<NDIM, double> > d_mu_interp_var;
#endif

    /*!
     * Functions resetting rho and mu if they are maintained by this integrator.
     */
    std::vector<ResetFluidPropertiesFcnPtr> d_reset_rho_fcns, d_reset_mu_fcns;
    std::vector<void *> d_reset_rho_fcns_ctx, d_reset_mu_fcns_ctx;

    /*!
     * Temporary storage variables that contain intermediate quantities
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, double> > d_temp_sc_var;
    int d_temp_sc_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, double> > d_temp_cc_var;
    int d_temp_cc_idx;

    /*
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int d_U_current_idx, d_U_new_idx, d_U_scratch_idx;
    int d_P_current_idx, d_P_new_idx, d_P_scratch_idx;
    int d_F_current_idx, d_F_new_idx, d_F_scratch_idx;
    int d_Q_current_idx, d_Q_new_idx, d_Q_scratch_idx;
    int d_N_old_current_idx, d_N_old_new_idx, d_N_old_scratch_idx;
    int d_rho_current_idx, d_rho_new_idx, d_rho_scratch_idx;
    int d_mu_current_idx, d_mu_new_idx, d_mu_scratch_idx;

    /*
     * Patch data descriptor indices for all "plot" variables managed by the
     * integrator.
     *
     * Plot variables have one context: current.
     */
    int d_U_cc_idx, d_F_cc_idx, d_Omega_idx, d_Div_U_idx, d_EE_idx;

    /*
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_Omega_Norm_idx, d_U_regrid_idx, d_U_src_idx, d_indicator_idx, d_F_div_idx;
    int d_velocity_C_idx, d_velocity_D_idx, d_velocity_D_cc_idx, d_pressure_D_idx;
    int d_velocity_rhs_C_idx, d_velocity_rhs_D_idx, d_pressure_rhs_D_idx;
    int d_rho_interp_idx, d_mu_interp_idx;
    int d_N_full_idx;

    /*
     * Variables to indicate if either rho or mu is constant
     */
    bool d_rho_is_const, d_mu_is_const;

    /*
     * Variable to indicate the type of interpolation to be done for rho and mu
     */
    IBTK::VCInterpType d_vc_interp_type;

    /*
     * Variable to indicate the scaling factor used for A, p and u_rhs
     */
    double d_A_scale;

    /*
     * Variable to set how often the preconditioner is reinitialized
     */
    int d_precond_reinit_interval;
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBAMR_VCINSStaggeredHierarchyIntegrator
