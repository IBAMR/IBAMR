// Filename: IBHierarchyIntegrator.h
// Created on 12 Jul 2004 by Boyce Griffith
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

#ifndef included_IBHierarchyIntegrator
#define included_IBHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSC INCLUDES
#include <petscsys.h>

// IBAMR INCLUDES
#include <ibamr/INSHierarchyIntegrator.h>
#include <ibamr/IBEulerianForceFunction.h>
#include <ibamr/IBEulerianSourceFunction.h>
#include <ibamr/IBInstrumentPanel.h>
#include <ibamr/IBLagrangianForceStrategy.h>
#include <ibamr/IBLagrangianSourceStrategy.h>
#include <ibamr/IBPostProcessStrategy.h>

// IBTK INCLUDES
#include <ibtk/LDataManager.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBHierarchyIntegrator is an implementation of a formally
 * second-order accurate, semi-implicit version of the immersed boundary method.
 */
class IBHierarchyIntegrator
    : public HierarchyIntegrator
{
public:
    /*!
     * The constructor for class IBHierarchyIntegrator sets some default values,
     * reads in configuration information from input and restart databases, and
     * registers the integrator object with the restart manager when requested.
     */
    IBHierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<INSHierarchyIntegrator> ins_hier_integrator,
        bool register_for_restart=true);

    /*!
     * The destructor for class IBHierarchyIntegrator unregisters the integrator
     * object with the restart manager when the object is so registered.
     */
    ~IBHierarchyIntegrator();

    /*!
     * Supply a body force (optional).
     */
    void
    registerBodyForceFunction(
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> F_fcn);

    /*!
     * Supply a Lagrangian force object.
     */
    void
    registerIBLagrangianForceFunction(
        SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> ib_force_fcn);

    /*!
     * Supply a Lagrangian source object.
     */
    void
    registerIBLagrangianSourceFunction(
        SAMRAI::tbox::Pointer<IBLagrangianSourceStrategy> ib_source_fcn);

    /*!
     * Supply a Lagrangian initialization object.
     */
    void
    registerLInitStrategy(
        SAMRAI::tbox::Pointer<IBTK::LInitStrategy> l_initializer);

    /*!
     * Free references to Lagrangian initialization objects.
     */
    void
    freeLInitStrategy();

    /*!
     * Supply a postprocessor object.
     */
    void
    registerIBPostProcessor(
        SAMRAI::tbox::Pointer<IBPostProcessStrategy> ib_postprocessor);

    /*!
     * Return a pointer to the Lagrangian data manager object.
     */
    IBTK::LDataManager*
    getLDataManager() const;

    /*!
     * Return a pointer to the instrumentation manager object.
     */
    SAMRAI::tbox::Pointer<IBInstrumentPanel>
    getIBInstrumentPanel() const;

    /*!
     * Register a Lagrangian Silo data writer so this class will write plot
     * files that may be postprocessed with the VisIt visualization tool.
     */
    void
    registerLSiloDataWriter(
        SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> silo_writer);

#if (NDIM == 3)
    /*!
     * Register a Lagrangian myocardial3D data writer so this class will write
     * plot files that may be postprocessed with the myocardial3D visualization
     * program.
     */
    void
    registerLM3DDataWriter(
        SAMRAI::tbox::Pointer<IBTK::LM3DDataWriter> m3D_writer);
#endif

    /*!
     * Register a load balancer for non-uniform load balancing.
     */
    void
    registerLoadBalancer(
        SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer);

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
    void
    initializePatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Returns the number of cycles to perform for the present time step.
     */
    int
    getNumberOfCycles();

    /*!
     * Prepare to advance the data from current_time to new_time.
     */
    void
    preprocessIntegrateHierarchy(
        double current_time,
        double new_time,
        int num_cycles=1);

    /*!
     * Synchronously advance each level in the hierarchy over the given time
     * increment.
     */
    void
    integrateHierarchy(
        double current_time,
        double new_time,
        int cycle_num=0);

    /*!
     * Clean up data following call(s) to integrateHierarchy().
     */
    void
    postprocessIntegrateHierarchy(
        double current_time,
        double new_time,
        bool skip_synchronize_new_state_data,
        int num_cycles=1);

    /*!
     * Regrid the hierarchy.
     */
    void
    regridHierarchy();

    /*!
     * Execute user-defined post-processing operations.
     */
    void
    postProcessData();

protected:
    /*!
     * Function to determine whether to regridding should occur at the current
     * time step.
     */
    bool
    atRegridPointSpecialized() const;

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     */
    void
    initializeLevelDataSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
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
    resetHierarchyConfigurationSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int coarsest_level,
        int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    void
    applyGradientDetectorSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int level_number,
        double error_data_time,
        int tag_index,
        bool initial_time,
        bool uses_richardson_extrapolation_too);

    /*!
     * Write out specialized object state to the given database.
     */
    void
    putToDatabaseSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Time stepping scheme.
     */
    TimesteppingType d_timestepping_scheme;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBHierarchyIntegrator(
        const IBHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBHierarchyIntegrator&
    operator=(
        const IBHierarchyIntegrator& that);

    /*!
     * Initialize the IBLagrangianForceStrategy object for the current
     * configuration of the curvilinear mesh.
     */
    void
    resetLagrangianForceFunction(
        double init_data_time,
        bool initial_time);

    /*!
     * Initialize the IBLagrangianSourceStrategy object for the current
     * configuration of the curvilinear mesh.
     */
    void
    resetLagrangianSourceFunction(
        double init_data_time,
        bool initial_time);

    /*!
     * Initialize the IBPostProcessStrategy object for the current configuration
     * of the curvilinear mesh.
     */
    void
    resetPostProcessor(
        double init_data_time,
        bool initial_time);

    /*!
     * Compute the flow rates and pressures in the internal flow meters and
     * pressure gauges.
     */
    void
    updateIBInstrumentationData(
        int timestep_num,
        double data_time);

    /*!
     * Set the elements of the Lagrangian vector to zero at anchored nodes of
     * the curvilinear mesh.
     */
    void
    resetAnchorPointValues(
        std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > V_data,
        int coarsest_ln,
        int finest_ln);

    /*!
     * Perform a single forward Euler step for the pIB state variables.
     */
    void
    pIBEulerStep(
        SAMRAI::tbox::Pointer<IBTK::LData> K_data,
        SAMRAI::tbox::Pointer<IBTK::LData> M_data,
        SAMRAI::tbox::Pointer<IBTK::LData> X_current_data,
        SAMRAI::tbox::Pointer<IBTK::LData> Y_current_data,
        SAMRAI::tbox::Pointer<IBTK::LData> Y_new_data,
        SAMRAI::tbox::Pointer<IBTK::LData> dY_dt_current_data,
        SAMRAI::tbox::Pointer<IBTK::LData> dY_dt_new_data,
        double dt);

    /*!
     * Perform a single explicit midpoint-rule step for the pIB state variables.
     */
    void
    pIBMidpointStep(
        SAMRAI::tbox::Pointer<IBTK::LData> K_data,
        SAMRAI::tbox::Pointer<IBTK::LData> M_data,
        SAMRAI::tbox::Pointer<IBTK::LData> X_half_data,
        SAMRAI::tbox::Pointer<IBTK::LData> Y_current_data,
        SAMRAI::tbox::Pointer<IBTK::LData> Y_new_data,
        SAMRAI::tbox::Pointer<IBTK::LData> dY_dt_current_data,
        SAMRAI::tbox::Pointer<IBTK::LData> dY_dt_new_data,
        double dt);

    /*!
     * Compute pIB forces for a midpoint-rule type time stepping scheme.
     */
    void
    pIBMidpointForcing(
        SAMRAI::tbox::Pointer<IBTK::LData> F_data,
        SAMRAI::tbox::Pointer<IBTK::LData> K_data,
        SAMRAI::tbox::Pointer<IBTK::LData> X_half_data,
        SAMRAI::tbox::Pointer<IBTK::LData> Y_current_data,
        SAMRAI::tbox::Pointer<IBTK::LData> Y_new_data);

    /*!
     * Compute pIB forces for a trapezoidal-rule type time stepping scheme.
     */
    void
    pIBTRForcing(
        SAMRAI::tbox::Pointer<IBTK::LData> F_data,
        SAMRAI::tbox::Pointer<IBTK::LData> K_data,
        SAMRAI::tbox::Pointer<IBTK::LData> X_data,
        SAMRAI::tbox::Pointer<IBTK::LData> Y_data);

    /*!
     * Set the values of the distributed internal sources/sinks on the Cartesian
     * grid hierarchy.
     *
     * \note This method computes the point source/sink strengths, spreads those
     * values to the Cartesian grid using the cosine delta function, and
     * synchronizes the data on the hierarchy.
     */
    void
    computeSourceStrengths(
        int coarsest_level,
        int finest_level,
        double data_time,
        const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& X_data);

    /*!
     * Get the values of the pressures at the positions of the distributed
     * internal sources/sinks.
     *
     * \note This method interpolates the \em new Cartesian grid pressure at the
     * given locations of the internal sources/sinks.
     */
    void
    computeSourcePressures(
        int coarsest_level,
        int finest_level,
        double data_time,
        const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& X_data);

    /*!
     * Read input values from a given database.
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        bool is_from_restart);

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
     * The INSHierarchyIntegrator is used to provide time integration
     * capability for the incompressible Navier-Stokes equations.
     */
    SAMRAI::tbox::Pointer<INSHierarchyIntegrator> d_ins_hier_integrator;

    /*
     * The LDataManager is used to coordinate the distribution of Lagrangian
     * data on the patch hierarchy.
     */
    IBTK::LDataManager* d_l_data_manager;
    std::string d_interp_delta_fcn, d_spread_delta_fcn;
    SAMRAI::hier::IntVector<NDIM> d_ghosts;

    /*
     * Visualization data writers.
     */
    SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> d_silo_writer;
#if (NDIM == 3)
    SAMRAI::tbox::Pointer<IBTK::LM3DDataWriter> d_m3D_writer;
#endif

    /*
     * The load balancer.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;

    /*
     * Instrumentation (flow meter and pressure gauge) algorithms and data
     * structures.
     */
    SAMRAI::tbox::Pointer<IBInstrumentPanel> d_instrument_panel;
    std::vector<double> d_total_flow_volume;

    /*
     * The specification and initialization information for the Lagrangian data
     * used by the integrator.
     */
    SAMRAI::tbox::Pointer<IBTK::LInitStrategy> d_l_initializer;

    /*
     * The force generators.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_body_force_fcn;
    SAMRAI::tbox::Pointer<IBEulerianForceFunction> d_eulerian_force_fcn;
    SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> d_ib_force_fcn;
    bool d_ib_force_fcn_needs_init;

    /*
     * The source/sink generators.
     */
    SAMRAI::tbox::Pointer<IBEulerianSourceFunction> d_eulerian_source_fcn;
    SAMRAI::tbox::Pointer<IBLagrangianSourceStrategy> d_ib_source_fcn;
    bool d_ib_source_fcn_needs_init;
    std::vector<std::vector<blitz::TinyVector<double,NDIM> > > d_X_src;
    std::vector<std::vector<double > > d_r_src, d_P_src, d_Q_src;
    std::vector<int> d_n_src;
    bool d_normalize_source_strength;

    /*
     * Post-processor object.
     */
    SAMRAI::tbox::Pointer<IBPostProcessStrategy> d_post_processor;
    bool d_post_processor_needs_init;

    /*
     * Parameters for the penalty IB method for boundaries with additional
     * boundary mass.
     */
    bool d_using_pIB_method;
    blitz::TinyVector<double,NDIM> d_gravitational_acceleration;

    /*
     * The regrid CFL interval indicates the number of meshwidths a particle may
     * move in any coordinate direction between invocations of the regridding
     * process.
     *
     * NOTE: Currently, when the CFL-based regrid interval is specified, it is
     * used instead of the fixed regrid interval.
     */
    double d_regrid_cfl_interval, d_regrid_cfl_estimate;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM,double> > d_hier_velocity_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_pressure_cc_data_ops;

    /*
     * Eulerian variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_V_var, d_F_var, d_Q_var;
    int d_V_idx, d_F_idx, d_F_current_idx, d_Q_idx;

    /*
     * Lagrangian variables.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_X_current_data, d_X_new_data, d_X_half_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_U_current_data, d_U_new_data, d_U_half_data, d_U_old_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_F_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_K_data, d_M_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_Y_current_data, d_Y_new_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_dY_dt_current_data, d_dY_dt_new_data;

    /*
     * List of local indices of local anchor points.
     *
     * NOTE: IB points are automatically considered to be anchored if they are
     * within 2.0*sqrt(epsilon_mach) of the physical boundary.
     */
    std::vector<std::set<int> > d_anchor_point_local_idxs;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBHierarchyIntegrator
