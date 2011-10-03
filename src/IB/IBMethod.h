// Filename: IBMethod.h
// Created on 21 Sep 2011 by Boyce Griffith
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

#ifndef included_IBMethod
#define included_IBMethod

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/IBInstrumentPanel.h>
#include <ibamr/IBLagrangianForceStrategy.h>
#include <ibamr/IBLagrangianSourceStrategy.h>
#include <ibamr/IBMethodStrategy.h>
#include <ibamr/IBPostProcessStrategy.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBMethod is an implementation of the abstract base class
 * IBMethodStrategy that provides functionality required by the standard IB
 * method.
 */
class IBMethod
    : public IBMethodStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    IBMethod(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        bool register_for_restart=true);

    /*!
     * \brief Destructor.
     */
    ~IBMethod();

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
     * Return the number of ghost cells required by the Lagrangian-Eulerian
     * interaction routines.
     */
    const SAMRAI::hier::IntVector<NDIM>&
    getMinimumGhostCellWidth() const;

    /*!
     * Method to prepare to advance data from current_time to new_time.
     */
    void
    preprocessIntegrateData(
        double current_time,
        double new_time,
        int num_cycles);

    /*!
     * Method to clean up data following call(s) to integrateHierarchy().
     */
    void
    postprocessIntegrateData(
        double current_time,
        double new_time,
        int num_cycles);

    /*!
     * Interpolate the Eulerian velocity to the curvilinear mesh at the
     * specified time within the current time interval.
     */
    void
    interpolateVelocity(
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        double data_time);

    /*!
     * Advance the positions of the Lagrangian structure using forward Euler.
     */
    void
    eulerStep(
        double current_time,
        double new_time);

    /*!
     * Advance the positions of the Lagrangian structure using the midpoint
     * rule.
     */
    void
    midpointStep(
        double current_time,
        double new_time);

    /*!
     * Compute the Lagrangian force at the specified time within the current
     * time interval.
     */
    void
    computeLagrangianForce(
        double data_time);

    /*!
     * Spread the Lagrangian force to the Cartesian grid at the specified time
     * within the current time interval.
     */
    void
    spreadForce(
        int f_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& f_prolongation_scheds,
        double data_time);

    /*!
     * Indicate whether there are any internal fluid sources/sinks.
     */
    virtual bool
    hasFluidSources() const;

    /*!
     * Compute the Lagrangian source/sink density at the specified time within
     * the current time interval.
     *
     * An empty default implementation is provided.
     */
    void
    computeLagrangianFluidSource(
        double data_time);

    /*!
     * Spread the Lagrangian source/sink density to the Cartesian grid at the
     * specified time within the current time interval.
     *
     * An empty default implementation is provided.
     */
    void
    spreadFluidSource(
        int q_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& q_prolongation_scheds,
        double data_time);

    /*!
     * Compute the pressures at the positions of any distributed internal fluid
     * sources or sinks.
     *
     * An empty default implementation is provided.
     */
    void
    interpolatePressure(
        int p_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& p_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& p_ghost_fill_scheds,
        double data_time);

    /*!
     * Execute user-defined post-processing operations.
     */
    void
    postprocessData(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy);

    /*!
     * Initialize Lagrangian data corresponding to the given AMR patch hierarchy
     * at the start of a computation.  If the computation is begun from a
     * restart file, data may be read from the restart databases.
     *
     * A patch data descriptor is provided for the Eulerian velocity in case
     * initialization requires interpolating Eulerian data.  Ghost cells for
     * Eulerian data will be filled upon entry to this function.
     */
    void
    initializePatchHierarchy(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg,
        int u_data_idx,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >& u_synch_scheds,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& u_ghost_fill_scheds,
        int integrator_step,
        double init_data_time,
        bool initial_time);

    /*!
     * Register a load balancer and work load patch data index with the IB
     * strategy object.
     */
    void
    registerLoadBalancer(
        SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer,
        int workload_data_idx);

    /*!
     * Update work load estimates on each level of the patch hierarchy.
     */
    void
    updateWorkloadEstimates(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        int workload_data_idx);

    /*!
     * Begin redistributing Lagrangian data prior to regridding the patch
     * hierarchy.
     */
    void
    beginDataRedistribution(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Complete redistributing Lagrangian data following regridding the patch
     * hierarchy.
     */
    void
    endDataRedistribution(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
     */
    void
    initializeLevelData(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int level_number,
        double init_data_time,
        bool can_be_refined,
        bool initial_time,
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
        bool allocate_data);

    /*!
     * Reset cached hierarchy dependent data.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration
     */
    void
    resetHierarchyConfiguration(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int coarsest_level,
        int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to user-supplied feature detection criteria.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::applyGradientDetector
     */
    void
    applyGradientDetector(
        SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        int level_number,
        double error_data_time,
        int tag_index,
        bool initial_time,
        bool uses_richardson_extrapolation_too);

    /*!
     * Write out object state to the given database.
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log;

    /*
     * Pointers to the patch hierarchy and gridding algorithm objects associated
     * with this object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;

    /*
     * The current time step interval.
     */
    double d_current_time, d_new_time, d_half_time;

    /*
     * Boolean values tracking whether certain quantities need to be
     * reinitialized.
     */
    bool d_X_current_needs_ghost_fill, d_X_new_needs_ghost_fill, d_X_half_needs_ghost_fill;
    bool d_F_current_needs_ghost_fill, d_F_new_needs_ghost_fill, d_F_half_needs_ghost_fill;
    bool d_X_half_needs_reinit, d_U_half_needs_reinit;
    /*
     * The LDataManager is used to coordinate the distribution of Lagrangian
     * data on the patch hierarchy.
     */
    IBTK::LDataManager* d_l_data_manager;
    std::string d_interp_delta_fcn, d_spread_delta_fcn;
    SAMRAI::hier::IntVector<NDIM> d_ghosts;

    /*
     * Lagrangian variables.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_X_current_data, d_X_new_data, d_X_half_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_U_current_data, d_U_new_data, d_U_half_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > d_F_current_data, d_F_new_data, d_F_half_data;

    /*
     * List of local indices of local anchor points.
     *
     * NOTE: IB points are automatically considered to be anchored if they are
     * within 2.0*sqrt(epsilon_mach) of the physical boundary.
     */
    std::vector<std::set<int> > d_anchor_point_local_idxs;

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
    SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> d_ib_force_fcn;
    bool d_ib_force_fcn_needs_init;

    /*
     * The source/sink generators.
     */
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
     * Visualization data writers.
     */
    SAMRAI::tbox::Pointer<IBTK::LSiloDataWriter> d_silo_writer;
#if (NDIM == 3)
    SAMRAI::tbox::Pointer<IBTK::LM3DDataWriter> d_m3D_writer;
#endif

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBMethod();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBMethod(
        const IBMethod& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBMethod&
    operator=(
        const IBMethod& that);

    /*!
     * Reset the Lagrangian force function object.
     */
    void
    resetLagrangianForceFunction(
        double init_data_time,
        bool initial_time);

    /*!
     * Reset the Lagrangian source function object.
     */
    void
    resetLagrangianSourceFunction(
        double init_data_time,
        bool initial_time);

    /*!
     * Reset the post-processor object.
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
     * Interpolate the current and new data to obtain values at the midpoint of
     * the time interval.
     */
    void
    reinitMidpointData(
        const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& current_data,
        const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& new_data,
        const std::vector<SAMRAI::tbox::Pointer<IBTK::LData> >& half_data);

    /*!
     * Set the elements of the Lagrangian vector to zero at anchored nodes of
     * the curvilinear mesh.
     */
    void
    resetAnchorPointValues(
        std::vector<SAMRAI::tbox::Pointer<IBTK::LData> > U_data,
        int coarsest_ln,
        int finest_ln);

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
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.
     */
    std::string d_object_name;

    /*
     * A boolean value indicating whether the class is registered with the
     * restart database.
     */
    bool d_registered_for_restart;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBMethod.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBMethod
