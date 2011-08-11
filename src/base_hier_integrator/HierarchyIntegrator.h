// Filename: HierarchyIntegrator.h
// Created on 10 Aug 2011 by Boyce Griffith
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

#ifndef included_HierarchyIntegrator
#define included_HierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class HierarchyIntegrator provides an abstract interface for a time
 * integrator for a system of equations defined on an AMR grid hierarchy, along
 * with basic data management for variables defined on that hierarchy.
 */
class HierarchyIntegrator
    : public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>,
      public SAMRAI::tbox::Serializable
{
public:
    /*!
     * The constructor for class HierarchyIntegrator registers the integrator
     * object with the restart manager when requested.
     */
    HierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::hier::Database> input_db,
        bool register_for_restart=true);
    
    /*!
     * The destructor for class HierarchyIntegrator unregisters the integrator
     * object with the restart manager when the object is so registered.
     */
    ~HierarchyIntegrator();
    
    /*!
     * Return the name of the hierarchy integrator object.
     */
    const std::string&
    getName() const;

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by a concrete time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.
     *
     * \note This method must be implemented so that it is safe to make multiple
     * calls to this method.  Specifically, it must be possible for users to be
     * able to make an explicit call to initializeHierarchyIntegrator() prior to
     * calling initializePatchHierarchy().
     *
     * \note Because this method is called prior to the construction of the
     * patch hierarchy, this method cannot initialize patch data associated with
     * the variables managed by the integrator object, nor can it initialize the
     * communications schedules associated with the integrator.
     */
    virtual void
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
     * Integrate data on all patches on all levels of the patch hierarchy over
     * the specified time increment.
     */
    void
    advanceHierarchy(
        const double dt,
        const int num_cycles=1);

    /*!
     * Return the maximum stable time step size for the state data associated
     * with the current, scratch, or new variable context.
     *
     * A default implementation is provided that returns
     * min(dt_max,dt_growth_factor*dt_current).  The growth condition prevents
     * excessive changes in the time step size as the computation progresses.
     */
    virtual double
    getStableTimestep(
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx);

    /*!
     * Synchronize data defined on the grid hierarchy.
     *
     * A default implementation is provided that synchronizes the state data
     * associated with either the current or new data contexts, depending on
     * which VariableContext is passed as an argument.
     */
    virtual void
    synchronizeHierarchyData(
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx);

    /*!
     * Reset the current data to equal the new data, update the time level of
     * the current data, and deallocate the scratch and new data.
     *
     * A default implementation is provided that swaps data between the current
     * and new variable contexts prior to deallocating the scratch and new data
     * contexts.
     */
    virtual void
    resetTimeDependentHierarchyData(
        const double new_time);

    /*!
     * Reset the hierarchy integrator to the state at the beginning of the
     * current time step.
     *
     * A default implementation is provided that deallocates the scratch and new
     * data contexts.
     */
    virtual void
    resetIntegratorToPreadvanceState();

    /*!
     * Regrid the patch hierarchy.
     *
     * A default implementation is provided that calls
     * GriddingAlgorithm::regidAllFinerLevels() to regrid the patch hierarchy.
     */
    virtual void
    regridHierarchy();

    /*!
     * Return true if the current step count indicates that regridding should
     * occur.
     *
     * A default implementation is provided that indicates that the hierarchy
     * should be regridded at a fixed integer interval of time steps.
     */
    virtual bool
    atRegridPoint() const;

    /*!
     * Return the current integration time.
     */
    double
    getIntegratorTime() const;

    /*!
     * Return the initial integration time.
     */
    double
    getStartTime() const;

    /*!
     * Return the final integration time.
     */
    double
    getEndTime() const;

    /*!
     * Return the number of time steps taken by the integrator.
     */
    int
    getIntegratorStep() const;

    /*!
     * Return the maximum number of time steps allowed by the integrator.
     */
    int
    getMaxIntegratorSteps() const;

    /*!
     * Return true if any time steps remain, false otherwise.
     */
    bool
    stepsRemaining() const;

    /*!
     * Return a pointer to the patch hierarchy managed by the integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
    getPatchHierarchy() const;

    /*!
     * Return a pointer to the gridding algorithm object managed by the integrator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
    getGriddingAlgorithm() const;

    /*!
     * Register a VisIt data writer so the integrator can output data files that
     * may be postprocessed with the VisIt visualization tool.
     */
    void
    registerVisItDataWriter(
        SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer);

    /*!
     * Prepare variables for plotting.
     *
     * \note Concrete time integrator objects may require that this function be
     * called immediately before writing out VisIt visualization data.
     *
     * An empty default implementation is provided.
     */
    virtual void
    setupVisItPlotData();

    ///
    ///  Routines to implement the time integration scheme.
    ///

    /*!
     * Prepare to advance the data from current_time to new_time.
     *
     * An empty default implementation is provided.
     */
    virtual void
    preprocessIntegrateHierarchy(
        const double current_time,
        const double new_time,
        const int num_cycles=1);

    /*!
     * Advance the data from current_time to new_time but do not synchronize the
     * data on the hierarchy.
     */
    virtual void
    integrateHierarchy(
        const double current_time,
        const double new_time,
        const int cycle_num=0) = 0;

    /*!
     * Clean up data following call to integrateHierarchy().
     *
     * An empty default implementation is provided.
     */
    virtual void
    postprocessIntegrateHierarchy(
        const double current_time,
        const double new_time,
        const int num_cycles=1);

    ///
    ///  Implementations of functions declared in the
    ///  SAMRAI::mesh::StandardTagAndInitStrategy abstract base class.
    ///

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     *
     * A default implementation is provided that uses any registered
     * initialization functions to set initial values at the initial time, and
     * that copies data from the old hierarchy to the new one at subsequent
     * times.
     *
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::initializeLevelData
     */
    virtual void
    initializeLevelData(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level=SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> >(NULL),
        const bool allocate_data=true);

    /*!
     * Reset cached hierarchy dependent data.
     *
     * A default implementation is provided that sets up the default
     * communication algorithms following a regridding operation.
     * 
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration
     */
    virtual void
    resetHierarchyConfiguration(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to user-supplied feature detection criteria.
     *
     * An empty default implementation is provided.
     * 
     * \see SAMRAI::mesh::StandardTagAndInitStrategy::applyGradientDetector
     */
    virtual void
    applyGradientDetector(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index,
        const bool initial_time,
        const bool uses_richardson_extrapolation_too);

    ///
    ///  Routines to access to the variable contexts maintained by the
    ///  integrator.
    ///

    /*!
     * Return a pointer to the "current" variable context used by integrator.
     * Current data corresponds to state data at the beginning of a time step,
     * or when a new level is initialized.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
    getCurrentContext() const;

    /*!
     * Return a pointer to the "new" variable context used by integrator.  New
     * data corresponds to advanced state data at the end of a time step.  The
     * data is one time step later than the "current" data.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
    getNewContext() const;

    /*!
     * Return a pointer to the "scratch" variable context used by integrator.
     * Scratch data typically corresponds to storage that user-routines in the
     * concrete GodunovAdvector object manipulate; in particular, scratch data
     * contains ghost cells.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
    getScratchContext() const;

    ///
    ///  Implementations of functions declared in the SAMRAI::tbox::Serializable
    ///  abstract base class.
    ///

    /*!
     * Write out object state to the given database.  The implementation of this
     * method also calls putToDatabaseSpecialized().
     */
    void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    /*!
     * Write out specialized object state to the given database.
     *
     * An empty default implementation is provided.
     */
    virtual void
    putToDatabaseSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*!
     * Register a state variable with the integrator.  When a refine operator is
     * specified, the data for the variable are automatically maintained as the
     * patch hierarchy evolves.
     *
     * All state variables are registered with three contexts: current, new, and
     * scratch.  The current context of a state variable is maintained from time
     * step to time step and, if the necessary coarsen and refine operators are
     * specified, as the patch hierarchy evolves.
     */
    void
    registerVariable(
        int& current_idx,
        int& new_idx,
        int& scratch_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
        const SAMRAI::hier::IntVector<NDIM>& scratch_ghosts=SAMRAI::hier::IntVector<NDIM>(0),
        const std::string& coarsen_name="NO_COARSEN",
        const std::string& refine_name="NO_REFINE",
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> init_fcn=SAMRAI::tbox::Pointer<IBTK::CartGridFunction>(NULL));

    /*!
     * Register a "scratch" variable with the integrator.  This data is not
     * maintained from time step to time step.
     *
     * All scratch variables are registered with the scratch context.
     */
    void
    registerVariable(
        int& scratch_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
        const SAMRAI::hier::IntVector<NDIM>& ghosts=SAMRAI::hier::IntVector<NDIM>(0));

    /*!
     * Read input values from a given database.  The implementation of this
     * method also calls getFromInputSpecialized().
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        const bool is_from_restart);

    /*!
     * Read specialized object input values from a given database.
     */
    virtual void
    getFromInputSpecialized(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        const bool is_from_restart);
    
    /*!
     * Read object state from the restart file and initialize class data
     * members.  The implementation of this method also calls
     * getFromRestartSpecialized().
     */
    void
    getFromRestart();
    
    /*!
     * Read specialized object state from the restart file and initialize class
     * data members.
     */
    virtual void
    getFromRestartSpecialized();

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
    
    /*
     * Pointers to the patch hierarchy and gridding algorithm objects associated
     * with this time integration object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;

    /*
     * Indicates whether the hierarchy has been initialized.
     */
    bool d_hierarchy_is_initialized;

    /*
     * The object used to write out data for postprocessing by the VisIt
     * visualization tool.
     */
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > d_visit_writer;

    /*
     * Time and time step size data read from input or set at initialization.
     */
    double d_start_time, d_end_time;
    double d_dt_previous, d_dt_max, d_dt_growth_factor;

    /*
     * The maximum number of time steps that may be taken by the integrator
     * object.
     */
    int d_max_integrator_steps;

    /*
     * The number of integration steps taken between invocations of the
     * regridding process.
     */
    int d_regrid_interval;

    /*
     * The regrid mode.  "Standard" regridding involves only one call to
     * SAMRAI::mesh::GriddingAlgorithm::regridAllFinerLevels().  This limits the
     * amount that the AMR grid hierarchy can change within a single call to
     * regridHierarchy().  "Agressive" regridding involes multiple calls to
     * SAMRAI::mesh::GriddingAlgorithm::regridAllFinerLevels(), effectively
     * allowing arbitrary changes to the grid hierarchy configuration within a
     * single call to regridHierarchy().
     */    
    RegridMode d_regrid_mode;

    /*
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log;

    /*
     * The number of cells on each level by which tagged cells will be buffered
     * after they have selected for refinement.  The tag buffer helps to
     * guarantee that refined cells near important features in the solution will
     * remain refined until the level is regridded next.
     */
    SAMRAI::tbox::Array<int> d_tag_buffer;

    /*
     * Cached communications algorithms, strategies, and schedules.
     */
    typedef std::map<std::string,SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > > RefineAlgorithmMap;
    typedef std::map<std::string,SAMRAI::xfer::RefinePatchStrategy<NDIM>* > RefinePatchStrategyMap;
    typedef std::map<std::string,std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > > RefineScheduleMap;

    RefineAlgorithmMap     d_refine_algs;
    RefinePatchStrategyMap d_refine_strategies;
    RefineScheduleMap      d_refine_scheds;

    typedef std::map<std::string,SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > > CoarsenAlgorithmMap;
    typedef std::map<std::string,SAMRAI::xfer::CoarsenPatchStrategy<NDIM>* > CoarsenPatchStrategyMap;
    typedef std::map<std::string,std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > > CoarsenScheduleMap;
    
    CoarsenAlgorithmMap     d_calgs;
    CoarsenPatchStrategyMap d_cstrategies;
    CoarsenScheduleMap      d_cscheds;

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_fill_after_regrid;

    /*
     * SAMRAI::hier::Variable lists and SAMRAI::hier::ComponentSelector objects
     * are used for data management.
     */
    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_state_variables;
    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_scratch_variables;

    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_copy_scratch_to_current_fast;
    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_copy_scratch_to_current_slow;

    SAMRAI::hier::ComponentSelector d_current_data, d_new_data, d_scratch_data;

    /*!
     * Variable contexts.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_current_context, d_new_context, d_scratch_context;
    
private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    HierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    HierarchyIntegrator(
        const HierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    HierarchyIntegrator&
    operator=(
        const HierarchyIntegrator& that);
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/HierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_HierarchyIntegrator
