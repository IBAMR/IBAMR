#ifndef included_ParabolicLevelIntegrator
#define included_ParabolicLevelIntegrator

// Filename: ParabolicLevelIntegrator.h
// Last modified: <16.Apr.2007 05:54:42 boyce@bigboy.nyconnect.com>
// Created on 09 Jan 2007 by Boyce Griffith (boyce@box221.cims.nyu.edu)

// NOTE: This implementation is directly derived from the SAMRAI
// HyperbolicLevelIntegrator source code.

/////////////////////////////// INCLUDES /////////////////////////////////////

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * XXXX
 */
class ParabolicLevelIntegrator
    : public SAMRAI::algs::TimeRefinementLevelStrategy<NDIM>,
      public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>,
      public tbox::Serializable
{
public:
    /*!
     * Enumerated type for the different ways in which variable
     * storage can be manipulated by the level integration algorithm.
     * See registerVariable(...) function for more details.
     *

     * - \b TIME_DEP      {Data that changes in time and needs more than one
     *                       time level to be stored.}
     * - \b INPUT         {Data that is set once and do not change during
     *                      the ghosts are never re-filled outside of
     *                      user-defined routines.}
     * - \b FLUX          {Face-centered double values used in conservative
     *                      difference and synchronization (i.e., refluxing)
     *                      process.  A corresponding variable to store flux
     *                      integral information is created for each FLUX
     *                      variable.}
     * - \b TEMPORARY     {Accessory values intended to live only for
     *                      computation on a single patch (i.e., they cannot
     *                      be assumed to exist between patch routine function
     *                      calls.)}
     *

     */
    enum PAR_VAR_TYPE{ TIME_DEP = 0,
                       INPUT = 1,
                       NO_FILL = 2,
                       FLUX = 3,
                       TEMPORARY = 4 };

    /*!
     * Constructor for ParabolicLevelIntegrator initializes
     * integration parameters to default values and constructs
     * standard communication algorithms.  Other data members are read
     * in from the specified input database or the restart database
     * corresponding to the specified object_name.  The constructor
     * also registers this object for restart using the specified
     * object name when the boolean argument is true.  Whether object
     * will write its state to restart files during program execution
     * is determined by this argument.  Note that it has a default
     * state of true.  This class is used by the time refinement
     * integrator for refined timestepping when the
     * using_time_refinement argument is true, and for synchronized
     * timestepping when the boolean is false.
     *
     * When assertion checking is active, passing in any null pointer
     * or an empty string will result in an unrecoverable exception.
     */
    ParabolicLevelIntegrator(
        const std::string& object_name,
        tbox::Pointer<tbox::Database> input_db,
        ParabolicPatchStrategy<NDIM>* patch_strategy,
        bool register_for_restart = true,
        bool using_time_refinement = true);

    /*!
     * The destructor for ParabolicLevelIntegrator unregisters the
     * integrator object with the restart manager when so registered.
     */
    virtual
    ~ParabolicLevelIntegrator();

    /*!
     * Initialize level integrator by by setting the number of time
     * levels of data needed based on specifications of the gridding
     * algorithm.
     *
     * This routine also invokes variable registration in the patch
     * strategy.
     *
     * Assertion checking will throw unrecoverable exceptions if
     * either pointer is null.
     */
    virtual void
    initializeLevelIntegrator(
        tbox::Pointer<mesh::BaseGriddingAlgorithm<NDIM> > base_gridding_alg);

    /*!
     * Determine time increment to advance data on level and return
     * that value.  The double dt_time argument is the simulation time
     * when the routine is called.  The initial_time boolean is true
     * if this routine is called during hierarchy initialization
     * (i.e., at the initial simulation time).  Otherwise, it is
     * false.  The recompute_dt option specifies whether to compute
     * the timestep using the current level data or to return the
     * value stored by the time integrator. The default true setting
     * means the timestep will be computed if no value is supplied.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the level pointer is null.
     */
    virtual double
    getLevelDt(
        const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
        const double dt_time,
        const bool initial_time);

    /*!
     * Return the maximum allowable time increment for the level with
     * the specified level number based on the time increment for the
     * next coarser level and the mesh refinement ratio between the
     * two levels.  For the common explicit integration methods for
     * hyperbolic conservation laws (constrained by a CFL limit), the
     * fine time increment is typically the coarse increment divided
     * by the refinement ratio.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the ratio vector is not acceptable (i.e., all
     * values > 0).
     */
    virtual double
    getMaxFinerLevelDt(
        const int finer_level_number,
        const double coarse_dt,
        const hier::IntVector<NDIM>& ratio_to_coarser);

    /*!
     * Integrate data on all patches on the given patch level from
     * current time (current_time) to new time (new_time).  This
     * routine is used to advance the solution on each level in the
     * hierarchy and during time-dependent regridding procedures, such
     * as Richardson extrapolation.  The boolean arguments are used to
     * determine the state of the algorithm and the data when the
     * advance routine is called.  The first_step and last_step
     * indicate whether the step is the first or last in the current
     * timestep sequence on the level.  Typically, the current
     * timestep sequence means each step on the level between advance
     * steps on a coarser level in the hierarchy, if one exists.  The
     * regrid_advance value is true when the advance is called as part
     * of a time-dependent regridding procedure.  Usually when this
     * happens, the results of the colution advance will be discarded.
     * So, for example, when this is true flux information is not
     * maintained and flux integrals are not updated.  The final
     * boolean argument indicates whether or not the level resides in
     * the hierarchy.  For example, during time-dependent regridding,
     * such as Richardson extrapolation, a temporary level that is not
     * in the hierarchy is created and advanced.  Then, a
     * communication schedule must be generated for the level before
     * the advance begins.
     *
     * This routine is called at two different points during time
     * integration: during the regular time advance sequence, and
     * possibly at the initial simulation time.  The second call
     * advances the solution on a coarser level ahead in time to
     * provide time-dependent boundary values for some finer level
     * when time-dependent regridding is used.  In the first case, the
     * values of the boolean flags are:
     *

     *    - \b  first_step
     *        = true for first step in level time step sequence; else, false.
     *    - \b  last_step
     *        = true for last step in level time step sequence; else, false.
     *    - \b  regrid_advance
     *        = false.
     *

     *
     * In the second case, the values of the boolean flags are:
     *

     *    - \b  first_step
     *        = true.
     *    - \b  last_step
     *        = false.
     *    - \b  regrid_advance
     *        = true.
     *

     *
     * When time-dependent regridding (i.e., Richardson extrapolation)
     * is used, the routine is called from two different points in
     * addition to those described above: to advance a temporary level
     * that is coarser than the hierarchy level on which error
     * estimation is performed, and to advance the hierarchy level
     * itself.  In the first case, the values of the boolean flags
     * are:
     *


     *    - \b  first_step
     *        = true.
     *    - \b  last_step
     *        = true.
     *    - \b  regrid_advance
     *        = true.
     *


     *
     * In the second case, the values of the boolean flags are:
     *


     *    - \b  first_step
     *      (when regridding during time integration sequence)
     *        = true when the level is not coarsest level to synchronize
     *          immediately before the regridding process; else, false.
     *      (when generating initial hierarchy construction)
     *        = true, even though there may be multiple advance steps.
     *    - \b  last_step
     *        = true when the advance is the last in the Richardson
     *          extrapolation step sequence; else false.
     *    - \b  regrid_advance
     *        = true.
     *


     *
     * When assertion checking is active, an unrecoverable exception
     * will result if either the level or hierarchy pointer is null,
     * or the new time is not greater than the given time.
     */
    virtual double
    advanceLevel(
        const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
        const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
        const double current_time,
        const double new_time,
        const bool first_step,
        const bool last_step,
        const bool regrid_advance = false);

    /*!
     * Synchronize data between given patch levels in patch hierarchy
     * according to the standard hyperbolic AMR flux correction
     * algorithm.  This routine synchronizes data between two levels
     * at a time from the level with index finest_level down to the
     * level with index coarsest_level.  The array of old time values
     * are used in the re-integration of the time-dependent data.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the hierarchy pointer is null, the level numbers
     * do not properly match existing levels in the hierarchy (either
     * coarsest_level > finest_level or some level is null), or all of
     * the old time values are less than the value of sync_time.
     */
    virtual void
    standardLevelSynchronization(
        const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const tbox::Array<double>& old_times);

    /*!
     * This overloaded version of standardLevelSynchronization
     * implements a routine used for synchronized timestepping.  Only
     * a single value for the old time is needed, since all levels
     * would have the same old time.
     */
    virtual void
    standardLevelSynchronization(
        const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const double old_time);

    /*!
     * Coarsen current solution data from finest hierarchy level
     * specified down through the coarsest hierarchy level specified,
     * if initial_time is true.  In this case, the hierarchy is being
     * constructed at the initial simulation time, After data is
     * coarsened, the application- specific initialization routine is
     * called to set data before that solution is further coarsened to
     * the next coarser level in the hierarchy.  This operation makes
     * the solution consistent between coarser levels and finer levels
     * that did not exist when the coarse levels where created and
     * initialized originally.
     *
     * When initial_time is false, this routine does nothing since the
     * standard hyperbolic AMR algorithm for conservation laws
     * requires no data synchronization after regridding beyond
     * interpolation of data from coarser levels in the hierarchy in
     * some conservative fashion.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the hierarchy pointer is null, the level numbers
     * do not properly match existing levels in the hierarchy (either
     * coarsest_level > finest_level or some level is null).
     */
    virtual void
    synchronizeNewLevels(
        const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const bool initial_time);

    /*!
     * Resets time-dependent data storage and update time for patch
     * level.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the level pointer is null.
     */
    virtual void
    resetTimeDependentData(
        const tbox::Pointer<hier::BasePatchLevel<NDIM> > level,
        const double new_time,
        const bool can_be_refined);

    /*!
     * Deallocate all new simulation data on the given level.  This
     * may be necessary during regridding, or setting up levels
     * initially.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the level pointer is null.
     */
    virtual void
    resetDataToPreadvanceState(
        const tbox::Pointer<hier::BasePatchLevel<NDIM> > level);

    /*!
     * Initialize data on a new level after it is inserted into an AMR
     * patch hierarchy by the gridding algorithm.  The level number
     * indicates that of the new level.  The old_level pointer
     * corresponds to the level that resided in the hierarchy before
     * the level with the specified number was introduced.  If the
     * pointer is null, there was no level in the hierarchy prior to
     * the call and the level data is set based on the user routines
     * and the simulation time.  Otherwise, the specified level
     * replaces the old level and the new level receives data from the
     * old level appropriately before it is destroyed.
     *
     * Typically, when data is set, it is interpolated from coarser
     * levels in the hierarchy.  If the data is to be set, the level
     * number must match that of the old level, if non-null.  If the
     * old level is non-null, then data is copied from the old level
     * to the new level on regions of intersection between those
     * levels before interpolation occurs.  Then, user-supplied patch
     * routines are called to further initialize the data if needed.
     * The boolean argument initial_time is passed into the user's
     * routines.
     *
     * The boolean argument initial_time indicates whether the level
     * is being introduced for the first time (i.e., at initialization
     * time), or after some regrid process during the calculation
     * beyond the initial hierarchy construction.  This information is
     * provided since the initialization of the data on a patch may be
     * different in each of those circumstances.  The can_be_refined
     * boolean argument indicates whether the level is the finest
     * level allowed in the hierarchy.  This may or may not affect the
     * data initialization process depending on the problem.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the hierarchy pointer is null, the level number
     * does not match any level in the hierarchy, or the old level
     * number does not match the level number (if the old level
     * pointer is non-null).
     */
    virtual void
    initializeLevelData(
        const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double init_data_time,
        const bool can_be_refined,
        const bool initial_time,
        const tbox::Pointer<hier::BasePatchLevel<NDIM> > old_level =
        tbox::Pointer<hier::BasePatchLevel<NDIM> >(NULL),
        const bool allocate_data = true);

    /*!
     * Reset cached communication schedules after the hierarchy has
     * changed (due to regidding, for example) and the data has been
     * initialized on the new levels.  The intent is that the cost of
     * data movement on the hierarchy will be amortized across
     * multiple communication cycles, if possible.  The level numbers
     * indicate the range of levels in the hierarchy that have
     * changed.  However, this routine updates communication schedules
     * every level finer than and including that indexed by the
     * coarsest level number given.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the hierarchy pointer is null, any pointer to a
     * level in the hierarchy that is coarser than the finest level is
     * null, or the given level numbers not specified properly; e.g.,
     * coarsest_level > finest_level.
     */
    virtual void
    resetHierarchyConfiguration(
        const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the
     * given level should occur according to some user-supplied
     * gradient criteria.  The double time argument is the regrid
     * time.  The integer "tag_index" argument is the patch descriptor
     * index of the cell-centered integer tag array on each patch in
     * the hierarchy.  The boolean argument initial_time indicates
     * whether the level is being subject to refinement at the initial
     * simulation time.  If it is false, then the error estimation
     * process is being invoked at some later time after the AMR
     * hierarchy was initially constructed.  The boolean argument
     * uses_richardson_extrapolation_too is true when Richardson
     * extrapolation error estimation is used in addition to the
     * gradient detector, and false otherwise.  This argument helps
     * the user to manage multiple regridding criteria.  This
     * information is passed along to the user's patch tagging
     * routines since the application of the gradient detector may be
     * different in each case.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the hierarchy pointer is null or the level
     * number does not match any existing level in the hierarchy.
     */
    virtual void
    applyGradientDetector(
        const tbox::Pointer<hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const double error_data_time,
        const int tag_index,
        const bool initial_time,
        const bool uses_richardson_extrapolation_too);

    /*!
     * Set integer tags to "one" where refinement onf the given level
     * should occur according to some user-supplied Richardson
     * extrapolation criteria.  The "error_data_time" argument is the
     * regrid time.  The "deltat" argument is the time increment to
     * advance the solution on the level to be refined.  Note that
     * that level is finer than the level in the argument list, in
     * general.  The ratio between the argument level and the actual
     * hierarchy level is given by the integer "coarsen ratio".
     *
     * The integer "tag_index" argument is the patch descriptor index
     * of the cell-centered integer tag array on each patch in the
     * hierarchy.
     *
     * The boolean argument initial_time indicates whether the level
     * is being subject to refinement at the initial simulation time.
     * If it is false, then the error estimation process is being
     * invoked at some later time after the AMR hierarchy was
     * initially constructed.  Typically, this information is passed
     * to the user's patch tagging routines since the application of
     * the Richardson extrapolation process may be different in each
     * case.
     *
     * The boolean uses_gradient_detector_too is true when a gradient
     * detector procedure is used in addition to Richardson
     * extrapolation, and false otherwise.  This argument helps the
     * user to manage multiple regridding criteria.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the level pointer is null.
     */
    virtual void
    applyRichardsonExtrapolation(
        const tbox::Pointer<hier::PatchLevel<NDIM> > level,
        const double error_data_time,
        const int tag_index,
        const double deltat,
        const int error_coarsen_ratio,
        const bool initial_time,
        const bool uses_gradient_detector_too);

    /*!
     * Coarsen solution data from level to coarse_level for Richardson
     * extrapolation.  Note that this routine will be called twice
     * during the Richardson extrapolation error estimation process.
     * The before_advance boolean argument indicates whether data is
     * set on the coarse level by coarsening the "old" time level
     * solution (i.e., before it has been advanced) or by coarsening
     * the "new" solution on the fine level (i.e., after it has been
     * advanced).
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if either level pointer is null.
     */
    virtual void
    coarsenDataForRichardsonExtrapolation(
        const tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
        const int level_number,
        const tbox::Pointer<hier::PatchLevel<NDIM> > coarse_level,
        const double coarsen_data_time,
        const bool before_advance);

    /*!
     * Register a variable with the parabolic integration algorithm.
     * The variable type must be one of the options defined by the
     * enumerated type defined above.  Typically, this routine is
     * called from the parabolic patch model when the variable
     * registration process is invoked by calling the function
     * initializeLevelIntegrator() above.  In fact, that function
     * should be called before this routine is called.
     *
     * When assertion checking is active, an unrecoverable exception
     * will result if the variable pointer or geometry pointer is
     * null.
     */
    virtual void
    registerVariable(
        const tbox::Pointer<hier::Variable<NDIM> > var,
        const hier::IntVector<NDIM> ghosts,
        const PAR_VAR_TYPE h_v_type,
        const tbox::Pointer<xfer::Geometry<NDIM> > transfer_geom,
        const std::string& coarsen_name = std::string(),
        const std::string& refine_name = std::string());

    /*!
     * Print class data representation for parabolic level integrator
     * object.  This is done automatically, when an unrecoverable
     * run-time exception is thrown within some member function of
     * this class.
     */
    virtual void printClassData(ostream& os) const;

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database point must be
     * non-null.
     */
    virtual void putToDatabase(tbox::Pointer<tbox::Database> db);

    /*!
     * Return pointer to "current" variable context used by
     * integrator.  Current data corresponds to state data at the
     * beginning of a timestep, or when a new level is initialized.
     */
    tbox::Pointer<hier::VariableContext>
    getCurrentContext() const;

    /*!
     * Return pointer to "new" variable context used by integrator.
     * New data corresponds to advanced state data at the end of a
     * timestep.  The data is one timestep later than the "current"
     * data.
     */
    tbox::Pointer<hier::VariableContext>
    getNewContext() const;

    /*!
     * Return pointer to "old" variable context used by integrator.
     * Old data corresponds to an extra time level of state data used
     * for Richardson extrapolation error estimation.  The data is one
     * timestep earlier than the "current" data.
     *
     * Note that only in certain cases when using time-dependent error
     * estimation, such as Richardson extrapolation, is the returned
     * pointer will non-null.  See contructor for more information.
     */
    tbox::Pointer<hier::VariableContext>
    getOldContext() const;

    /*!
     * Return pointer to "scratch" variable context used by
     * integrator.  Scratch data typically corresponds to storage that
     * user-routines in the concrete ParabolicPatchStrategy object
     * manipulate; in particular, scratch data contains ghost cells.
     */
    tbox::Pointer<hier::VariableContext>
    getScratchContext() const;

    /*!
     * Return pointer to variable context used for plotting.  This
     * context corresponds to the data storage that should be written
     * to plot files.  Typically, this is the same as the "current"
     * context.
     */
    tbox::Pointer<hier::VariableContext>
    getPlotContext() const;

    /*!
     * Return true if this class has been constructed to use refined
     * timestepping and false if it has been constructed to use
     * synchronized timestepping.
     */
    bool
    usingRefinedTimestepping() const;

protected:
    /*!
     * Read input values, indicated above, from given database.  The
     * boolean argument is_from_restart should be set to true if the
     * simulation is beginning from restart.  Otherwise it should be
     * set to false.
     *
     * When assertion checking is active, the database pointer must be
     * non-null.
     */
    virtual void
    getFromInput(
        tbox::Pointer<tbox::Database> db,
        bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class
     * data members.  The database from which the restart data is read
     * is determined by the object_name specified in the constructor.
     *
     * Unrecoverable Errors:
     *

     *
     *    -
     *        The database corresponding to object_name is not found
     *        in the restart file.
     *
     *    -
     *        The class version number and restart version number do not
     *        match.
     *

     *
     */
    virtual void
    getFromRestart();

    /*
     * Pre-process flux storage before advancing solution on level
     * from cur_time to new_time.  The boolean flags are used to
     * determine how flux and flux integral storage is allocated and
     * initialized.  These are needed since the advanceLevel() routine
     * is used for both level integration and time-dependent error
     * estimation.
     *
     * When assertion checking is active, the level and schedule
     * pointers must be non-null and the current time must be less
     * than the new time.
     */
    virtual void
    preprocessFluxData(
        const tbox::Pointer<hier::PatchLevel<NDIM> > level,
        const double cur_time,
        const double new_time,
        const bool regrid_advance,
        const bool first_step,
        const bool last_step);

    /*
     * Post-process flux storage after advancing solution on level.
     * The boolean flag is used to determine how flux and flux
     * integral storage is copied and de-allocated. This is needed
     * since the advanceLevel() routine is used for both level
     * integration and time-dependent error estimation.
     *
     * When assertion checking is active, the level pointer must be
     * non-null.
     */
    virtual void
    postprocessFluxData(
        const tbox::Pointer<hier::PatchLevel<NDIM> > level,
        const bool regrid_advance,
        const bool first_step,
        const bool last_step);

    /*
     * Copy time-dependent data from source space to destination
     * space.
     *
     * When assertion checking is active, the level and context
     * pointers must be non-null.
     */
    virtual void
    copyTimeDependentData(
        const tbox::Pointer<hier::PatchLevel<NDIM> > level,
        const tbox::Pointer<hier::VariableContext> src_context,
        const tbox::Pointer<hier::VariableContext> dst_context);

    /*!
     * Apply the standard AMR hyperbolic flux synchronization process
     * preserve conservation properties in the solution between the
     * fine level and the coarse level.  The sync_time argument
     * indicates the time at which the solution on the two levels is
     * being synchronized.  The variable coarse_sim_time indicates the
     * previous simulation time on the coarse level (recall the
     * conservative difference will be repeated on the coarse level
     * during the synchronization process).  After the
     * synchronization, the flux and flux integral data storage is
     * reset on the levels.
     *
     * When assertion checking is turned on, an unrecoverable
     * exception will result if either level pointer is null, the
     * levels are not consecutive in the AMR hierarchy, or the coarse
     * sim time is not less than the sync time.
     */
    virtual void
    synchronizeLevelWithCoarser(
        const tbox::Pointer<hier::PatchLevel<NDIM> > fine,
        const tbox::Pointer<hier::PatchLevel<NDIM> > coarse,
        const double sync_time,
        const double coarse_sim_time);

private:
    /*
     * The patch strategy supplies the application-specific operations
     * needed to treat data on patches in the AMR hierarchy.
     */
    ParabolicPatchStrategy<NDIM>* d_patch_strategy;

    /*
     * The object name is used as a handle to databases stored in
     * restart files and for error reporting purposes.  The boolean is
     * used to control restart file writing operations.
     */
    std::string d_object_name;
    bool d_using_time_refinement;
    bool d_registered_for_restart;

    /*
     * Courant-Friedrichs-Levy parameters for time increment
     * selection.
     */
    double d_cfl;
    double d_cfl_init;

    /*
     * Boolean flags for algorithm variations during time integration.
     *
     * d_lag_dt_computation   indicates when time increment is
     *                        computed for next step on a level.  A
     *                        value of true means that the current
     *                        solution values will be used to compute
     *                        dt.  A value of false means that dt will
     *                        be computed after the current solution
     *                        is advanced and the new solution is used
     *                        to compute dt. The default value is
     *                        true.
     *
     * d_using_ghosts_for_dt  indicates whether the time increment
     *                        computation on a patch requires ghost
     *                        cell data (e.g., if boundary conditions
     *                        are needed).  This value must be
     *                        consistent with the numerical routines
     *                        used in the parabolic patch strategy
     *                        object to calculate the time step size.
     *                        The default is true.
     */
    bool d_lag_dt_computation;
    bool d_using_ghosts_for_dt;

    /*
     * Boolean flags for indicated whether face or side data types are
     * used for fluxes (choice is determined by numerical routines in
     * parabolic patch model).
     */
    bool d_flux_is_face;
    bool d_flux_face_registered;
    bool d_flux_side_registered;

/*
 * The following communication algorithms and schedules are created
 * and maintained to manage inter-patch communication during AMR
 * integration.  The algorithms are created in the class constructor.
 * They are initialized when variables are "registered" are registered
 * with the integrator.
 */

    /*
     * The "advance" schedule is used prior to advancing a level and
     * prior to computing dt at initialization. It must be reset each
     * time a level is regridded. All ghosts are filled with TIME_DEP
     * and INPUT data at specified time. TIME_DEP data in patch
     * interiors will be filled with CURRENT_VAR values.
     */
    tbox::Pointer<xfer::RefineAlgorithm<NDIM> > d_bdry_fill_advance;
    tbox::Array<tbox::Pointer<xfer::RefineSchedule<NDIM> > > d_bdry_sched_advance;

    /*
     * The "advance new" schedule can be used twice during a time
     * integration cycle. The first is when ghost cell data is
     * required during the conservative difference process (i.e.,
     * d_using_ghosts_for_cons_diff is true).  If this is the case,
     * ghosts must be refilled before the conservative difference on a
     * coarser level during the refluxing process can take place.  See
     * synchronizeLevelWithCoarser in class ParabolicLevelIntegrator
     * second occurs when the dt calculation is not lagged and the
     * physical boundary conditions are needed to compute dt (i.e.,
     * (!d_lag_dt_computation && d_using_ghosts_for_dt_computation) is
     * true).  In either case, all ghosts are filled with TIME_DEP and
     * INPUT data at specified time.  TIME_DEP data in patch interiors
     * will be filled with values corresponding to NEW descriptor
     * indices.  See notes accompanying
     * ParabolicLevelIntegrator::advanceLevel.
     */
    tbox::Pointer<xfer::RefineAlgorithm<NDIM> > d_bdry_fill_advance_new;
    tbox::Array<tbox::Pointer<xfer::RefineSchedule<NDIM> > >
    d_bdry_sched_advance_new;

    /*
     * The "advance old" algorithm is used to fill ghosts using time
     * interpolated data from OLD_VAR and NEW_VAR on the coarser
     * hierarchy level.  It is currently only used for advancing data
     * on a temporary level during the Richardson extrapolation
     * algorithm. Use of OLD_VAR data is required only when three time
     * levels are used (i.e. d_number_time_data_levels=3).
     */
    tbox::Pointer<xfer::RefineAlgorithm<NDIM> > d_bdry_fill_advance_old;

    /*
     * Coarsen algorithms for conservative data synchronization (e.g.,
     * flux correction or refluxing).
     */
    tbox::Pointer<xfer::CoarsenAlgorithm<NDIM> > d_coarsen_fluxsum;
    tbox::Pointer<xfer::CoarsenAlgorithm<NDIM> > d_coarsen_sync_data;
    tbox::Pointer<xfer::CoarsenAlgorithm<NDIM> > d_sync_initial_data;

    /*
     * Coarsen algorithms for Richardson extrapolation.
     */
    tbox::Pointer<xfer::CoarsenAlgorithm<NDIM> > d_coarsen_rich_extrap_init;
    tbox::Pointer<xfer::CoarsenAlgorithm<NDIM> > d_coarsen_rich_extrap_final;

    /*
     * Algorithm for filling a new patch level in the hierarchy.
     */
    tbox::Pointer<xfer::RefineAlgorithm<NDIM> > d_fill_new_level;

    /*
     * Number of levels of time-dependent data that must be maintained
     * on each patch level.  This value is used to coordinate the
     * needs of the time integration and the regridding process with
     * the patch data types and descriptor indices.
     */
    int d_number_time_data_levels;

    /*
     * hier::Variable contexts and lists of variables used for data
     * management.  The contexts are set in the constructor.  Note
     * that they must be consistent with those defined by the concrete
     * subclass of the ParabolicPatchStrategy object.  The variable
     * lists and component selectors are set in the registerVariable()
     * function.
     */

    tbox::Pointer<hier::VariableContext> d_scratch;
    tbox::Pointer<hier::VariableContext> d_current;
    tbox::Pointer<hier::VariableContext> d_new;
    tbox::Pointer<hier::VariableContext> d_old;
    tbox::Pointer<hier::VariableContext> d_plot_context;

    tbox::List<tbox::Pointer<hier::Variable<NDIM> > > d_all_variables;
    tbox::List<tbox::Pointer<hier::Variable<NDIM> > > d_time_dep_variables;
    tbox::List<tbox::Pointer<hier::Variable<NDIM> > > d_flux_variables;
    tbox::List<tbox::Pointer<hier::Variable<NDIM> > > d_fluxsum_variables;

    /*
     * SCRATCH descriptor indices for (non-TEMPORARY) variables (i.e.,
     * TIME_DEP, INPUT, FLUX). Note that these are used to create
     * scratch space before ghost cells are filled on level prior to
     * advancing the data.
     */
    hier::ComponentSelector d_saved_var_scratch_data;

    /*
     * SCRATCH descriptor indices for TEMPORARY variables.  Note that
     * these are used to create scratch space on a patch-by-patch
     * basis.
     */
    hier::ComponentSelector d_temp_var_scratch_data;

    /*
     * CURRENT descriptor indices for TIME_DEP, INPUT, NO_FILL
     * variables.  Note that these are used to create storage for
     * quantities when new patches are made (e.g., during hierachy
     * initialization, before error estimation during regridding,
     * after regridding new patch levels, etc.).
     */
    hier::ComponentSelector d_new_patch_init_data;

    /*
     * NEW descriptor indices for TIME_DEP variables.  Note that these
     * are used to create space for new data before patch level is
     * advanced.
     */
    hier::ComponentSelector d_new_time_dep_data;

    /*
     * Descriptor indices for FLUX quantities and integrals of fluxes
     * (used to accumulate flux information around fine patch
     * boundaries).  Also, a boolean flag to track flux storage on
     * level 0.
     */
    hier::ComponentSelector d_flux_var_data;
    hier::ComponentSelector d_fluxsum_data;
    bool d_have_flux_on_level_zero;

    /*
     * OLD descriptor indices for TIME_DEP variables.  Note that these
     * are used only when three time levels of data are used.
     */
    hier::ComponentSelector d_old_time_dep_data;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/ParabolicIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_ParabolicLevelIntegrator
