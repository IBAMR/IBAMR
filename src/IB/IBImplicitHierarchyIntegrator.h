// Filename: IBImplicitHierarchyIntegrator.h
// Created on 08 May 2008 by Boyce Griffith
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

#ifndef included_IBImplicitHierarchyIntegrator
#define included_IBImplicitHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petsc.h>

// IBAMR INCLUDES
#include <ibamr/IBLagrangianForceStrategy.h>
#include <ibamr/IBImplicitJacobian.h>
#include <ibamr/IBImplicitModHelmholtzOperator.h>
#include <ibamr/IBImplicitOperator.h>
#include <ibamr/IBImplicitSFSstarOperator.h>
#include <ibamr/IBImplicitSJSstarOperator.h>
#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibamr/INSCoefs.h>
#include <ibamr/INSStaggeredBlockFactorizationPreconditioner.h>
#include <ibamr/INSStaggeredBoxRelaxationFACOperator.h>
#include <ibamr/INSStaggeredPPMConvectiveOperator.h>
#include <ibamr/INSStaggeredPhysicalBoundaryHelper.h>
#include <ibamr/INSStaggeredProjectionPreconditioner.h>
#include <ibamr/INSStaggeredStokesOperator.h>

// IBTK INCLUDES
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonFACOperator.h>
#include <ibtk/CartGridFunction.h>
#if (NDIM == 3)
#include <ibtk/LagM3DDataWriter.h>
#endif
#include <ibtk/LagSiloDataWriter.h>
#include <ibtk/FACPreconditioner.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/PETScMFFDJacobianOperator.h>
#include <ibtk/PETScNewtonKrylovSolver.h>
#include <ibtk/SCLaplaceOperator.h>
#include <ibtk/SCPoissonFACOperator.h>
#include <ibtk/SideDataSynchronization.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <ComponentSelector.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <StandardTagAndInitStrategy.h>
#include <VariableContext.h>
#include <VisItDataWriter.h>
#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Serializable.h>

// C++ STDLIB INCLUDES
#include <list>
#include <map>
#include <ostream>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBImplicitHierarchyIntegrator is a \em HIGHLY \em EXPERIMENTAL
 * implementation of a formally second-order accurate, implicit version of the
 * immersed boundary method.
 */
class IBImplicitHierarchyIntegrator
    : public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>,
      public SAMRAI::tbox::Serializable
{
public:
    friend class IBImplicitSFSstarOperator;
    friend class IBImplicitSJSstarOperator;

    /*!
     * The constructor for IBImplicitHierarchyIntegrator sets some default
     * values, reads in configuration information from input and restart
     * databases, and registers the integrator object with the restart manager
     * when requested.
     *
     * When assertion checking is active, passing in any null pointer or an
     * empty string will result in an unrecoverable exception.
     */
    IBImplicitHierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        bool register_for_restart=true);

    /*!
     * The destructor for IBImplicitHierarchyIntegrator unregisters the
     * integrator object with the restart manager when so registered.
     */
    virtual
    ~IBImplicitHierarchyIntegrator();

    /*!
     * Return the name of the hierarchy integrator object.
     */
    const std::string&
    getName() const;

    /*!
     * Supply initial conditions for the (side centered) velocity.
     */
    void
    registerVelocityInitialConditions(
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> u_init);

    /*!
     * Supply physical boundary conditions for the (side centered) velocity.
     */
    void
    registerVelocityPhysicalBcCoefs(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs);

    /*!
     * Supply initial conditions for the (cell centered) pressure.
     *
     * \note These initial conditions are used for output purposes only.  They
     * are not actually used in the computation.
     */
    void
    registerPressureInitialConditions(
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> p_init);

    /*!
     * Supply an optional side centered body force specification object.
     */
    void
    registerBodyForceSpecification(
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> f_fcn);

    /*!
     * Supply an optional cell centered source/sink specification object.
     */
    void
    registerSourceSpecification(
        SAMRAI::tbox::Pointer<IBTK::CartGridFunction> q_fcn);

    /*!
     * Supply an optional advection-diffusion solver object.
     */
    void
    registerAdvDiffHierarchyIntegrator(
        SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator);

    /*!
     * Register a concrete strategy object with the integrator that specifies
     * the initial configuration of the curvilinear mesh nodes.
     */
    void
    registerLNodeInitStrategy(
        SAMRAI::tbox::Pointer<IBTK::LNodeInitStrategy> lag_init_strategy);

    /*!
     * Free the concrete initialization strategy object.
     *
     * \note Be sure to call this method only once the initialization object is
     * no longer needed.
     */
    void
    freeLNodeInitStrategy();

    /*!
     * Register a concrete strategy object with the integrator that specifies
     * the forces generated by the curvilinear mesh nodes.
     */
    void
    registerIBLagrangianForceStrategy(
        SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> lag_force_strategy);

    /*!
     * Register a VisIt data writer so this object will write plot files that
     * may be postprocessed with the VisIt visualization tool.
     */
    void
    registerVisItDataWriter(
        SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer);

    /*!
     * Register a Lagrangian Silo data writer so this class will write plot
     * files that may be postprocessed with the VisIt visualization tool.
     */
    void
    registerLagSiloDataWriter(
        SAMRAI::tbox::Pointer<IBTK::LagSiloDataWriter> silo_writer);

#if (NDIM == 3)
    /*!
     * Register a Lagrangian myocardial3D data writer so this class will write
     * plot files that may be postprocessed with the myocardial3D visualization
     * program.
     */
    void
    registerLagM3DDataWriter(
        SAMRAI::tbox::Pointer<IBTK::LagM3DDataWriter> m3D_writer);
#endif

    /*!
     * Register a load balancer for non-uniform load balancing.
     */
    void
    registerLoadBalancer(
        SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > load_balancer);

    ///
    ///  The following routines:
    ///
    ///      registerRegridHierarchyCallback(),
    ///      registerApplyGradientDetectorCallback()
    ///
    ///  allow for the registration of callback functions that are executed by
    ///  the hierarchy integrator.
    ///

    /*!
     * \brief Callback registration function.
     */
    void
    registerRegridHierarchyCallback(
        void (*callback)(const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy, const double regrid_data_time, const bool initial_time, void* ctx),
        void* ctx);

    /*!
     * \brief Callback registration function.
     */
    void
    registerApplyGradientDetectorCallback(
        void (*callback)(const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy, const int level_number, const double error_data_time, const int tag_index, const bool initial_time, const bool uses_richardson_extrapolation_too, void* ctx),
        void* ctx);

    ///
    ///  The following routines:
    ///
    ///      getHierarchyMathOps(),
    ///      setHierarchyMathOps(),
    ///      isManagingHierarchyMathOps()
    ///
    ///  allow for the sharing of a single HierarchyMathOps object between
    ///  multiple HierarchyIntegrator objects.
    ///

    /*!
     * Return a pointer to the HierarchyMathOps object being used by this
     * integrator.
     */
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps>
    getHierarchyMathOps() const;

    /*!
     * Set the HierarchyMathOps object being used by this integrator.
     *
     * When manage_ops is true, the HierarchyMathOps object is managed by the
     * integrator.  In particular, the integrator is responsible for invoking
     * HierarchyMathOps::setPatchHierarchy() and HierarchyMathOps::resetLevels()
     * following any changes to the configuration of the patch hierarchy.
     */
    void
    setHierarchyMathOps(
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
        const bool manage_ops=false);

    /*!
     * Returns whether this integrator is managing the state of its
     * HierarchyMathOps object.
     *
     * When the integrator is managing the state of its HierarchyMathOps object,
     * the integrator is responsible for invoking
     * HierarchyMathOps::setPatchHierarchy() and HierarchyMathOps::resetLevels()
     * following any changes to the configuration of the patch hierarchy.
     */
    bool
    isManagingHierarchyMathOps() const;

    ///
    ///  The following routines:
    ///
    ///      initializeHierarchyIntegrator(),
    ///      initializeHierarchy(),
    ///      advanceHierarchy(),
    ///      getStableTimestep(),
    ///      atRegridPoint(),
    ///      getIntegratorTime(),
    ///      getStartTime(),
    ///      getEndTime(),
    ///      getIntegratorStep(),
    ///      getMaxIntegratorSteps(),
    ///      stepsRemaining(),
    ///      getPatchHierarchy(),
    ///      getGriddingAlgorithm(),
    ///      getLDataManager()
    ///
    ///  allow the IBImplicitHierarchyIntegrator to be used as a hierarchy
    ///  integrator.
    ///

    /*!
     * Initialize the variables and communications algorithms managed and used
     * by the integrator.
     *
     * This method must be called prior to any calls to initializeHierarchy() or
     * advanceHierarchy().  Otherwise, when assertion checking is active an
     * unrecoverable exception will occur.
     */
    virtual void
    initializeHierarchyIntegrator(
        SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg);

    /*!
     * Set AMR patch hierarchy configuration and data at start of simulation.
     * If the computation is begun from a restart file, the hierarchy and data
     * are read from the hierarchy database.  Otherwise, the hierarchy and data
     * are initialized by the gridding algorithm data member.  In this case, the
     * coarsest level is constructed and initialized.  Then, error estimation is
     * performed to determine if and where it should be refined.  Successively
     * finer levels are created and initialized until the maximum allowable
     * number of levels is achieved or no further refinement is needed.  The
     * double return value is the time increment for the first data advance
     * step.
     *
     * This function assumes that the hierarchy exists, but that it contains no
     * patch levels, when it is called.  On return from this function, the
     * initial hierarchy configuration and simulation data is set properly for
     * the advanceHierarchy() function to be called.  In particular, on each
     * level constructed only the data needed for initialization exists.
     */
    virtual double
    initializeHierarchy();

    /*!
     * Integrate data on all patches on all levels of the patch hierarchy from
     * current time (current_time) to new time (new_time).
     *
     * When assertion checking is active, an unrecoverable exception will result
     * if the new time is not greater than the given time.
     */
    virtual double
    advanceHierarchy(
        const double dt);

    /*!
     * Returns the maximum stable timestep according to the hyperbolic CFL
     * condition and a growth condition.  The growth condition is imposed to
     * prevent excessive changes in the maximum stable timestep as the
     * computation progresses.
     */
    virtual double
    getStableTimestep(
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx);

    /*!
     * Return true if the current step count indicates that regridding should
     * occur.
     */
    bool
    atRegridPoint() const;

    /*!
     * Return the current integration time for coarsest hierarchy level.
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
     * Return the integration step count for entire hierarchy (i.e., number of
     * steps taken on the coarsest level).
     */
    int
    getIntegratorStep() const;

    /*!
     * Return the maximum number of integration steps allowed for entire
     * hierarchy (i.e., steps allowed on coarsest level).
     */
    int
    getMaxIntegratorSteps() const;

    /*!
     * Return true if any integration steps remain, false otherwise.
     */
    bool
    stepsRemaining() const;

    /*!
     * Return a const pointer to the patch hierarchy managed by integrator.
     */
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
    getPatchHierarchy() const;

    /*!
     * Return a pointer to the gridding algorithm object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
    getGriddingAlgorithm() const;

    /*!
     * Return a pointer to the Lagrangian data manager object.
     */
    IBTK::LDataManager*
    getLDataManager() const;

    ///
    ///  The following routines:
    ///
    ///      regridHierarchy(),
    ///      integrateHierarchy_initialize(),
    ///      integrateHierarchy(),
    ///      integrateHierarchy_finalize(),
    ///      synchronizeHierarchy(),
    ///      synchronizeNewLevels(),
    ///      resetTimeDependentHierData(),
    ///      resetHierDataToPreadvanceState()
    ///
    ///  allow the IBImplicitHierarchyIntegrator to provide data management
    ///  for a time integrator which makes use of this class.
    ///

    /*!
     * Regrid the hierarchy.
     */
    virtual void
    regridHierarchy();

    /*!
     * Prepare to advance the data from current_time to new_time.
     *
     * Note that hierarchy data ARE NOT synchronized by this routine.
     */
    virtual void
    integrateHierarchy_initialize(
        const double current_time,
        const double new_time);

    /*!
     * Advance the data from current_time to new_time but do not synchronize the
     * data on the hierarchy.  This function assumes that all data on each level
     * in the hierarchy has been set and that only the data need for
     * initialization exists on each level (as opposed to both current and new
     * data, for example).
     */
    virtual void
    integrateHierarchy(
        const double current_time,
        const double new_time);

    /*!
     * Clean up data following call to integrateHierarchy().
     *
     * Note that data IS NOT synchronized by this routine.
     */
    virtual void
    integrateHierarchy_finalize(
        const double current_time,
        const double new_time);

    /*!
     * Synchronize the hierarchy.
     */
    virtual void
    synchronizeHierarchy();

    /*!
     * Coarsen current solution data from finest hierarchy level specified down
     * through the coarsest hierarchy level specified, if initial_time is true.
     * In this case, the hierarchy is being constructed at the initial
     * simulation time, After data is coarsened, the application- specific
     * initialization routine is called to set data before that solution is
     * further coarsened to the next coarser level in the hierarchy.  This
     * operation makes the solution consistent between coarser levels and finer
     * levels that did not exist when the coarse levels where created and
     * initialized originally.
     *
     * When initial_time is false, this routine does nothing since the standard
     * hyperbolic AMR algorithm for conservation laws requires no data
     * synchronization after regridding beyond interpolation of data from
     * coarser levels in the hierarchy in some conservative fashion.
     *
     * When assertion checking is active, an unrecoverable exception will result
     * if the hierarchy pointer is null, the level numbers do not properly match
     * existing levels in the hierarchy (either coarsest_level > finest_level or
     * some level is null).
     */
    virtual void
    synchronizeNewLevels(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level,
        const double sync_time,
        const bool initial_time);

    /*!
     * Reset time dependent data.
     */
    virtual void
    resetTimeDependentHierData(
        const double new_time);

    /*!
     * Deallocate all new simulation data.
     */
    virtual void
    resetHierDataToPreadvanceState();

    ///
    ///  The following routines:
    ///
    ///      initializeLevelData(),
    ///      resetHierarchyConfiguration(),
    ///      applyGradientDetector()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::mesh::StandardTagAndInitStrategy abstract base class.
    ///

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.  The level number indicates that of
     * the new level.  The old_level pointer corresponds to the level that
     * resided in the hierarchy before the level with the specified number was
     * introduced.  If the pointer is null, there was no level in the hierarchy
     * prior to the call and the level data is set based on the user routines
     * and the simulation time.  Otherwise, the specified level replaces the old
     * level and the new level receives data from the old level appropriately
     * before it is destroyed.
     *
     * Typically, when data is set, it is interpolated from coarser levels in
     * the hierarchy.  If the data is to be set, the level number must match
     * that of the old level, if non-null.  If the old level is non-null, then
     * data is copied from the old level to the new level on regions of
     * intersection between those levels before interpolation occurs.  Then,
     * user-supplied patch routines are called to further initialize the data if
     * needed.  The boolean argument initial_time is passed into the user's
     * routines.
     *
     * The boolean argument initial_time indicates whether the level is being
     * introduced for the first time (i.e., at initialization time), or after
     * some regrid process during the calculation beyond the initial hierarchy
     * construction.  This information is provided since the initialization of
     * the data on a patch may be different in each of those circumstances.  The
     * can_be_refined boolean argument indicates whether the level is the finest
     * level allowed in the hierarchy.  This may or may not affect the data
     * initialization process depending on the problem.
     *
     * When assertion checking is active, an unrecoverable exception will result
     * if the hierarchy pointer is null, the level number does not match any
     * level in the hierarchy, or the old level number does not match the level
     * number (if the old level pointer is non-null).
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
     * Reset cached communication schedules after the hierarchy has changed (for
     * example, due to regridding) and the data has been initialized on the new
     * levels.  The intent is that the cost of data movement on the hierarchy
     * will be amortized across multiple communication cycles, if possible.  The
     * level numbers indicate the range of levels in the hierarchy that have
     * changed.  However, this routine updates communication schedules every
     * level finer than and including that indexed by the coarsest level number
     * given.
     *
     * When assertion checking is active, an unrecoverable exception will result
     * if the hierarchy pointer is null, any pointer to a level in the hierarchy
     * that is coarser than the finest level is null, or the given level numbers
     * not specified properly; e.g., coarsest_level > finest_level.
     */
    virtual void
    resetHierarchyConfiguration(
        const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
        const int coarsest_level,
        const int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to some user-supplied gradient criteria.  The
     * double time argument is the regrid time.  The integer "tag_index"
     * argument is the patch descriptor index of the cell centered integer tag
     * array on each patch in the hierarchy.  The boolean argument initial_time
     * indicates whether the level is being subject to refinement at the initial
     * simulation time.  If it is false, then the error estimation process is
     * being invoked at some later time after the AMR hierarchy was initially
     * constructed.  The boolean argument uses_richardson_extrapolation_too is
     * true when Richardson extrapolation error estimation is used in addition
     * to the gradient detector, and false otherwise.  This argument helps the
     * user to manage multiple regridding criteria.  This information is passed
     * along to the user's patch tagging routines since the application of the
     * gradient detector may be different in each case.
     *
     * When assertion checking is active, an unrecoverable exception will result
     * if the hierarchy pointer is null or the level number does not match any
     * existing level in the hierarchy.
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
    ///  The following routines:
    ///
    ///      getVelocityVar(),
    ///      getPressureVar(),
    ///      getExtrapolatedPressureVar(),
    ///      getForceVar(),
    ///      getSourceVar()
    ///
    ///  allows access to the various state variables maintained by the
    ///  integrator.
    ///

    /*!
     * Return a pointer to the fluid velocity state variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >
    getVelocityVar();

    /*!
     * Return a pointer to the fluid pressure state variable.
     *
     * \note The pressure state variable is defined at time level n-1/2.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
    getPressureVar();

    /*!
     * Return a pointer to the fluid pressure variable extrapolated forward in
     * time to time level n.
     *
     * \note The pressure state variable is defined at time level n-1/2.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
    getExtrapolatedPressureVar();

    /*!
     * Return a pointer to the body force variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >
    getForceVar();

    /*!
     * Return a pointer to the source strength variable.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >
    getSourceVar();

    ///
    ///  The following routines:
    ///
    ///      getCurrentContext(),
    ///      getNewContext(),
    ///      getScratchContext()
    ///
    ///  allow access to the various variable contexts maintained by the
    ///  integrator.
    ///

    /*!
     * Return pointer to "current" variable context used by integrator.  Current
     * data corresponds to state data at the beginning of a timestep, or when a
     * new level is initialized.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
    getCurrentContext() const;

    /*!
     * Return pointer to "new" variable context used by integrator.  New data
     * corresponds to advanced state data at the end of a timestep.  The data is
     * one timestep later than the "current" data.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
    getNewContext() const;

    /*!
     * Return pointer to "scratch" variable context used by integrator.  Scratch
     * data typically corresponds to storage that user-routines in the concrete
     * GodunovAdvector object manipulate; in particular, scratch data contains
     * ghost cells.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
    getScratchContext() const;

    ///
    /// The following routines:
    ///
    ///      reinterpolateVelocity(),
    ///      reinterpolateForce()
    ///
    /// are miscellaneous utility functions.

    /*!
     * Re-interpolate the staggered velocity from cell faces to cell centers.
     */
    void
    reinterpolateVelocity(
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx);

    /*!
     * Re-interpolate the staggered body force from cell faces to cell centers.
     */
    void
    reinterpolateForce(
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx);

    ///
    ///  The following routines:
    ///
    ///      putToDatabase()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::tbox::Serializable abstract base class.
    ///

    /*!
     * Write out object state to the given database.
     *
     * When assertion checking is active, database pointer must be non-null.
     */
    virtual void
    putToDatabase(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

protected:
    /*!
     * Register a "state" variable with the integrator.  When a refine operator
     * is specified, the data for the variable is automatically maintained as
     * the patch hierarchy evolves.
     *
     * All state variables are registered with three contexts: current, new, and
     * scratch.  The current context of a state variable is maintained from
     * timestep to timestep and, optionally, as the patch hierarchy evolves.
     *
     * When a coarsen operator is specified, at the end of each timestep refined
     * regions of the new context are re-filled with the underlying fine data.
     * Whether or not a coarsen operation occurs, data in the current context is
     * then overwritten by data in the new context.
     *
     * \note If a refine operator is not specified, the data for the variable is
     * UNDEFINED following any changes to the hierarchy configuration.
     *
     * \todo Add variable registration to a HierarchyIntegrator base class.
     */
    void
    registerVariable(
        int& current_idx,
        int& new_idx,
        int& scratch_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
        const SAMRAI::hier::IntVector<NDIM>& scratch_ghosts=SAMRAI::hier::IntVector<NDIM>(0),
        const std::string& coarsen_name="NO_COARSEN",
        const std::string& refine_name="NO_REFINE");

    /*!
     * Register a "scratch" variable with the integrator.  This data IS NOT
     * maintained as the patch hierarchy evolves.
     *
     * All scratch variables are registered with the scratch context.
     *
     * \todo Add variable registration to a HierarchyIntegrator base class.
     */
    void
    registerVariable(
        int& scratch_idx,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > variable,
        const SAMRAI::hier::IntVector<NDIM>& ghosts=SAMRAI::hier::IntVector<NDIM>(0));

private:
    typedef std::map<std::string,SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > >              RefineAlgMap;
    typedef std::map<std::string,SAMRAI::xfer::RefinePatchStrategy<NDIM>* >                                 RefinePatchStrategyMap;
    typedef std::map<std::string,std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > > RefineSchedMap;

    typedef std::map<std::string,SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > >              CoarsenAlgMap;
    typedef std::map<std::string,SAMRAI::xfer::CoarsenPatchStrategy<NDIM>* >                                 CoarsenPatchStrategyMap;
    typedef std::map<std::string,std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > > CoarsenSchedMap;

    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBImplicitHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBImplicitHierarchyIntegrator(
        const IBImplicitHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBImplicitHierarchyIntegrator&
    operator=(
        const IBImplicitHierarchyIntegrator& that);

    /*!
     * Project the velocity field following a regridding operation.
     */
    void
    regridProjection();

    /*!
     * (Re-)initialize the operators and solvers used by the hierarchy
     * integrator.
     */
    void
    initializeOperatorsAndSolvers(
        const double current_time,
        const double new_time);

    /*!
     * Determine the largest stable timestep on an individual patch level.
     */
    double
    getLevelDt(
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level,
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx) const;

    /*!
     * Determine the largest stable timestep on an individual patch.
     */
    double
    getPatchDt(
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch,
        SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> ctx) const;

    /*!
     * Initialize the IBLagrangianForceStrategy object for the current
     * configuration of the curvilinear mesh.
     */
    void
    resetLagrangianForceStrategy(
        const double init_data_time,
        const bool initial_time);

    /*!
     * Set the elements of the Lagrangian vector to zero at anchored nodes of
     * the curvilinear mesh.
     */
    void
    resetAnchorPointValues(
        std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > V_data,
        const int coarsest_ln,
        const int finest_ln);

    /*!
     * Read input values, indicated above, from given database.  The boolean
     * argument is_from_restart should be set to true if the simulation is
     * beginning from restart.  Otherwise it should be set to false.
     *
     * When assertion checking is active, the database pointer must be non-null.
     */
    void
    getFromInput(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
        bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.  The database from which the restart data is read is determined
     * by the object_name specified in the constructor.
     *
     * Unrecoverable Errors:
     *
     *    -   The database corresponding to object_name is not found in the
     *        restart file.
     *
     *    -   The class version number and restart version number do not match.
     *
     */
    void
    getFromRestart();

    /*
     * The object name is used as a handle to databases stored in restart files
     * and for error reporting purposes.  The boolean is used to control restart
     * file writing operations.
     */
    std::string d_object_name;
    bool d_registered_for_restart;

    /*
     * Pointers to the patch hierarchy and gridding algorithm objects associated
     * with this time integration object.
     *
     * The gridding algorithm provides grid generation and regridding routines
     * for the AMR hierarchy.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > d_gridding_alg;

    /*
     * We cache pointers to the visualization data writer to register plot
     * variables.
     *
     * Double precision values are (optional) factors used to rescale the
     * velocity, pressure, and force for plotting.
     *
     * Boolean values indicates whether to output various quantities for
     * visualization.
     */
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > d_visit_writer;
    double d_u_scale, d_p_scale, d_f_scale, d_q_scale;
    bool d_output_u, d_output_p, d_output_f, d_output_q, d_output_omega, d_output_div_u;
    SAMRAI::tbox::Pointer<IBTK::LagSiloDataWriter> d_silo_writer;
#if (NDIM == 3)
    SAMRAI::tbox::Pointer<IBTK::LagM3DDataWriter> d_m3D_writer;
#endif

    /*
     * We cache a pointer to the load balancer.
     */
    SAMRAI::tbox::Pointer<SAMRAI::mesh::LoadBalancer<NDIM> > d_load_balancer;

    /*
     * Integrator data read from input or set at initialization.
     */
    double d_start_time;
    double d_end_time;
    double d_grow_dt;
    int d_max_integrator_steps;

    /*
     * The number of cycles of fixed-point iteration to use per timestep.
     */
    int d_num_cycles;

    /*
     * The regrid interval indicates the number of integration steps taken
     * between invocations of the regridding process.
     */
    int d_regrid_interval;

    /*
     * The tag buffer indicates the number of cells on each level by which
     * tagged cells will be buffered after they have selected for refinement.
     * These values are passed into the gridding algorithm routines during
     * hierarchy construction and regridding.  The tag buffer helps to guarantee
     * that refined cells near important features in the solution will remain
     * refined until the level is regridded next.
     */
    bool d_using_default_tag_buffer;
    SAMRAI::tbox::Array<int> d_tag_buffer;

    /*
     * This boolean value determines whether the advection term is computed
     * using conservative or non-conservative differencing.
     */
    bool d_conservation_form;

    /*
     * Tag cells based on the relative and absolute magnitudes of the local
     * vorticity.
     */
    bool d_using_vorticity_tagging;
    SAMRAI::tbox::Array<double> d_omega_rel_thresh, d_omega_abs_thresh;
    double d_omega_max;

    /*
     * This boolean value determines whether the pressure is normalized to have
     * zero mean (i.e., discrete integral) at the end of each timestep.
     */
    bool d_normalize_pressure;

    /*
     * Integrator data that evolves during time integration and maintains the
     * state of the timestep sequence over the levels in the AMR hierarchy.
     */
    double d_old_dt, d_op_and_solver_init_dt;
    double d_integrator_time;
    int    d_integrator_step;

    /*
     * The CFL number.
     */
    double d_cfl;

    /*
     * A maximum timestep constraint over the specified time interval.
     */
    double d_dt_max;
    double d_dt_max_time_max;
    double d_dt_max_time_min;

    /*
     * Indicates whether the integrator has been initialized.
     */
    bool d_is_initialized;

    /*
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log;

    /*
     * The fluid density (rho), dynamic viscosity (mu), kinematic viscosity
     * (nu), and (optional) drag coefficient (lambda).
     *
     * \note rho_water = 1.00 g cm^-3
     *       mu_water  = 0.01 g cm^-1 s^-1
     *       nu_water  = 0.01 cm^2 s^-1
     */
    double d_rho, d_mu, d_lambda;
    INSCoefs* d_problem_coefs;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;
    bool d_is_managing_hier_math_ops;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_wgt_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_wgt_sc_var;
    double d_volume;

    /*
     * Hierarchy operators and solvers.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_u_scratch_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_u_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_u_half_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_n_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_p_scratch_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_p_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_sol_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_rhs_vec;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_nul_vec;

    bool d_stokes_op_needs_init;
    SAMRAI::tbox::Pointer<INSStaggeredStokesOperator> d_stokes_op;

    bool d_convective_op_needs_init;
    SAMRAI::tbox::Pointer<INSStaggeredPPMConvectiveOperator> d_convective_op;

    bool d_helmholtz_solver_needs_init;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>          d_helmholtz_hypre_pc_db, d_helmholtz_fac_pc_db;
    SAMRAI::tbox::Pointer<IBTK::SCLaplaceOperator>         d_helmholtz_op;
    SAMRAI::tbox::Pointer<IBImplicitModHelmholtzOperator>  d_mod_helmholtz_op;
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

    bool d_projection_pc_needs_init;
    SAMRAI::tbox::Pointer<INSStaggeredProjectionPreconditioner> d_projection_pc;

    bool d_vanka_pc_needs_init;
    SAMRAI::tbox::Pointer<INSStaggeredBoxRelaxationFACOperator> d_vanka_fac_op;
    SAMRAI::tbox::Pointer<IBTK::FACPreconditioner>              d_vanka_fac_pc;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>               d_vanka_fac_pc_db;

    bool d_block_pc_needs_init;
    SAMRAI::tbox::Pointer<INSStaggeredBlockFactorizationPreconditioner> d_block_pc;

    bool d_ib_op_needs_init;
    SAMRAI::tbox::Pointer<IBImplicitSFSstarOperator> d_ib_SFSstar_op;
    SAMRAI::tbox::Pointer<IBTK::JacobianOperator> d_ib_SJSstar_op;
    SAMRAI::tbox::Pointer<IBImplicitOperator> d_ib_op;

    bool d_ib_solver_needs_init;
    SAMRAI::tbox::Pointer<IBTK::PETScNewtonKrylovSolver> d_ib_solver;
    SAMRAI::tbox::Pointer<IBImplicitJacobian>            d_ib_jac_op;

    bool d_needs_regrid_projection;
    double d_regrid_max_div_growth_factor;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>        d_regrid_projection_fac_pc_db;
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>        d_regrid_projection_bc_coef;
    SAMRAI::tbox::Pointer<IBTK::CCLaplaceOperator>       d_regrid_projection_op;
    SAMRAI::solv::PoissonSpecifications*                 d_regrid_projection_spec;
    SAMRAI::tbox::Pointer<IBTK::CCPoissonFACOperator>    d_regrid_projection_fac_op;
    SAMRAI::tbox::Pointer<IBTK::FACPreconditioner>       d_regrid_projection_fac_pc;
    SAMRAI::tbox::Pointer<IBTK::PETScKrylovLinearSolver> d_regrid_projection_solver;

    /*
     * Communications algorithms, patch strategies, and schedules.
     */
    RefineAlgMap           d_ralgs;
    RefinePatchStrategyMap d_rstrategies;
    RefineSchedMap         d_rscheds;

    CoarsenAlgMap           d_calgs;
    CoarsenPatchStrategyMap d_cstrategies;
    CoarsenSchedMap         d_cscheds;

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_fill_after_regrid;

    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_u_bdry_bc_fill_op, d_no_fill_op;

    SAMRAI::tbox::Pointer<IBTK::SideDataSynchronization> d_side_synch_op;

    /*
     * Objects to set initial conditions (note that the initial value of the
     * pressure is only used for visualization), boundary conditions, and body
     * forcing.
     */
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_u_init, d_p_init;
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* d_default_u_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_u_bc_coefs;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_u_star_bc_coefs;
    SAMRAI::tbox::Pointer<INSStaggeredPhysicalBoundaryHelper> d_u_bc_helper;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_p_bc_coef;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_phi_bc_coef;
    SAMRAI::tbox::Pointer<IBTK::CartGridFunction> d_f_fcn, d_q_fcn;

    /*
     * Optional advection-diffusion solver which will use the computed
     * incompressible velocity field to advect and diffuse associated
     * quantities.
     */
    SAMRAI::tbox::Pointer<AdvDiffHierarchyIntegrator> d_adv_diff_hier_integrator;

    /*
     * SAMRAI::hier::Variable lists and SAMRAI::hier::ComponentSelector objects
     * are used for data management.
     */
    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_state_variables;
    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_scratch_variables;

    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_copy_scratch_to_current_fast;
    std::list<SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > > d_copy_scratch_to_current_slow;

    SAMRAI::hier::ComponentSelector d_current_data;
    SAMRAI::hier::ComponentSelector d_new_data;
    SAMRAI::hier::ComponentSelector d_scratch_data;
    SAMRAI::hier::ComponentSelector d_fill_after_regrid_bc_idxs;

    /*!
     * Variable contexts.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_current_context, d_new_context, d_scratch_context;

    /*!
     * State and scratch variables.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_u_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_u_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > d_u_fc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_p_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_p_extrap_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_f_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_f_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_q_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_omega_var;
#if (NDIM == 3)
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_omega_norm_var;
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_div_u_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_phi_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_u_half_ib_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_u_regrid_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_u_src_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_indicator_var;

    /*
     * Patch data descriptor indices for all "state" variables managed by the
     * integrator.
     *
     * State variables have three contexts: current, scratch, and new.
     */
    int          d_u_current_idx,          d_u_new_idx,          d_u_scratch_idx;
    int       d_u_cc_current_idx,       d_u_cc_new_idx,       d_u_cc_scratch_idx;
    int          d_p_current_idx,          d_p_new_idx,          d_p_scratch_idx;
    int   d_p_extrap_current_idx,   d_p_extrap_new_idx,   d_p_extrap_scratch_idx;
    int          d_f_current_idx,          d_f_new_idx,          d_f_scratch_idx;
    int       d_f_cc_current_idx,       d_f_cc_new_idx,       d_f_cc_scratch_idx;
    int          d_q_current_idx,          d_q_new_idx,          d_q_scratch_idx;
    int      d_omega_current_idx,      d_omega_new_idx,      d_omega_scratch_idx;
#if (NDIM == 3)
    int d_omega_norm_current_idx, d_omega_norm_new_idx, d_omega_norm_scratch_idx;
#endif
    int      d_div_u_current_idx,      d_div_u_new_idx,      d_div_u_scratch_idx;

    /*
     * Patch data descriptor indices for all "scratch" variables managed by the
     * integrator.
     *
     * Scratch variables have only one context: scratch.
     */
    int d_phi_idx, d_u_half_ib_idx, d_u_regrid_idx, d_u_src_idx, d_indicator_idx;

    /*
     * Patch data descriptors for all variables managed by the HierarchyMathOps
     * class.
     *
     * Such variables have only one context.
     */
    int d_wgt_cc_idx, d_wgt_sc_idx;

    /*!
     * \brief Callback function pointers and callback contexts.
     */
    std::vector<void (*)(const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy, const double regrid_data_time, const bool initial_time, void* ctx)> d_regrid_hierarchy_callbacks;
    std::vector<void*> d_regrid_hierarchy_callback_ctxs;

    std::vector<void (*)(const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy, const int level_number, const double error_data_time, const int tag_index, const bool initial_time, const bool uses_richardson_extrapolation_too, void* ctx)> d_apply_gradient_detector_callbacks;
    std::vector<void*> d_apply_gradient_detector_callback_ctxs;

    /*
     * The name of the discrete delta function to employ for interpolation and
     * spreading.
     */
    std::string d_delta_fcn;
    SAMRAI::hier::IntVector<NDIM> d_ghosts;

    /*
     * The LDataManager is used to coordinate the distribution of Lagrangian
     * data on the patch hierarchy.
     */
    IBTK::LDataManager* d_lag_data_manager;

    /*
     * The specification and initialization information for the Lagrangian data
     * used by the integrator.
     */
    SAMRAI::tbox::Pointer<IBTK::LNodeInitStrategy> d_lag_init_strategy;

    /*
     * The force generators.
     */
    SAMRAI::tbox::Pointer<IBLagrangianForceStrategy> d_lag_force_strategy;
    bool d_lag_force_strategy_needs_init;

    /*
     * Lagrangian data.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > d_X_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > d_X_mid_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > d_X_new_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > d_X_half_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > d_U_half_data;
    std::vector<SAMRAI::tbox::Pointer<IBTK::LNodeLevelData> > d_F_half_data;

    /*
     * List of local indices of local anchor points.
     *
     * NOTE: IB points are automatically considered to be anchored if they are
     * within 2.0*sqrt(epsilon_mach) of the physical boundary.
     */
    std::vector<std::set<int > > d_anchor_point_local_idxs;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBImplicitHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBImplicitHierarchyIntegrator
