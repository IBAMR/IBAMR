#ifndef included_AdvDiffHierarchyIntegrator
#define included_AdvDiffHierarchyIntegrator

// Filename: AdvDiffHierarchyIntegrator.h
// Last modified: <12.Aug.2009 18:21:55 griffith@boyce-griffiths-mac-pro.local>
// Created on 16 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petsc.h>

// IBAMR INCLUDES
#include <ibamr/AdvDiffHypPatchOps.h>
#include <ibamr/GodunovAdvector.h>

// IBTK INCLUDES
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonFACOperator.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/KrylovLinearSolver.h>
#include <ibtk/SetDataStrategy.h>

// SAMRAI INCLUDES
#include <CellVariable.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenSchedule.h>
#include <FACPreconditioner.h>
#include <FaceVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HyperbolicLevelIntegrator.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <RobinBcCoefStrategy.h>
#include <SAMRAIVectorReal.h>
#include <StandardTagAndInitStrategy.h>
#include <VariableContext.h>
#include <VisItDataWriter.h>
#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>
#include <tbox/Serializable.h>

// C++ STDLIB INCLUDES
#include <ostream>
#include <map>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class AdvDiffHierarchyIntegrator manages the spatial discretization
 * and time integration of scalar- and vector-valued quantities whose dynamics
 * are governed by the advection-diffusion equation.
 *
 * Each quantity \f$ Q \f$ managed by the integrator may have a unique diffusion
 * coefficient \f$ \mu \f$ and drag coefficient \f$ \lambda \f$, and may
 * optionally have a forcing term \f$ F \f$.  Only one advection velocity \f$
 * \vec{u}^{\mbox{\scriptsize ADV}} \f$ may be registered with the integrator.
 *
 * This integrator employs adaptive local spatial refinement.  All levels of the
 * patch hierarchy are synchronously integrated in time.  In particular,
 * subcycling in time is \em not performed.
 *
 * The trapezoidal rule is employed for the implicit treatment of the diffusive
 * terms.  The advective terms are discretized by the GodunovAdvector object
 * supplied to the constructor.
 *
 * \see GodunovAdvector
 * \see AdvDiffHypPatchOps
 * \see SAMRAI::algs::HyperbolicLevelIntegrator
 * \see SAMRAI::mesh::StandardTagAndInitStrategy
 * \see SAMRAI::algs::TimeRefinementIntegrator
 * \see SAMRAI::algs::TimeRefinementLevelStrategy
 */
class AdvDiffHierarchyIntegrator
    : public SAMRAI::mesh::StandardTagAndInitStrategy<NDIM>,
      public virtual SAMRAI::tbox::Serializable
{
public:
    /*!
     * The constructor for AdvDiffHierarchyIntegrator sets some default values,
     * reads in configuration information from input and restart databases, and
     * registers the integrator object with the restart manager when requested.
     *
     * When assertion checking is active, passing in any null pointer or an
     * empty std::string will result in an unrecoverable exception.
     */
    AdvDiffHierarchyIntegrator(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
        SAMRAI::tbox::Pointer<GodunovAdvector> explicit_predictor,
        bool register_for_restart=true);

    /*!
     * The destructor for AdvDiffHierarchyIntegrator unregisters the integrator
     * object with the restart manager when so registered.
     */
    virtual
    ~AdvDiffHierarchyIntegrator();

    /*!
     * Return the name of the hierarchy integrator object.
     */
    const std::string&
    getName() const;

    /*!
     * Return the type of viscous time integration scheme being employed by the
     * advection-diffusion solver.
     *
     * At the present time, valid time integration schemes include:
     *
     *    - "TGA"
     *    - "CRANK-NICOLSON"
     *    - "BACKWARD-EULER"
     *
     * \note The choice of time integration scheme employed by the solver is set
     * via the input database provided to the class constructor.
     */
    const std::string&
    getViscousTimesteppingType() const;

    /*!
     * Register a VisIt data writer so this class will write plot files that may
     * be postprocessed with the VisIt visualization tool.
     */
    void
    registerVisItDataWriter(
        SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer);

    ///
    ///  The following routines:
    ///
    ///      registerAdvectedAndDiffusedQuantity(),
    ///      registerAdvectedAndDiffusedQuantityWithSourceTerm(),
    ///      registerAdvectionVelocity()
    ///
    ///  allow the specification of quantities to be advected and diffused.
    ///

    /*!
     * Register a scalar-valued cell-centered quantity to be advected and
     * diffused according to the specified advection velocity and diffusion
     * coefficient.
     *
     * Conservative differencing is employed in evaluating the advective term
     * when conservation_form is true.  Otherwise, non-conservative differencing
     * is used to update the quantity.
     *
     * Optional concrete IBTK::SetDataStrategy and
     * SAMRAI::solv::RobinBcCoefStrategy objects allow for the specification of
     * initial and boundary data for the advected and diffused quantity Q.  If
     * an initialization object is not specified, Q is initialized to zero.  If
     * a boundary condition object is not specified for Q, it is necessary that
     * the computational domain have only periodic boundaries.  (I.e. the domain
     * can have no "physical" boundaries.)
     *
     * When the advected and diffused quantity Q is an incompressible velocity
     * field, an optional face centered gradient may be specified that
     * approximately enforces the incompressibility constraint.  The gradient is
     * subtracted from the predicted face centered and time centered values
     * prior to the computation of the advective fluxes.
     */
    void
    registerAdvectedAndDiffusedQuantity(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        const double Q_mu,
        const double Q_lambda,
        const bool conservation_form=true,
        SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> Q_init=NULL,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const Q_bc_coef=NULL,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var=NULL);

    /*!
     * Register a vector-valued cell-centered quantity to be advected and
     * diffused according to the specified advection velocity and diffusion
     * coefficient.
     *
     * Conservative differencing is employed in evaluating the advective term
     * when conservation_form is true.  Otherwise, non-conservative differencing
     * is used to update the quantity.
     *
     * Optional concrete IBTK::SetDataStrategy and
     * SAMRAI::solv::RobinBcCoefStrategy objects allow for the specification of
     * initial and boundary data for the advected and diffused quantity Q.  If
     * an initialization object is not specified, Q is initialized to zero.  If
     * a boundary condition object is not specified for Q, it is necessary that
     * the computational domain have only periodic boundaries.  (I.e. the domain
     * can have no "physical" boundaries.)
     *
     * When the advected and diffused quantity Q is an incompressible velocity
     * field, an optional face centered gradient may be specified that
     * approximately enforces the incompressibility constraint.  The gradient is
     * subtracted from the predicted face centered and time centered values
     * prior to the computation of the advective fluxes.
     */
    void
    registerAdvectedAndDiffusedQuantity(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        const double Q_mu,
        const double Q_lambda,
        const bool conservation_form=true,
        SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> Q_init=NULL,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& Q_bc_coefs=std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(),
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var=NULL);

    /*!
     * Register a scalar-valued cell-centered quantity to be advected and
     * diffused according to the specified advection velocity, diffusion
     * coefficient, and source term.
     *
     * Conservative differencing is employed in evaluating the advective term
     * when conservation_form is true.  Otherwise, non-conservative differencing
     * is used to update the quantity.
     *
     * Optional concrete IBTK::SetDataStrategy and
     * SAMRAI::solv::RobinBcCoefStrategy objects allow for the specification of
     * initial and boundary data for the advected and diffused quantity Q.  If
     * an initialization object is not specified, Q is initialized to zero.  If
     * a boundary condition object is not specified for Q, it is necessary that
     * the computational domain have only periodic boundaries.  (I.e. the domain
     * can have no "physical" boundaries.)
     *
     * The value of the source term is determined by an (optional)
     * IBTK::SetDataStrategy object.  This allows for the specification of
     * either a constant or a time-dependent source term.  If this object is not
     * provided, the source term is initialized to zero.
     *
     * When the advected and diffused quantity Q is an incompressible velocity
     * field, an optional face centered gradient may be specified that
     * approximately enforces the incompressibility constraint.  The gradient is
     * subtracted from the predicted face centered and time centered values
     * prior to the computation of the advective fluxes.
     */
    void
    registerAdvectedAndDiffusedQuantityWithSourceTerm(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        const double Q_mu,
        const double Q_lambda,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var,
        const bool conservation_form=true,
        SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> Q_init=NULL,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const Q_bc_coef=NULL,
        SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> F_set=NULL,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var=NULL);

    /*!
     * Register a vector-valued cell-centered quantity to be advected and
     * diffused according to the specified advection velocity, diffusion
     * coefficient, and source term.
     *
     * Conservative differencing is employed in evaluating the advective term
     * when conservation_form is true.  Otherwise, non-conservative differencing
     * is used to update the quantity.
     *
     * Optional concrete IBTK::SetDataStrategy and
     * SAMRAI::solv::RobinBcCoefStrategy objects allow for the specification of
     * initial and boundary data for the advected and diffused quantity Q.  If
     * an initialization object is not specified, Q is initialized to zero.  If
     * a boundary condition object is not specified for Q, it is necessary that
     * the computational domain have only periodic boundaries.  (I.e. the domain
     * can have no "physical" boundaries.)
     *
     * The value of the source term is determined by an (optional)
     * IBTK::SetDataStrategy object.  This allows for the specification of
     * either a constant or a time-dependent source term.  If this object is not
     * provided, the source term is initialized to zero.
     *
     * When the advected and diffused quantity Q is an incompressible velocity
     * field, an optional face centered gradient may be specified that
     * approximately enforces the incompressibility constraint.  The gradient is
     * subtracted from the predicted face centered and time centered values
     * prior to the computation of the advective fluxes.
     */
    void
    registerAdvectedAndDiffusedQuantityWithSourceTerm(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
        const double Q_mu,
        const double Q_lambda,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var,
        const bool conservation_form=true,
        SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> Q_init=NULL,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& Q_bc_coefs=std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(),
        SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> F_set=NULL,
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var=NULL);

    /*!
     * Register a face centered advection velocity, used by the integrator to
     * advect the cell centered quantities registered with the integrator.
     *
     * An optional IBTK::SetDataStrategy object allows for the specification
     * of a constant or time-dependent advection velocity.  If this object is
     * not provided, the advection velocity is initialized to zero.
     *
     * The value of u_is_div_free determines whether the patch strategy is able
     * to assume that the discrete divergence of the advection velocity is zero.
     */
    void
    registerAdvectionVelocity(
        SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > u_var,
        const bool u_is_div_free,
        SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> u_set=NULL);

    ///
    ///  The following routines:
    ///
    ///      getHierarchyMathOps(),
    ///      setHierarchyMathOps(),
    ///      isManagingHierarchyMathOps()
    ///
    ///  allow for the sharing of a single HierarchyMathOps object between
    ///  mutiple HierarchyIntegrator objects.
    ///

    /*!
     * Return a pointer to the HierarchyMathOps object being used by this
     * integrator.
     *
     * The HierarchyMathOps object supplies discrete differential operations on
     * the patch hierarchy as well as cell weights used in computing discrete
     * norms of quantities defined on the patch hierarchy.
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
     * IBTK::HierarchyMathOps object.
     *
     * When the integrator is managing the state of its IBTK::HierarchyMathOps
     * object, the integrator is responsible for invoking
     * IBTK::HierarchyMathOps::setPatchHierarchy() and
     * IBTK::HierarchyMathOps::resetLevels() following any changes to the
     * configuration of the patch hierarchy.
     */
    bool
    isManagingHierarchyMathOps() const;

    ///
    ///  The following routines:
    ///
    ///      initializeHierarchyIntegrator(),
    ///      initializeHierarchy(),
    ///      advanceHierarchy(),
    ///      atRegridPoint(),
    ///      getIntegratorTime(),
    ///      getStartTime(),
    ///      getEndTime(),
    ///      getIntegratorStep(),
    ///      getMaxIntegratorSteps(),
    ///      stepsRemaining(),
    ///      getPatchHierarchy(),
    ///      getGriddingAlgorithm(),
    ///      getHyperbolicLevelIntegrator(),
    ///      getHyperbolicPatchStrategy()
    ///
    ///  allow the AdvDiffHierarchyIntegrator to be used as a hierarchy
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
     * Synchronously advance each level in the hierarchy through the given time
     * increment and return an appropriate time increment for subsequent
     * advances.  The boolean argument indicates whether the coarsest hierarchy
     * level (i.e., level 0) should be load balanced before the levels are
     * advanced.  In general, the problem domain (determined by the union of
     * patches on level 0) does not change once set.  However, the boolean flag
     * here allows one to reconfigure the patches on the coarsest level which
     * constitute this union.  This may be required depending on a dynamic
     * change of the work load.  By default, the level will not be subject to
     * load balancing.
     *
     * This function assumes that all data on each level in the hierarchy has
     * been set and that only the data need for initialization exists on each
     * level (as opposed to both current and new data, for example).  Upon
     * return from this function, the simulation data on each hierarchy levels
     * is advanced through the time increment dt.  In addition, data on all
     * hierarchy levels has been synchronized so that it is consistent at the
     * new simulation time.  Thus, the data is set properly for any subsequent
     * calls to this function.
     *
     * When assertion checking is active, an unrecoverable exception will result
     * if the new time is not greater than the given time.
     */
    virtual double
    advanceHierarchy(
        const double dt);

    /*!
     * Return true if the current step count indicates that regridding should
     * occur.  In particular, true is returned if both the coarsest level allows
     * refinement and the step count is an integer multiple of the regrid step
     * interval.  Otherwise, false is returned.
     */
    bool
    atRegridPoint() const;

    /*!
     * Return the current integration time for the coarsest hierarchy level.
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
     * Return the integration step count for the entire hierarchy (i.e., number
     * of steps taken on the coarsest level).
     */
    int
    getIntegratorStep() const;

    /*!
     * Return the maximum number of integration steps allowed for the entire
     * hierarchy (i.e., steps allowed on coarsest level).
     */
    int
    getMaxIntegratorSteps() const;

    /*!
     * Return true if any steps remain in current step sequence.  Return false
     * otherwise.
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
     * Return a pointer to the SAMRAI::algs::HyperbolicLevelIntegrator being
     * used to integrate the advective terms.
     */
    SAMRAI::tbox::Pointer<SAMRAI::algs::HyperbolicLevelIntegrator<NDIM> >
    getHyperbolicLevelIntegrator() const;

    /*!
     * Return a pointer to the SAMRAI::algs::HyperbolicPatchStrategy being used
     * to specify the numerical routines used to integrate the advective terms.
     */
    SAMRAI::tbox::Pointer<AdvDiffHypPatchOps>
    getHyperbolicPatchStrategy() const;

    ///
    ///  The following routines:
    ///
    ///      regridHierarchy(),
    ///      integrateHierarchy(),
    ///      synchronizeHierarchy(),
    ///      synchronizeNewLevels(),
    ///      resetTimeDependentHierData(),
    ///      resetHierDataToPreadvanceState()
    ///
    ///  allow the AdvDiffHierarchyIntegrator to provide data management for a
    ///  time integrator which making use of this class.
    ///

    /*!
     * Regrid the hierarchy according to the error estimator specified by the
     * patch strategy.
     */
    virtual void
    regridHierarchy();

    /*!
     * Advance the data from current_time to new_time but do not synchronize the
     * data on the hierarchy.  This function assumes that all data on each level
     * in the hierarchy has been set and that only the data need for
     * initialization exists on each level (as opposed to both current and new
     * data, for example).  Upon return from this function, the simulation data
     * on each hierarchy levels is advanced through the specified time
     * increment.
     *
     * Note that data IS NOT synchronized by this routine.
     */
    virtual double
    integrateHierarchy(
        const double current_time,
        const double new_time);

    /*!
     * Coarsen new solution data on all levels of the hierarchy to synchronize
     * the data.  This operation makes the solution consistent between coarser
     * levels and finer levels.
     */
    virtual void
    synchronizeHierarchy();

    /*!
     * Coarsen current solution data from finest hierarchy level specified down
     * through the coarsest hierarchy level specified, if initial_time is true.
     * In this case, the hierarchy is being constructed at the initial
     * simulation time, After data is coarsened, initialization routines are
     * called to set data before that solution is further coarsened to the next
     * coarser level in the hierarchy.  This operation makes the solution
     * consistent between coarser levels and finer levels that did not exist
     * when the coarse levels where created and initialized originally.
     *
     * When initial_time is false, this routine does nothing since the standard
     * hyperbolic AMR algorithm for conservation laws, here adapted to the
     * advection-diffusion equation, requires no data synchronization after
     * regridding beyond interpolation of data from coarser levels in the
     * hierarchy in some conservative fashion.
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
     * Reset time-dependent data storage on the patch hierarchy.  This routine
     * is called when the current level data is no longer needed and it is
     * appropriate to replace the current data with the new data on the
     * hierarchy, if such data exists.
     */
    virtual void
    resetTimeDependentHierData(
        const double new_time);

    /*!
     * Reset data on the patch hierarchy to its state before the time advance.
     * This is needed, for example, when the integrator is embedded in a
     * nonlinear solver.  This routine is called to discard the new solution
     * data so that subsequent calls to advance are provided proper data at the
     * correct time.
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
     * example, due to regidding) and the data has been initialized on the new
     * levels.  The intent is that the cost of data movement on the hierarchy
     * will be amortized across multiple communication cycles, if possible.  The
     * level numbers indicate the range of levels in the hierarchy that have
     * changed.  However, this routine updates communication schedules every
     * level finer than and including that indexed by the coarsest level number
     * given.  When the integrator is managing the state of its
     * IBTK::HierarchyMathOps object, the integrator also invokes
     * IBTK::HierarchyMathOps::setPatchHierarchy() and
     * IBTK::HierarchyMathOps::resetLevels().
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
     * should occur according to some gradient criteria specified by the
     * GodunovAdvector object.  The double time argument is the regrid time.
     * The integer "tag_index" argument is the patch descriptor index of the
     * cell centered integer tag array on each patch in the hierarchy.  The
     * boolean argument initial_time indicates whether the level is being
     * subject to refinement at the initial simulation time.  If it is false,
     * then the error estimation process is being invoked at some later time
     * after the AMR hierarchy was initially constructed.  The boolean argument
     * uses_richardson_extrapolation_too is true when Richardson extrapolation
     * error estimation is used in addition to the gradient detector, and false
     * otherwise.  This argument helps the user to manage multiple regridding
     * criteria.  This information is passed along to the user's patch tagging
     * routines since the application of the gradient detector may be different
     * in each case.
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
    ///      getCurrentContext(),
    ///      getNewContext(),
    ///      getOldContext(),
    ///      getScratchContext(),
    ///      getPlotContext()
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
     * Return pointer to "old" variable context used by integrator.  Old data
     * corresponds to an extra time level of state data used for Richardson
     * extrapolation error estimation.  The data is one timestep earlier than
     * the "current" data.
     *
     * Note that only in certain cases when using time-dependent error
     * estimation, such as Richardson extrapolation, is the returned pointer
     * will non-null.  See contructor for more information.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
    getOldContext() const;

    /*!
     * Return pointer to "scratch" variable context used by integrator.  Scratch
     * data typically corresponds to storage that user-routines in the concrete
     * GodunovAdvector object manipulate; in particular, scratch data contains
     * ghost cells.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
    getScratchContext() const;

    /*!
     * Return pointer to variable context used for plotting.  This context
     * corresponds to the data storage that should be written to plot files.
     * Typically, this is the same as the "current" context.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
    getPlotContext() const;

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
     * String indicating the time integration employed for the implicit
     * discretization of the viscous terms.
     *
     * Choices are:
     *     - BACKWARD_EULER
     *     - CRANK_NICOLSON
     *     - TGA
     */
    std::string d_viscous_timestepping_type;

    /*!
     * Advected and diffused quantities Q, source terms F (possibly NULL),
     * advection source terms Psi = F + mu*L*Q, and optional face centered
     * gradient terms to enforce incompressibility.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_Q_vars;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_F_vars;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > > d_Psi_vars;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > > d_grad_vars;

    /*!
     * Objects to set initial and boundary conditions as well as forcing terms
     * for each advected and diffused quantity.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> > d_Q_inits;
    std::vector<std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> > d_Q_bc_coefs;
    std::vector<SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> > d_F_sets;

    /*!
     * The diffusivity and drag coefficients associated with each advected and
     * diffused quantity.
     */
    std::vector<double> d_Q_mus, d_Q_lambdas;

    /*!
     * The advection velocity.
     */
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > d_u_var;
    SAMRAI::tbox::Pointer<IBTK::SetDataStrategy> d_u_set;
    bool d_u_is_div_free;

private:
    typedef std::map<std::string,SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > >              CoarsenAlgMap;
    typedef std::map<std::string,SAMRAI::xfer::CoarsenPatchStrategy<NDIM>* >                                 CoarsenPatchStrategyMap;
    typedef std::map<std::string,std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > > CoarsenSchedMap;

    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    AdvDiffHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    AdvDiffHierarchyIntegrator(
        const AdvDiffHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    AdvDiffHierarchyIntegrator&
    operator=(
        const AdvDiffHierarchyIntegrator& that);

    /*!
     * Read input values, indicated below, from given database.  The boolean
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
     * The SAMRAI::algs::HyperbolicLevelIntegrator supplies generic operations
     * needed to handle the explicit integration of advection terms.
     */
    SAMRAI::tbox::Pointer<SAMRAI::algs::HyperbolicLevelIntegrator<NDIM> > d_hyp_level_integrator;

    /*
     * The advection patch strategy supplies the advection specific operations
     * needed to treat data on patches in the AMR hierarchy.
     */
    SAMRAI::tbox::Pointer<AdvDiffHypPatchOps> d_hyp_patch_ops;

    /*
     * Integrator data read from input or set at initialization.
     */
    double d_start_time;
    double d_end_time;
    double d_grow_dt;
    int d_max_integrator_steps;

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
     * Integrator data that evolves during time integration and maintains the
     * state of the timestep sequence over the levels in the AMR hierarchy.
     */
    double d_old_dt;
    double d_integrator_time;
    int    d_integrator_step;

    /*
     * Indicates whether the integrator has been initialized.
     */
    bool d_is_initialized;

    /*
     * Indicates whether the integrator should output logging messages.
     */
    bool d_do_log;

    /*
     * Hierarchy operations objects.
     */
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;
    bool d_is_managing_hier_math_ops;

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_wgt_var;
    int d_wgt_idx;

    /*
     * Variable context used for temporary storage.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_temp_context;

    /*
     * Communications algorithms, patch strategies, and schedules.
     */
    CoarsenAlgMap           d_calgs;
    CoarsenPatchStrategyMap d_cstrategies;
    CoarsenSchedMap         d_cscheds;

    std::vector<SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> > d_hier_bdry_fill_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_no_fill_op;

    /*
     * Linear solvers and associated data including Poisson specifications,
     * boundary conditions, and solver configuation databases.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > > d_sol_vecs, d_rhs_vecs;

    int d_max_iterations;
    double d_abs_residual_tol, d_rel_residual_tol;
    bool d_using_FAC;

    std::vector<SAMRAI::tbox::Pointer<IBTK::CCLaplaceOperator> >                d_helmholtz_ops;
    std::vector<SAMRAI::solv::PoissonSpecifications>                            d_helmholtz_specs;
    std::vector<SAMRAI::tbox::Pointer<IBTK::KrylovLinearSolver> >               d_helmholtz_solvers;
    std::vector<SAMRAI::tbox::Pointer<IBTK::CCPoissonFACOperator> >             d_helmholtz_fac_ops;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::FACPreconditioner<NDIM> > > d_helmholtz_fac_pcs;

    std::vector<bool> d_helmholtz_solvers_need_init;
    int d_coarsest_reset_ln, d_finest_reset_ln;

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_fac_op_db, d_fac_pc_db;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/AdvDiffHierarchyIntegrator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_AdvDiffHierarchyIntegrator
