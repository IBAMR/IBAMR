#ifndef included_INSStaggeredBoxRelaxationFACOperator
#define included_INSStaggeredBoxRelaxationFACOperator

// Filename: INSStaggeredBoxRelaxationFACOperator.h
// Last modified: <11.Jun.2010 18:25:30 griffith@boyce-griffiths-mac-pro.local>
// Created on 11 Jun 2010 by Boyce Griffith (griffith@boyce-griffiths-mac-pro.local)

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscmat.h>

// IBAMR INCLUDES
#include <ibamr/INSCoefs.h>

// IBTK INCLUDES
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/CoarseFineBoundaryRefinePatchStrategy.h>

// SAMRAI INCLUDES
#include <CoarsenAlgorithm.h>
#include <FACOperatorStrategy.h>
#include <FACPreconditioner.h>
#include <LocationIndexRobinBcCoefs.h>
#include <RefineAlgorithm.h>

// C++ STDLIB INCLUDES
#include <map>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredBoxRelaxationFACOperator is a concrete
 * SAMRAI::solv::FACOperatorStrategy implementing a box relaxation (Vanka-type)
 * smoother for use as a multigrid preconditioner.
 */
class INSStaggeredBoxRelaxationFACOperator
    : public SAMRAI::solv::FACOperatorStrategy<NDIM>
{
public:
    /*!
     * \brief Constructor.
     */
    INSStaggeredBoxRelaxationFACOperator(
        const std::string& object_name,
        const INSCoefs& problem_coefs,
        const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db=NULL);

    /*!
     * \brief Virtual destructor.
     */
    virtual ~INSStaggeredBoxRelaxationFACOperator();

    /*!
     * \name Functions for specifying the problem coefficients.
     */
    //\{

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a U_bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be NULL; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param U_bc_coefs  Vector of pointers to objects that can set the Robin boundary condition coefficients
     * \param P_bc_coef   Pointers to object that can set the Robin boundary condition coefficients
     */
    void
    setPhysicalBcCoefs(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef);

    /*!
     * \brief Specify whether the boundary conditions are homogeneous.
     */
    void
    setHomogeneousBc(
        const bool homogeneous_bc);

    /*!
     * \brief Set the current time interval.
     */
    void
    setTimeInterval(
        const double current_time,
        const double new_time);

    /*!
     * \brief Set the SAMRAI::solv::FACPreconditioner that is using this
     * concrete SAMRAI::solv::FACOperatorStrategy object.
     *
     * \param preconditioner  Pointer to the FAC preconditioner that is using this concrete FAC strategy
     */
    void
    setPreconditioner(
        const SAMRAI::solv::FACPreconditioner<NDIM>* preconditioner);

    //\}

    /*!
     * \name Functions for configuring the solver.
     */
    //\{

    /*!
     * \brief Specify the levels that need to be reset the next time the
     * operator is re-initialized.
     *
     * When the operator is initialized, then only the specified range of levels
     * are reset in the operator state the next time that the operator is
     * initialized.  If the operator is not initialized, this method has no
     * effect.
     *
     * To ensure the range of levels that is reset includes all levels in the
     * patch hierarchy, use \a coarsest_ln = \a finest_ln = \p -1.
     *
     * \note This function is used to save some unnecessary computations when
     * the hierarchy is regridded.  The range of levels specified must include
     * all levels which need to be reset by
     * SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration().
     * Any data residing outside of this range of levels will not be reset.
     * This \b is \b not what you want to have happen if, for instance, the
     * Poisson specifications changes.
     */
    void
    setResetLevels(
        const int coarsest_ln,
        const int finest_ln);

    /*!
     * \brief Specify the ghost cell width for \em both the solution and the
     * right-hand side patch data.
     */
    void
    setGhostCellWidth(
        const SAMRAI::hier::IntVector<NDIM>& ghost_cell_width);

    /*!
     * \brief Specify the smoother type.
     *
     * Select from:
     * - \c "additive"
     *
     * \note The smoother is always additive between processors ("processor
     * block Gauss-Seidel").
     */
    void
    setSmootherChoice(
        const std::string& smoother_choice);

    /*!
     * \brief Specify the coarse level solver.
     *
     * Select from:
     * - \c "block_jacobi"
     */
    void
    setCoarsestLevelSolverChoice(
        const std::string& coarse_solver_choice);

    /*!
     * \brief Set tolerance for coarse level solve.
     *
     * If the coarse level solver requires a tolerance, the specified value is
     * used.
     *
     * \note Currently, this value is not used.
     */
    void
    setCoarsestLevelSolverTolerance(
        double coarse_solver_tol);

    /*!
     * \brief Set the maximum number of iterations for the coarsest level solve.
     */
    void
    setCoarsestLevelSolverMaxIterations(
        int coarse_solver_max_its);

    /*!
     * \brief Set the name of the prolongation methods.
     */
    void
    setProlongationMethods(
        const std::string& U_prolongation_method,
        const std::string& P_prolongation_method);

    /*!
     * \brief Set the name of the restriction methods.
     */
    void
    setRestrictionMethods(
        const std::string& U_restriction_method,
        const std::string& P_restriction_method);

    /*!
     * \brief Set the maximum number of FAC cycles used by the
     * SAMRAI::solv::FACPreconditioner that employs this concrete
     * SAMRAI::solv::FACOperatorStrategy.
     */
    void
    setFACPreconditionerMaxCycles(
        int fac_max_cycles);

    /*!
     * \brief Set whether the SAMRAI::solv::FACPreconditioner that employs this
     * concrete SAMRAI::solv::FACOperatorStrategy is employing presmoothing.
     */
    void
    setFACPreconditionerUsesPresmoothing(
        bool fac_uses_presmoothing);

    /*!
     * \brief Set whether the SAMRAI::solv::FACPreconditioner that employs this
     * concrete SAMRAI::solv::FACOperatorStrategy uses a nonzero initial guess.
     */
    void
    setFACPreconditionerInitialGuessNonzero(
        bool fac_initial_guess_nonzero);

    /*!
     * \brief Set whether to skip restricting the solution.
     *
     * If this operator is being used as a preconditioner, it may be possible to
     * avoid restricting the solution.
     */
    void
    setSkipRestrictSolution(
        bool skip_restrict_sol);

    /*!
     * \brief Set whether to skip restricting the residual.
     *
     * If this operator is being used as a preconditioner, it may be possible to
     * avoid restricting the residual.
     */
    void
    setSkipRestrictResidual(
        bool skip_restrict_residual);

    //\}

    ///
    ///  The following routines:
    ///
    ///      restrictSolution(),
    ///      restrictResidual(),
    ///      prolongErrorAndCorrect(),
    ///      smoothError(),
    ///      solveCoarsestLevel(),
    ///      computeCompositeResidualOnLevel(),
    ///      computeResidualNorm(),
    ///      initializeOperatorState(),
    ///      deallocateOperatorState()
    ///
    ///  are concrete implementations of functions declared in the
    ///  SAMRAI::solv::FACOperatorStrategy abstract base class.
    ///

    /*!
     * \name Implementation of SAMRAI::solv::FACOperatorStrategy interface.
     */
    //\{

    /*!
     * \brief Restrict the solution quantity to the specified level from the
     * next finer level.
     *
     * Restrict the residual data to level dst_ln in the destination vector dst,
     * from level dst_ln+1 in the source vector src.
     *
     * Can assume:
     * - dst_ln is not the finest level in the range being solved.
     * - corresponding solution has been computed on level dst_ln+1.
     * - the source and destination residual vectors (src and dst) may or may
     *   not be the same.  (This function must work in either case.)
     *
     * Upon return from this function, the solution on the refined region of the
     * coarse level will represent the coarsened version of the fine solution in
     * a manner that is consistent with the linear system approximation on the
     * composite grid.  This function must not change the solution values
     * anywhere except on level dst_ln of the destination vector.
     *
     * The source and destination vectors may be the same.
     *
     * \param src source solution
     * \param dst destination solution
     * \param dst_ln destination level number
     */
    virtual void
    restrictSolution(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& src,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dst,
        int dst_ln);

    /*!
     * \brief Restrict the residual quantity to the specified level from the
     * next finer level.
     *
     * Restrict the residual data to level dst_ln in the destination vector dst,
     * from level dst_ln+1 in the source vector src.
     *
     * Can assume:
     * - dst_ln is not the finest level in the range being solved.
     * - corresponding residual has been computed on level dst_ln+1.
     * - the source and destination residual vectors (src and dst) may or may
     *   not be the same.  (This function must work in either case.)
     *
     * Upon return from this function, the residual on the refined region of the
     * coarse level will represent the coarsened version of the fine residual in
     * a manner that is consistent with the linear system approximation on the
     * composite grid.  This function must not change the residual values
     * anywhere except on level dst_ln of the destination vector.
     *
     * The source and destination vectors may be the same.
     *
     * \param src source residual
     * \param dst destination residual
     * \param dst_ln destination level number
     */
    virtual void
    restrictResidual(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& src,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dst,
        int dst_ln);

    /*!
     * \brief Prolong the error quantity to the specified level from the next
     * coarser level and apply the correction to the fine-level error.
     *
     * On the part of the coarse level that does \em not overlap the fine level,
     * the error is the correction to Au=f.
     *
     * On the part of the coarse level that \em does overlap the fine level, the
     * error is the correction to Ae=r of the fine level.
     *
     * This function should apply the coarse-level correction to the fine level,
     * that is \f$ e^{\mbox{\scriptsize fine}} \leftarrow e^{\mbox{\scriptsize
     * fine}} + I^{\mbox{\scriptsize fine}}_{\mbox{\scriptsize coarse}}
     * e^{\mbox{\scriptsize coarse}} \f$
     *
     * \b Note: You probably have to store the refined error in a temporary
     * location before adding it to the current error.
     *
     * The array of boundary information contains a description of the
     * coarse-fine level boundary for each patch on the level; the boundary
     * information for patch N is obtained as the N-th element in the array,
     * coarse_fine_boundary[N].
     *
     * Upon return from this function, the error on the fine level must
     * represent the correction to the solution on that level.  Also, this
     * function must not change the error values on the coarse level.
     *
     * The source and destination vectors may be the same.
     *
     * \param src source error vector
     * \param dst destination error vector
     * \param dst_ln destination level number of data transfer
     */
    virtual void
    prolongErrorAndCorrect(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& src,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dst,
        int dst_ln);

    /*!
     * \brief Perform a given number of relaxations on the error.
     *
     * Relax the residual equation Ae=r by applying the given number of
     * smoothing sweeps on the specified level.  The relaxation may ignore the
     * possible existence of finer levels on a given level.
     *
     * \param error error vector
     * \param residual residual vector
     * \param level_num level number
     * \param num_sweeps number of sweeps to perform
     */
    virtual void
    smoothError(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& error,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        int level_num,
        int num_sweeps);

    /*!
     * \brief Solve the residual equation Ae=r on the coarsest level in the FAC
     * iteration.
     *
     * This routine must fill boundary values for given solution quantity on all
     * patches on the specified level before the solve is performed.
     *
     * \param error error vector
     * \param residual residual vector
     * \param coarsest_ln coarsest level number
     */
    virtual int
    solveCoarsestLevel(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& error,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        int coarsest_ln);

    /*!
     * \brief Compute composite grid residual on a single level.
     *
     * For the specified level number level_num, compute the \em composite
     * residual r=f-Au, where f is the right hand side and u is the solution.
     * Note that the composite residual is not a one-level residual.  It must
     * take into account the composite grid stencil around the coarse-fine grid
     * interface.
     *
     * May assume:
     * - Composite residual on next finer level, level_num+1, has been computed
     *   already.
     * - If any intermediately computed data is needed from level level_num+1,
     *   it has been done and stored on that level.
     * - Residual computations for the original equation and the error equations
     *   will not be intermingled within one FAC cycle.
     *
     * Steps:
     * -# Fill boundary ghosts.
     * -# If needed, coarsen intermediate data from level level_num+1.
     * -# Compute residual \f$ r^{ln} \leftarrow f - A u^{ln} \f$.
     *
     * Final step before leaving function:
     * - If any intermediately computed data is needed in at level level_num-1,
     *   it must be computed and stored before leaving this function.
     *
     * \b Important: Do not restrict residual from finer levels.  (However, you
     * must write the function restrictResidual() to do this.)
     *
     * \b Important: This function must also work when the right-hand-side and
     * the residual are identical.  In that case, it should effectively do \f$ r
     * \leftarrow r - A u \f$.
     *
     * \param residual residual vector
     * \param solution solution vector
     * \param rhs source (right hand side) vector
     * \param level_num level number
     * \param error_equation_indicator flag stating whether u is an error
     * vector or a solution vector
     *
     * \b Note: The residual needs to be computed in two different cases:
     *
     * -# before each FAC sweep commences, and
     * -# after performing any presmoothing sweeps
     *
     * If the FAC preconditioner
     *
     * -# (a) does not use presmoothing,
     * -# (b) uses a zero initial guess, and
     * -# (c) only employs one FAC sweep (as is often the case, for instance,
     *        when the preconditioner is being used in conjunction with a Krylov
     *        subspace method),
     *
     * then we simply set the residual equal to the right hand side.  This
     * avoids some unnecessary computation as well as parallel communication.
     *
     * \b IMPORTANT: We assume that the FAC algorithm being used does not use
     * the composite residual in postsweeps.
     */
    virtual void
    computeCompositeResidualOnLevel(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs,
        int level_num,
        bool error_equation_indicator);

    /*!
     * \brief Compute the norm of the residual.
     *
     * Compute norm of the given residual on the given range of hierarchy
     * levels.  The residual vector is computed already and you should \em not
     * change it.  The only purpose of this function to allow you to choose how
     * to define the norm.
     *
     * The norm value is used during the FAC iteration to determine
     * convergence of the composite grid linear system.
     *
     * Residual values that lie under a finer level should not be counted.
     *
     * \param residual residual vector
     * \param fine_ln finest level number
     * \param coarse_ln coarsest level number
     *
     * \return norm value of residual vector, which should be non-negative
     */
    virtual double
    computeResidualNorm(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        int fine_ln,
        int coarse_ln);

    /*!
     * \brief Compute hierarchy-dependent data.
     *
     * Note that although the vector arguments given to other methods in this
     * class may not necessarily be the same as those given to this method,
     * there will be similarities, including:
     *
     * - hierarchy configuration (hierarchy pointer and level range)
     * - number, type and alignment of vector component data
     * - ghost cell width of data in the solution (or solution-like) vector
     *
     * \param solution solution vector u
     * \param rhs right hand side vector f
     */
    virtual void
    initializeOperatorState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs);

    /*!
     * \brief Remove all hierarchy-dependent data.
     *
     * Remove all hierarchy-dependent data set by initializeOperatorState().
     *
     * \see initializeOperatorState
     */
    virtual void
    deallocateOperatorState();

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredBoxRelaxationFACOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredBoxRelaxationFACOperator(
        const INSStaggeredBoxRelaxationFACOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredBoxRelaxationFACOperator& operator=(
        const INSStaggeredBoxRelaxationFACOperator& that);

    /*!
     * \name For executing, caching and resetting communication schedules.
     */
    //\{

    /*!
     * \brief Execute a refinement schedule for prolonging data.
     */
    void
    xeqScheduleProlongation(
        const int dst_idx,
        const int src_idx,
        const int dst_ln,
        const bool homogeneous_bc);

    /*!
     * \brief Execute schedule for restricting solution to the specified level.
     */
    void
    xeqScheduleURestriction(
        const int dst_idx,
        const int src_idx,
        const int dst_ln);

    /*!
     * \brief Execute schedule for restricting residual to the specified level.
     */
    void
    xeqScheduleRRestriction(
        const int dst_idx,
        const int src_idx,
        const int dst_ln);

    /*!
     * \brief Execute schedule for filling ghosts on the specified level.
     */
    void
    xeqScheduleGhostFill(
        const int dst_idx,
        const int dst_ln,
        const bool homogeneous_bc);

    /*!
     * \brief Execute schedule for filling ghosts on the specified level.
     */
    void
    xeqScheduleGhostFillNoCoarse(
        const int dst_idx,
        const int dst_ln,
        const bool homogeneous_bc);

    /*!
     * \brief Execute schedule for synchronizing data on the specified level.
     */
    void
    xeqScheduleSideDataSynch(
        const int dst_idx,
        const int dst_ln);

    //\}

    /*!
     * \brief Construct a matrix corresponding to the linear operator restricted
     * to a single patch.
     */
    static void
    buildPatchOperator(
        Mat& A,
        const INSCoefs& problem_coefs,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const int component_axis,
        const SAMRAI::hier::IntVector<NDIM>& ghost_cell_width);

    /*!
     * \brief Check to make sure that all of the options make sense.
     */
    void
    sanityCheck();

    /*
     * The object name is used for error reporting purposes.
     *
     * The boolean indicates whether this object has been initialized.
     */
    std::string d_object_name;
    bool d_is_initialized;

    /*!
     * \name Hierarchy-dependent objects.
     */
    //\{

    /*
     * Solution and rhs vectors.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_solution, d_rhs;

    /*
     * Ghost cell width.
     */
    SAMRAI::hier::IntVector<NDIM> d_gcw;

    /*
     * Mappings from patch indices to patch operators.
     */
    std::vector<std::map<int,std::vector<Vec> > > d_patch_vec_e, d_patch_vec_f;
    std::vector<std::map<int,std::vector<Mat> > > d_patch_mat;
    std::vector<std::map<int,std::vector<SAMRAI::hier::BoxList<NDIM> > > > d_patch_bc_box_overlap;

    /*
     * Reference patch hierarchy and range of levels involved in the solve.
     *
     * This variable is non-null between the initializeOperatorState() and
     * deallocateOperatorState() calls.  It is not truly needed, because the
     * hierarchy is obtainable through variables in most function argument
     * lists.  We use it to enforce working on one hierarchy at a time.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    /*
     * Range of levels to be reset the next time the operator is initialized.
     */
    bool d_in_initialize_operator_state;
    int d_coarsest_reset_ln, d_finest_reset_ln;

    //\}

    /*!
     * \name Private state variables for solution process.
     */
    //\{

    /*
     * Problem coefficient specifications.
     */
    const INSCoefs& d_problem_coefs;

    /*
     * The kind of smoothing to perform.
     */
    std::string d_smoother_choice;

    /*
     * The name of the refinement operators used to prolong the coarse grid
     * correction.
     */
    std::string d_U_prolongation_method, d_P_prolongation_method;

    /*
     * The name of the coarsening operators used to restrict the fine grid error
     * or residual.
     */
    std::string d_U_restriction_method, d_P_restriction_method;

    /*
     * SAMRAI::tbox::Pointer to the SAMRAI::solv::FACPreconditioner that is
     * using this operator.
     *
     * Integer indicates the maximum number of FAC cycles.  When this number is
     * 1, some computations and parallel communication can be avoided.
     *
     * Booleans that indicate whether the FAC preconditioner associated with
     * this concrete SAMRAI::solv::FACOperatorStrategy is performing
     * presmoothing or using a nonzero initial guess.
     *
     * If the SAMRAI::solv::FACPreconditioner is:
     *
     *     1) being used as a preconditioner for a Krylov subspace linear
     *        solver, and
     *     2) not performing presmoothing (which doesn't seem to work anyway)
     *
     * then we can save a substantial amount of work.
     */
    const SAMRAI::solv::FACPreconditioner<NDIM>* d_preconditioner;
    int d_fac_max_cycles;
    bool d_fac_uses_presmoothing, d_fac_initial_guess_nonzero;

    /*
     * In some cases, it may be desirable to skip restricting the solution or
     * the residual.
     */
    bool d_skip_restrict_sol, d_skip_restrict_residual;

    /*
     * Coarse level solver parameters.
     */
    std::string d_coarse_solver_choice;
    double d_coarse_solver_tol;
    int d_coarse_solver_max_its;

    //\}

    /*!
     * \name Internal context and scratch data.
     */
    //\{

    /*
     * Variable context for internally maintained hierarchy data.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;

    /*
     * Patch descriptor indices for scratch data.
     */
    int d_side_scratch_idx, d_cell_scratch_idx;

    //\}

    /*!
     * \name Boundary condition handling objects.
     */
    //\{

    SAMRAI::tbox::Pointer<IBTK::CartSideRobinPhysBdryOp> d_U_bc_op;
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* const d_default_U_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_U_bc_coefs;

    SAMRAI::tbox::Pointer<IBTK::CartCellRobinPhysBdryOp> d_P_bc_op;
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* const d_default_P_bc_coef;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef;

    bool d_homogeneous_bc;
    double d_current_time, d_new_time, d_dt;

    //\}

    /*!
     * \name Various refine and coarsen objects.
     */
    //\{

    /*
     * Variable fill pattern object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_U_op_stencil_fill_pattern,  d_P_op_stencil_fill_pattern, d_U_synch_fill_pattern;

    /*
     * Error prolongation (refinement) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_U_prolongation_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_P_prolongation_refine_operator;
    SAMRAI::tbox::Pointer<IBTK::CoarseFineBoundaryRefinePatchStrategy> d_U_cf_bdry_op;
    SAMRAI::tbox::Pointer<IBTK::CoarseFineBoundaryRefinePatchStrategy> d_P_cf_bdry_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_prolongation_refine_patch_strategy;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_prolongation_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_prolongation_refine_schedules;

    /*
     * Solution restriction (coarsening) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_U_urestriction_coarsen_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_P_urestriction_coarsen_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_urestriction_coarsen_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_urestriction_coarsen_schedules;

    /*
     * Residual restriction (coarsening) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_U_rrestriction_coarsen_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_P_rrestriction_coarsen_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_rrestriction_coarsen_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_rrestriction_coarsen_schedules;

    /*
     * Refine operator for side data from coarser level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_U_ghostfill_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_P_ghostfill_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_refine_schedules;

    /*
     * Refine operator for side data from same level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_U_ghostfill_nocoarse_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_P_ghostfill_nocoarse_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_nocoarse_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_nocoarse_refine_schedules;

    /*
     * Operator for side data synchronization on same level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_U_synch_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_U_synch_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_U_synch_refine_schedules;

    //\}
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/INSStaggeredBoxRelaxationFACOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredBoxRelaxationFACOperator
