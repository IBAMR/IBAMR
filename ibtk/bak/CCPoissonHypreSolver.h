#ifndef included_CCPoissonHypreSolver
#define included_CCPoissonHypreSolver

// Filename: CCPoissonHypreSolver.h
// Last modified: <07.Sep.2006 23:32:43 boyce@bigboy.nyconnect.com>
// Created on 20 Jan 2006 by Boyce Griffith (boyce@boyce.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// HYPRE INCLUDES
#ifndef included_HYPRE_sstruct_ls
#define included_HYPRE_sstruct_ls
#include <HYPRE_sstruct_ls.h>
#endif

// IBTK INCLUDES
#include <ibtk/LinearSolver.h>

// SAMRAI INCLUDES
#include <ArrayData.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <CellData.h>
#include <CoarseFineBoundary.h>
#include <GhostCellRobinBcCoefs.h>
#include <OutersideVariable.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PoissonSpecifications.h>
#include <RobinBcCoefStrategy.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SimpleCellRobinBcCoefs.h>
#include <Variable.h>
#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <ostream>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Concrete LinearSolver for solving elliptic equations of the
 * form \f$ \mbox{$L u$} = \mbox{$(C I + \nabla \cdot D \nabla) u$} =
 * f \f$ on a range of levels in a SAMRAI::hier::PatchHierarchy object
 * using \em hypre.
 *
 * Class CCPoissonHypreSolver uses the \em hypre library to solve
 * linear equations of the form \f$ \nabla \cdot D \nabla u + C u = f
 * \f$, where \f$C\f$ is a cell-centered array, \f$D\f$ is a
 * side-centered array, and \f$u\f$ and \f$f\f$ are cell-centered
 * arrays (see class SAMRAI::solv::PoissonSpecifications).  The
 * discretization is globally second-order accurate away from
 * coarse-fine interfaces and first-order accurate in the vicinity of
 * coarse-fine interfaces.
 *
 * The user must perform the following steps to use class
 * CCPoissonHypreLevelSolver:
 * - Create a CCPoissonHypreLevelSolver object.
 * - Set the problem specification via setPoissonSpecifications(),
 *   setPhysicalBcCoefObject(), and setHomogeneousBc().
 * - Initialize CCPoissonHypreLevelSolver object using the function
 *   initializeSolverState().
 * - Solve the linear system, passing in
 *   SAMRAI::solv::SAMRAIVectorReal objects corresponding to \f$u\f$
 *   and \f$f\f$.
 *
 * Sample parameters for initialization from database (and their
 * default values):
 * \verbatim
 *     enable_logging = FALSE       // Whether to print some data for debugging
 *     max_iterations = 10          // Max iterations used by Hypre
 *     residual_tol = 1.0e-8        // Residual tolerance used by Hypre
 *     num_pre_relax_steps = 1      // # of presmoothing steps used by Hypre
 *     num_post_relax_steps = 1     // # of postsmoothing steps used by Hypre
 * \endverbatim
 * \todo Update the above values.
 *
 * \note The interface for class CCPoissonHypreSolver is identical to
 * the less general CCPoissonHypreLevelSolver.  Use of the level
 * solver is strongly encouranged for problems involving only a single
 * level in the patch hierarchy, whereas class CCPoissonHypreSolver
 * must be used for problems which involve data on multiple levels of
 * the patch hierarchy.
 */
class CCPoissonHypreSolver
    : public LinearSolver
{
public:
    /*!
     * \brief Constructor.
     *
     * \param object_name name of object.
     * \param input_db optional input Database for input.
     */
    CCPoissonHypreSolver(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db=NULL);

    /*!
     * The destructor releases all internally managed data.
     */
    virtual ~CCPoissonHypreSolver();

    /*!
     * \name Functions for specifying the Poisson problem.
     */
    //\{

    /*!
     * \brief Set the scalar Poisson equation specifications.
     *
     * \note
     * - Changes to the solver do not take effect until the next time
     *   that the solver is (re-)initialized.
     * - The SAMRAI::solv::PoissonSpecifications pointer may be NULL.
     *   In this case, the diffusion coefficient D is set to -1.0 and
     *   the damping parameter C is set to 0.0, so that the solver is
     *   configured to solve the Poisson problem, \f$-\Delta u = f\f$.
     */
    void setPoissonSpecifications(
        const SAMRAI::solv::PoissonSpecifications* poisson_spec);

    /*!
     * \brief Specify boundary condition through the use of a Robin
     * boundary condition object.
     *
     * The Robin boundary condition object is used when setting the
     * matrix coefficients and when solving the system.  If your
     * boundary conditions are fixed values at ghost cell centers, use
     * the SAMRAI::solv::GhostCellRobinBcCoefs implementation of the
     * SAMRAI::solv::RobinBcCoefStrategy strategy.
     *
     * \param bc_coef pointer to an object that can set the Robin
     *        boundary condition coefficients.
     * \param bc_var SAMRAI::hier::Variable pointer to be passed
     *        to SAMRAI::solv::RobinBcCoefStrategy::setBcCoefs(), but
     *        otherwise unused by this class.
     */
    void setPhysicalBcCoefObject(
        const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > bc_var=NULL);

    /*!
     * \brief Specify whether the boundary conditions are homogeneous.
     */
    void setHomogeneousBc(
        bool homogeneous_bc);

    /*!
     * \brief Set the depth of the solution data used by the solver.
     *
     * Changing the depth after setting up the matrix is permissible,
     * as the solution data does not affect the matrix.
     */
    void setSolDataDepth(
        const int depth);

    /*!
     * \brief Set the depth of the rhs data used by the solver.
     *
     * Changing the depth after setting up the matrix is permissible,
     * as the rhs data does not affect the matrix.
     */
    void setRhsDataDepth(
        const int depth);

    /*!
     * \brief Specify the levels that need to be reset the next time
     * the solver is re-initialized.
     *
     * When the solver is initialized, only the specified range of
     * levels are reset in the solver state the next time that the
     * solver is initialized.  If the solver is not initialized, this
     * method has no effect.
     *
     * Set \a coarsest_ln and \a finest_ln to -1 to reset all levels
     * in the patch hierarchy, use coarsest_ln = finest_ln = -1.
     *
     * \note This method may be used to avoid some unnecessary
     * computations when the hierarchy is regridded.  The range of
     * levels specified must include all levels which need to be reset
     * by
     * SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration().
     */
    void setResetLevels(
        const int coarsest_ln,
        const int finest_ln);

    //\}

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Attempt to solve the linear system Ax=b for x.
     *
     * The return value is true if the solver converges to the
     * specified tolerances and false otherwise.
     *
     * Before calling this function, the form of the solution and
     * right-hand-side quantities should be set properly by the user
     * on all patch interiors on the range of levels covered by the
     * solver.  All data in these vectors should be allocated.  The
     * user is responsible for managing the storage for the solution
     * and right-hand-side.
     *
     * Conditions on arguments:
     * - vectors x and b must have same hierarchy
     * - vectors x and b must have same variables
     *
     * The initial value of the solution vector x is used as the
     * initial guess to the solution.  Upon return from this function,
     * x will contain the result of the solve.
     *
     * See initializeSolverState() and deallocateSolverState() for
     * opportunities to save overhead when using multiple consecutive
     * solves.
     *
     * NOTE: The solver need not be initialized prior to calling
     * solveSystem().
     *
     * \see initializeSolverState
     *
     * \param x solution vector x
     * \param b right-hand-side vector b
     *
     * \return indicates whether the solver converged to the
     * specified tolerances
     */
    virtual bool solveSystem(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Compute hierarchy dependent data required for solving
     * Ax=b.
     *
     * By default, the solveSystem() method computes some required
     * hierarchy dependent data before solving and removes that data
     * after the solve.  For multiple solves using the same hierarchy
     * configuration, it is more efficient to use
     * initializeSolverState() to manually initialize the solver and
     * deallocateSolverState() to remove the hierarchy dependent data.
     *
     * If solveSystem() detects that the solver state is already
     * initialized, it will \em NOT change the state.
     *
     * The vector arguments for solveSystem() need not match those for
     * initializeSolverState().  However, there must be a certain
     * degree of similarity, including
     * - hierarchy configuration (hierarchy pointer and level range)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the solution x and
     *   right-hand-side b vectors
     *
     * NOTE: It is generally necessary to reinitialize the solver
     * state when the hierarchy configuration changes.
     *
     * It is safe to call initializeSolverState() when the state is
     * already initialized.  In this case, the solver state is first
     * deallocated and then reinitialized.
     *
     * Conditions on arguments:
     * - vectors x and b must have same hierarchy
     * - vectors x and b must have same structure, depth, etc.
     *
     * Call deallocateSolverState() to remove any data allocated
     * by this method.
     *
     * \see deallocateSolverState
     *
     * \param x solution vector x
     * \param b right-hand-side vector b
     */
    virtual void initializeSolverState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * Remove all hierarchy dependent data set by
     * initializeSolverState().  It is safe to call
     * deallocateSolverState() when the solver state is already
     * deallocated.
     *
     * \see initializeSolverState
     */
    virtual void deallocateSolverState();

    //\}

    /*!
     * \name Functions to access solving parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    virtual void setInitialGuessNonzero(
        bool initial_guess_nonzero=true);

    /*!
     * \brief Get whether the initial guess is non-zero.
     */
    virtual bool getInitialGuessNonzero() const;

    /*!
     * \brief Set the max number of iterations to use per solve.
     */
    virtual void setMaxIterations(
        int max_iterations);

    /*!
     * \brief Get the max number of iterations to use per solve.
     */
    virtual int getMaxIterations() const;

    /*!
     * \brief Set the absolute residual tolerance for stopping.
     *
     * NOTE: This member function presently is a no-op, since there is
     * no way to specify an absolute stopping criteria to most of the
     * hypre solvers.
     */
    virtual void setAbsoluteTolerance(
        double abs_residual_tol);

    /*!
     * \brief Get the absolute residual tolerance for stopping.
     */
    virtual double getAbsoluteTolerance() const;

    /*!
     * \brief Set the relative residual tolerance for stopping.
     */
    virtual void setRelativeTolerance(
        double rel_residual_tol);

    /*!
     * \brief Get the relative residual tolerance for stopping.
     */
    virtual double getRelativeTolerance() const;

    //\}

    /*!
     * \name Functions to access data on last solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent linear
     * solve.
     */
    virtual int getNumIterations() const;

    /*!
     * \brief Return residual norm from the just-completed
     * iteration.
     *
     * The latest computed norm is the one returned.
     */
    virtual double getResidualNorm() const;

    //\}

    /*!
     * \name Logging functions.
     */
    //\{

    /*!
     * \brief Enable or disable logging.
     *
     * \param enabled logging state: true=on, false=off
     */
    virtual void enableLogging(
        bool enabled=true);

    /*!
     * \brief Print data members for debugging.
     *
     * NOTE: The default implementation is empty.
     */
    virtual void printClassData(
        std::ostream& os) const;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     */
    CCPoissonHypreSolver();

    /*!
     * \brief Copy constructor.
     *
     * NOTE: This constructor is not implemented and should not be
     * used.
     *
     * \param from The value to copy to this object.
     */
    CCPoissonHypreSolver(
        const CCPoissonHypreSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * NOTE: This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CCPoissonHypreSolver& operator=(
        const CCPoissonHypreSolver& that);

    /*!
     * \brief Trim a boundary box so that it does not stick out past a
     * patch domain in direction transverse to the boundary normal.
     *
     * The supplied boundary box must be of type 1 (see
     * SAMRAI::hier::BoundaryBox::getBoundaryType()).
     */
    SAMRAI::hier::BoundaryBox<NDIM> trimBoundaryBox(
        const SAMRAI::hier::BoundaryBox<NDIM>& boundary_box,
        const SAMRAI::hier::Patch<NDIM>& patch) const;

    /*!
     * \brief Return box describing the index space of surfaces
     * defined by a boundary box.
     *
     * The supplied boundary_box must be type 1 (see
     * SAMRAI::hier::BoundaryBox::getBoundaryType()).
     */
    SAMRAI::hier::Box<NDIM> makeSideBoundaryBox(
        const SAMRAI::hier::BoundaryBox<NDIM>& boundary_box) const;

    /*!
     * \brief Modify the rhs to include boundary conditions.
     */
    void modifyBoundaryRHSEntries(
        const SAMRAI::hier::Patch<NDIM>& patch,
        const SAMRAI::tbox::Array<SAMRAI::hier::BoundaryBox<NDIM> >& bdry_boxes,
        const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
        SAMRAI::pdat::CellData<NDIM,double>& rhs);

    /*!
     * \brief Adjust matrix entries to account for boundary
     * conditions.
     */
    void adjustBoundaryMatrixEntries(
        SAMRAI::pdat::CellData<NDIM,double>& diagonal,
        SAMRAI::pdat::ArrayData<NDIM,double>& Ak0_data,
        const SAMRAI::pdat::SideData<NDIM,double>& off_diagonal,
        const SAMRAI::hier::Box<NDIM>& patch_box,
        const SAMRAI::pdat::ArrayData<NDIM,double>& acoef_data,
        const SAMRAI::hier::Box<NDIM>& bccoef_box,
        const SAMRAI::hier::BoundaryBox<NDIM>& trimmed_boundary_box,
        const double dx[NDIM]);

    /*!
     * \brief Compute diagonal matrix entries when C is non-constant.
     */
    void computeDiagonalEntries(
        SAMRAI::pdat::CellData<NDIM,double>& diagonal,
        const SAMRAI::pdat::CellData<NDIM,double>& C_data,
        const SAMRAI::pdat::SideData<NDIM,double>& off_diagonal,
        const SAMRAI::hier::Box<NDIM>& patch_box);

    /*!
     * \brief Compute diagonal matrix entries when C is constant.
     */
    void computeDiagonalEntries(
        SAMRAI::pdat::CellData<NDIM,double>& diagonal,
        const double C,
        const SAMRAI::pdat::SideData<NDIM,double>& off_diagonal,
        const SAMRAI::hier::Box<NDIM>& patch_box);

    /*!
     * \brief Compute diagonal matrix entries when C is zero.
     */
    void computeDiagonalEntries(
        SAMRAI::pdat::CellData<NDIM,double>& diagonal,
        const SAMRAI::pdat::SideData<NDIM,double>& off_diagonal,
        const SAMRAI::hier::Box<NDIM>& patch_box);

    void setupHypreSolver();
    void allocateHypreData();
    void setMatrixCoefficients();
    int solveSystem(
        const int u,
        const int f);
    void deallocateHypreData();
    void copyToHypre(
        HYPRE_SStructVector vector,
        SAMRAI::pdat::CellData<NDIM,double>& src,
        const int depth,
        const SAMRAI::hier::Box<NDIM>& box,
        const int part);
    void copyToHypre(
        HYPRE_SStructVector vector,
        const double scale,
        SAMRAI::pdat::CellData<NDIM,double>& src,
        const int depth,
        const SAMRAI::hier::Box<NDIM>& box,
        const int part);
    void copyFromHypre(
        SAMRAI::pdat::CellData<NDIM,double>& dst,
        const int depth,
        HYPRE_SStructVector vector,
        const SAMRAI::hier::Box<NDIM>& box,
        const int part);
    void copyFromHypre(
        const double scale,
        SAMRAI::pdat::CellData<NDIM,double>& dst,
        const int depth,
        HYPRE_SStructVector vector,
        const SAMRAI::hier::Box<NDIM>& box,
        const int part);

    /*
     * The object name is used for error reporting purposes.
     *
     * The boolean indicates whether this object has been initialized.
     */
    std::string d_object_name;
    bool d_is_initialized;

    /*!
     * \name Internal context and scratch data.
     */
    //\{

    /*
     * Variable context for internally maintained hierarchy data.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;

    /*
     * Patch descriptor indices.
     */
    const int d_is_refined_idx, d_scratch_idx;

    //\}

    /*!
     * \name Hierarchy-dependent objects.
     */
    //\{

    /*
     * Reference patch hierarchy and range of levels involved in the
     * solve.
     *
     * This variable is non-null between the initializeSolverState()
     * and deallocateSolverState() calls.  It is not truly needed,
     * because the hierarchy is obtainable through variables in most
     * function argument lists.  We use it to enforce working on one
     * hierarchy at a time.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    /*
     * Description of coarse-fine boundaries.
     *
     * There is one coarse-fine boundary object for each level.
     * d_coarse_fine_boundary[i] is the description of the coarse-fine
     * boundary between level i and level i-1.  The coarse-fine
     * boundary does not exist at the coarsest level, although the
     * SAMRAI::hier::CoarseFineBoundary object still exists (it should
     * not contain any boxes).
     *
     * This array is initialized in initializeSolverState() and
     * deallocated in deallocateSolverState().  When allocated, it is
     * allocated for the index range [0,d_coarsest_ln], though the
     * range [0,d_coarsest_ln-1] is not used.  This is okay because
     * SAMRAI::SAMRAI::hier::CoarseFineBoundary is a lightweight
     * object before it is set for a level.
     */
    std::vector<SAMRAI::hier::CoarseFineBoundary<NDIM> > d_cf_boundary;

    /*
     * Range of levels to be reset the next time the solver is
     * initialized.
     */
    bool d_in_initialize_solver_state;
    int d_coarsest_reset_ln, d_finest_reset_ln;

    //\}

    /*
     * Scalar Poisson equations specifications.
     */
    SAMRAI::solv::PoissonSpecifications d_poisson_spec;

    /*
     * Depths of the solution and rhs variables.
     */
    int d_sol_depth, d_rhs_depth;

    /*!
     * \name Hypre objects.
     */
    //\{
    HYPRE_SStructGrid    d_grid;
    HYPRE_SStructGraph   d_graph;
    HYPRE_SStructStencil d_stencil;
    HYPRE_SStructMatrix  d_matrix;
    HYPRE_SStructVector  d_rhs_vec, d_sol_vec;
    HYPRE_SStructSolver  d_solver, d_precond;

    int d_nparts;
    int *d_plevels;
    int (*d_rfactors)[3];

    std::string d_solver_type, d_precond_type;
    int d_max_iterations;
    double d_rel_residual_tol;
    bool d_initial_guess_nonzero;
    int d_rel_change;
    int d_num_pre_relax_steps, d_num_post_relax_steps;
    int d_memory_use;
    int d_relax_type;
    int d_skip_relax;
    int d_coarse_solver_type;
    int d_two_norm;

    int d_current_its;
    double d_current_residual_norm;
    //\}

    /*
     * Solution restriction (coarsening) operator.
     */
    bool d_synch_soln;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_urestriction_coarsen_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_urestriction_coarsen_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_urestriction_coarsen_schedules;

    /*!
     * \name Variables for debugging and analysis.
     */
    //\{

    /*
     * Flag to print solver info
     *
     * See setPrintSolverInfo().
     */
    bool d_print_solver_info, d_enable_logging;

    //\}
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include "CCPoissonHypreSolver.I"

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CCPoissonHypreSolver
