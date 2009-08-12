#ifndef included_INSStaggeredBlockFactorizationPreconditioner
#define included_INSStaggeredBlockFactorizationPreconditioner

// Filename: INSStaggeredBlockFactorizationPreconditioner.h
// Last modified: <12.Aug.2009 18:21:11 griffith@boyce-griffiths-mac-pro.local>
// Created on 22 Sep 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/INSCoefs.h>

// IBTK INCLUDES
#include <ibtk/LinearSolver.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredBlockFactorizationPreconditioner is a concrete
 * IBTK::LinearSolver which implements a staggered grid (MAC) block
 * factorization (approximate Schur complement) preconditioner for the
 * incompressible Navier-Stokes equations.
 *
 * This class is intended to be used with an iterative (Krylov or Newton-Krylov)
 * incompressible flow solver.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class INSStaggeredBlockFactorizationPreconditioner
    : public IBTK::LinearSolver
{
public:
    /*!
     * \brief Class constructor
     */
    INSStaggeredBlockFactorizationPreconditioner(
        const INSCoefs& problem_coefs,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef,
        const bool normalize_pressure,
        SAMRAI::tbox::Pointer<IBTK::LinearSolver> velocity_helmholtz_solver,
        SAMRAI::tbox::Pointer<IBTK::LinearSolver> pressure_poisson_solver,
        SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > hier_cc_data_ops,
        SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops,
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSStaggeredBlockFactorizationPreconditioner();

    /*!
     * \brief Set the current time interval.
     */
    void
    setTimeInterval(
        const double current_time,
        const double new_time);

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Compute the action of the preconditioner.
     */
    virtual bool
    solveSystem(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     *
     * <b>Conditions on Parameters:</b>
     * - vectors \a x and \a b must have same patch hierarchy
     * - vectors \a x and \a b must have same structure, depth, etc.
     *
     * \note It is safe to call initializeSolverState() when the solver state is
     * already initialized.
     *
     * \see deallocateSolverState
     *
     * \note A default implementation is provided which does nothing.
     */
    virtual void
    initializeSolverState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * \note It is safe to call deallocateSolverState() when the solver state is
     * already deallocated.
     *
     * \see initializeSolverState
     *
     * \note A default implementation is provided which does nothing.
     */
    virtual void
    deallocateSolverState();

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    virtual void
    setInitialGuessNonzero(
        bool initial_guess_nonzero=true);

    /*!
     * \brief Get whether the initial guess is non-zero.
     */
    virtual bool
    getInitialGuessNonzero() const;

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    virtual void
    setMaxIterations(
        int max_iterations);

    /*!
     * \brief Get the maximum number of iterations to use per solve.
     */
    virtual int
    getMaxIterations() const;

    /*!
     * \brief Set the absolute residual tolerance for convergence.
     */
    virtual void
    setAbsoluteTolerance(
        double abs_residual_tol);

    /*!
     * \brief Get the absolute residual tolerance for convergence.
     */
    virtual double
    getAbsoluteTolerance() const;

    /*!
     * \brief Set the relative residual tolerance for convergence.
     */
    virtual void
    setRelativeTolerance(
        double rel_residual_tol);

    /*!
     * \brief Get the relative residual tolerance for convergence.
     */
    virtual double
    getRelativeTolerance() const;

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent linear solve.
     */
    virtual int
    getNumIterations() const;

    /*!
     * \brief Return the residual norm from the most recent iteration.
     */
    virtual double
    getResidualNorm() const;

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
    virtual void
    enableLogging(
        bool enabled=true);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredBlockFactorizationPreconditioner();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredBlockFactorizationPreconditioner(
        const INSStaggeredBlockFactorizationPreconditioner& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredBlockFactorizationPreconditioner&
    operator=(
        const INSStaggeredBlockFactorizationPreconditioner& that);

    // Indicates whether the preconditioner should output logging messages.
    bool d_do_log;

    // Whether the operator is initialized.
    bool d_is_initialized;

    // The simulation time.
    double d_current_time, d_new_time, d_dt;

    // Problem coefficients.
    const INSCoefs& d_problem_coefs;
    SAMRAI::solv::PoissonSpecifications d_pressure_helmholtz_spec;

    // Normalize the pressure when necessary.
    const bool d_normalize_pressure;

    // Linear solver functionality.
    SAMRAI::tbox::Pointer<IBTK::LinearSolver> d_velocity_helmholtz_solver;
    SAMRAI::tbox::Pointer<IBTK::LinearSolver> d_pressure_poisson_solver;

    // Math objects.
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_wgt_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_wgt_sc_var;
    int d_wgt_cc_idx, d_wgt_sc_idx;
    double d_volume;

    // Boundary condition objects.
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_P_bdry_fill_op, d_no_fill_op;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_U_var;
    int d_U_scratch_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_P_var;
    int d_P_scratch_idx;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSStaggeredBlockFactorizationPreconditioner.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredBlockFactorizationPreconditioner
