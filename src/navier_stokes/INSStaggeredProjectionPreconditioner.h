#ifndef included_INSStaggeredProjectionPreconditioner
#define included_INSStaggeredProjectionPreconditioner

// Filename: INSStaggeredProjectionPreconditioner.h
// Last modified: <11.May.2008 18:56:22 griffith@box230.cims.nyu.edu>
// Created on 29 Mar 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/HierarchyProjector.h>

// IBTK INCLUDES
#include <ibtk/LinearSolver.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredProjectionPreconditioner is a concrete
 * IBTK::LinearSolver which implements a staggered grid (MAC) projection solver
 * for the incompressible Stokes operator.
 *
 * This class is intended to be used with an iterative (Krylov or Newton-Krylov)
 * incompressible flow solver.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class INSStaggeredProjectionPreconditioner
    : public IBTK::LinearSolver
{
public:
    /*!
     * \brief Class constructor
     */
    INSStaggeredProjectionPreconditioner(
        const double rho,
        const double mu,
        const double lambda,
        const std::string projection_type,
        const bool normalize_pressure,
        SAMRAI::tbox::Pointer<IBTK::LinearSolver> helmholtz_solver,
        SAMRAI::tbox::Pointer<HierarchyProjector> hier_projector,
        SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > hier_cc_data_ops,
        SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops,
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops)
        : d_do_log(false),
          d_is_initialized(false),
          d_current_time(std::numeric_limits<double>::quiet_NaN()),
          d_new_time(std::numeric_limits<double>::quiet_NaN()),
          d_dt(std::numeric_limits<double>::quiet_NaN()),
          d_rho(rho),
          d_mu(mu),
          d_lambda(lambda),
          d_pressure_helmholtz_spec("INSStaggeredProjectionPreconditioner::pressure_helmholtz_spec"),
          d_projection_type(projection_type),
          d_normalize_pressure(normalize_pressure),
          d_helmholtz_solver(helmholtz_solver),
          d_hier_projector(hier_projector),
          d_hier_cc_data_ops(hier_cc_data_ops),
          d_hier_sc_data_ops(hier_sc_data_ops),
          d_hier_math_ops(hier_math_ops),
          d_wgt_cc_var(NULL),
          d_wgt_sc_var(NULL),
          d_wgt_cc_idx(-1),
          d_wgt_sc_idx(-1),
          d_volume(std::numeric_limits<double>::quiet_NaN()),
          d_P_bdry_fill_op(SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>(NULL)),
          d_no_fill_op(SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>(NULL)),
          d_x_scratch(NULL),
          d_b_scratch(NULL),
          d_hierarchy(NULL),
          d_coarsest_ln(-1),
          d_finest_ln(-1)
        {
            // intentionally blank
            return;
        }// INSStaggeredProjectionPreconditioner

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSStaggeredProjectionPreconditioner()
        {
            deallocateSolverState();
            return;
        }// ~INSStaggeredProjectionPreconditioner

    /*!
     * \brief Set the current time interval.
     */
    void
    setTimeInterval(
        const double current_time,
        const double new_time)
        {
            d_current_time = current_time;
            d_new_time = new_time;
            d_dt = d_new_time-d_current_time;
            d_pressure_helmholtz_spec.setCConstant(1.0+0.5*d_dt*d_lambda/d_rho);
            d_pressure_helmholtz_spec.setDConstant(   -0.5*d_dt*d_mu    /d_rho);
            return;
        }// setTimeInterval

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
        bool initial_guess_nonzero=true)
        {
            // intentionally blank
            return;
        }// setInitialGuessNonzero


    /*!
     * \brief Get whether the initial guess is non-zero.
     */
    virtual bool
    getInitialGuessNonzero() const
        {
            // intentionally blank
            return true;
        }// getInitialGuessNonzero

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    virtual void
    setMaxIterations(
        int max_iterations)
        {
            // intentionally blank
            return;
        }// setMaxIterations

    /*!
     * \brief Get the maximum number of iterations to use per solve.
     */
    virtual int
    getMaxIterations() const
        {
            // intentionally blank
            return 1;
        }// getMaxIterations

    /*!
     * \brief Set the absolute residual tolerance for convergence.
     */
    virtual void
    setAbsoluteTolerance(
        double abs_residual_tol)
        {
            // intentionally blank
            return;
        }// setAbsoluteTolerance

    /*!
     * \brief Get the absolute residual tolerance for convergence.
     */
    virtual double
    getAbsoluteTolerance() const
        {
            // intentionally blank
            return 0.0;
        }// getAbsoluteTolerance

    /*!
     * \brief Set the relative residual tolerance for convergence.
     */
    virtual void
    setRelativeTolerance(
        double rel_residual_tol)
        {
            // intentionally blank
            return;
        }// setRelativeTolerance

    /*!
     * \brief Get the relative residual tolerance for convergence.
     */
    virtual double
    getRelativeTolerance() const
        {
            // intentionally blank
            return 0.0;
        }// getRelativeTolerance

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent linear solve.
     */
    virtual int
    getNumIterations() const
        {
            // intentionally blank
            return 0;
        }// getNumIterations

    /*!
     * \brief Return the residual norm from the most recent iteration.
     */
    virtual double
    getResidualNorm() const
        {
            return 0.0;
        }// getResidualNorm

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
        bool enabled=true)
        {
            d_do_log = enabled;
            return;
        }// enableLogging

    /*!
     * \brief Print out internal class data for debugging.
     */
    virtual void
    printClassData(
        std::ostream& os) const
        {
            // intentionally blank
            return;
        }// printClassData

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredProjectionPreconditioner();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredProjectionPreconditioner(
        const INSStaggeredProjectionPreconditioner& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredProjectionPreconditioner&
    operator=(
        const INSStaggeredProjectionPreconditioner& that);

    // Indicates whether the preconditioner should output logging messages.
    bool d_do_log;

    // Whether the operator is initialized.
    bool d_is_initialized;

    // The simulation time.
    double d_current_time, d_new_time, d_dt;

    // Problem coefficients.
    const double d_rho;
    const double d_mu;
    const double d_lambda;
    SAMRAI::solv::PoissonSpecifications d_pressure_helmholtz_spec;

    // The type of projection to perform ("pressure_increment" or
    // "pressure_update").
    const std::string d_projection_type;

    // Normalize the pressure when necessary.
    const bool d_normalize_pressure;

    // Helmholtz solver functionality.
    SAMRAI::tbox::Pointer<IBTK::LinearSolver> d_helmholtz_solver;

    // Projection functionality.
    SAMRAI::tbox::Pointer<HierarchyProjector> d_hier_projector;

    // Math objects.
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > d_hier_sc_data_ops;
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > d_wgt_cc_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > d_wgt_sc_var;
    int d_wgt_cc_idx, d_wgt_sc_idx;
    double d_volume;

    // Boundary condition objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_P_bdry_fill_op, d_no_fill_op;

    // Scratch data objects.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_x_scratch, d_b_scratch;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSStaggeredProjectionPreconditioner.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredProjectionPreconditioner
