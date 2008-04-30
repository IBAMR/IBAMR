#ifndef included_INSStaggeredProjectionPreconditioner
#define included_INSStaggeredProjectionPreconditioner

// Filename: INSStaggeredProjectionPreconditioner.h
// Last modified: <29.Apr.2008 15:37:04 griffith@box230.cims.nyu.edu>
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
        const std::string projection_type,
        const double rho,
        const double mu,
        const double lambda,
        const double current_time,
        const double new_time,
        const bool normalize_pressure,
        SAMRAI::tbox::Pointer<IBTK::LinearSolver> helmholtz_solver,
        SAMRAI::tbox::Pointer<HierarchyProjector> hier_projector,
        SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > hier_cc_data_ops,
        SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops,
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> P_bdry_fill_op)
        : d_projection_type(projection_type),
          d_rho(rho),
          d_mu(mu),
          d_lambda(lambda),
          d_current_time(current_time),
          d_new_time(new_time),
          d_dt(new_time-current_time),
          d_normalize_pressure(normalize_pressure),
          d_helmholtz_solver(helmholtz_solver),
          d_hier_projector(hier_projector),
          d_hier_cc_data_ops(hier_cc_data_ops),
          d_hier_sc_data_ops(hier_sc_data_ops),
          d_hier_math_ops(hier_math_ops),
          d_P_bdry_fill_op(P_bdry_fill_op),
          d_no_fill_op(SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>(NULL))
        {
            // Get the control volume weight variables and patch data descriptor
            // indices.
            d_wgt_cc_var = d_hier_math_ops->getCellWeightVariable();
            d_wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

            d_wgt_sc_var = d_hier_math_ops->getSideWeightVariable();
            d_wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();

            // Get the volume of the physical domain.
            d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();
            return;
        }// INSStaggeredProjectionPreconditioner

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSStaggeredProjectionPreconditioner()
        {
            // intentionally blank
            return;
        }// ~INSStaggeredProjectionPreconditioner

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Compute the action of the preconditioner.
     */
    virtual bool
    solveSystem(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x);

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
            // intentionally blank
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

    // The type of projection to perform ("pressure_increment" or
    // "pressure_update").
    const std::string d_projection_type;

    // Problem coefficients.
    const double d_rho;
    const double d_mu;
    const double d_lambda;

    // The simulation time.
    const double d_current_time, d_new_time, d_dt;

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
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSStaggeredProjectionPreconditioner.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_HierarchyProjector
