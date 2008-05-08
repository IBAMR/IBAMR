#ifndef included_INSStaggeredLinearOperator
#define included_INSStaggeredLinearOperator

// Filename: INSStaggeredLinearOperator.h
// Last modified: <29.Apr.2008 14:27:57 griffith@box230.cims.nyu.edu>
// Created on 29 Mar 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LinearOperator.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredLinearOperator is a concrete IBTK::LinearOperator
 * which implements a staggered grid (MAC) discretization of the incompressible
 * Stokes operator.
 *
 * This class is intended to be used with an iterative (Krylov or Newton-Krylov)
 * incompressible flow solver.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class INSStaggeredLinearOperator
    : public IBTK::LinearOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    INSStaggeredLinearOperator(
        const double rho,
        const double mu,
        const double lambda,
        const double current_time,
        const double new_time,
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> U_P_bdry_fill_op)
        : d_rho(rho),
          d_mu(mu),
          d_lambda(lambda),
          d_current_time(current_time),
          d_new_time(new_time),
          d_dt(d_new_time-d_current_time),
          d_hier_math_ops(hier_math_ops),
          d_U_P_bdry_fill_op(U_P_bdry_fill_op),
          d_no_fill_op(SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>(NULL))
        {
            // intentionally blank
            return;
        }// INSStaggeredLinearOperator


    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSStaggeredLinearOperator()
        {
            // intentionally blank
            return;
        }// ~INSStaggeredLinearOperator

    /*!
     * \name Linear operator functionality.
     */
    //\{

    /*!
     * \brief Compute y=Ax.
     *
     * Before calling this function, the form of the vectors x and y should be
     * set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * Thus, the user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must
     * have enough ghost cells for computation of Ax).
     *
     * \note In general, the vectors x and y \em cannot be the same.
     *
     * Upon return from this function, the y vector will contain the result of
     * the application of A to x.
     *
     * initializeOperatorState must be called prior to any calls to
     * applyOperator.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=Ax
     */
    virtual void
    apply(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

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
    INSStaggeredLinearOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredLinearOperator(
        const INSStaggeredLinearOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredLinearOperator&
    operator=(
        const INSStaggeredLinearOperator& that);

    // Problem coefficients.
    const double d_rho;
    const double d_mu;
    const double d_lambda;

    // The simulation time.
    const double d_current_time, d_new_time, d_dt;

    // Math objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;

    // Boundary condition objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_U_P_bdry_fill_op, d_no_fill_op;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSStaggeredLinearOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_HierarchyProjector
