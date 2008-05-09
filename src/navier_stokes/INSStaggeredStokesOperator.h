#ifndef included_INSStaggeredStokesOperator
#define included_INSStaggeredStokesOperator

// Filename: INSStaggeredStokesOperator.h
// Last modified: <08.May.2008 18:24:56 griffith@box230.cims.nyu.edu>
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
 * \brief Class INSStaggeredStokesOperator is a concrete IBTK::LinearOperator
 * which implements a staggered grid (MAC) discretization of the incompressible
 * Stokes operator.
 *
 * This class is intended to be used with an iterative (Krylov or Newton-Krylov)
 * incompressible flow solver.
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class INSStaggeredStokesOperator
    : public IBTK::LinearOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    INSStaggeredStokesOperator(
        const double rho,
        const double mu,
        const double lambda,
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
        SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> U_P_bdry_fill_op)
        : d_is_initialized(false),
          d_current_time(std::numeric_limits<double>::quiet_NaN()),
          d_new_time(std::numeric_limits<double>::quiet_NaN()),
          d_dt(std::numeric_limits<double>::quiet_NaN()),
          d_rho(rho),
          d_mu(mu),
          d_lambda(lambda),
          d_helmholtz_spec("INSStaggeredStokesOperator::helmholtz_spec"),
          d_hier_math_ops(hier_math_ops),
          d_U_P_bdry_fill_op(U_P_bdry_fill_op),
          d_no_fill_op(SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>(NULL)),
          d_x_scratch(NULL)
        {
            // intentionally blank
            return;
        }// INSStaggeredStokesOperator

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSStaggeredStokesOperator()
        {
            deallocateOperatorState();
            return;
        }// ~INSStaggeredStokesOperator

    /*!
     * \brief Set the current time interval.
     */
    void
    setCurrentTimeInterval(
        const double current_time,
        const double new_time)
        {
            d_current_time = current_time;
            d_new_time = new_time;
            d_dt = d_new_time-d_current_time;
            d_helmholtz_spec.setCConstant((d_rho/d_dt)+0.5*d_lambda);
            d_helmholtz_spec.setDConstant(            -0.5*d_mu    );
            return;
        }// setCurrentTimeInterval

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
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y)
        {
            // Initialize the operator (if necessary).
            const bool deallocate_at_completion = !d_is_initialized;
            if (!d_is_initialized) initializeOperatorState(x,y);

            // Get the vector components.
            const int U_out_idx      =            y.getComponentDescriptorIndex(0);
            const int P_out_idx      =            y.getComponentDescriptorIndex(1);
            const int U_scratch_idx  = d_x_scratch->getComponentDescriptorIndex(0);
            const int P_scratch_idx  = d_x_scratch->getComponentDescriptorIndex(1);

            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_out_var = y.getComponentVariable(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_out_var = y.getComponentVariable(1);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_out_sc_var = U_out_var;
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_scratch_var = d_x_scratch->getComponentVariable(0);
            const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_scratch_var = d_x_scratch->getComponentVariable(1);

            SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_scratch_sc_var = U_scratch_var;
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_scratch_cc_var = P_scratch_var;

            d_x_scratch->copyVector(SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >(&x,false));

            // Type of coarsening to perform prior to setting coarse-fine
            // boundary and physical boundary ghost cell values.
            static const std::string DATA_COARSEN_TYPE = "CONSERVATIVE_COARSEN";

            // Type of extrapolation to use at physical boundaries.
            static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

            // Whether to enforce consistent interpolated values at Type 2
            // coarse-fine interface ghost cells.
            static const bool CONSISTENT_TYPE_2_BDRY = false;

            // Reset the interpolation operators and fill the data.
            typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
            InterpolationTransactionComponent U_component(U_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);  // XXXX
            InterpolationTransactionComponent P_component(P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, NULL);

            std::vector<InterpolationTransactionComponent> transaction_comps(2);
            transaction_comps[0] = U_component;
            transaction_comps[1] = P_component;

            d_U_P_bdry_fill_op->resetTransactionComponents(transaction_comps);
            d_U_P_bdry_fill_op->fillData(d_new_time);

            // Compute the action of the operator:
            //      A*[u;p] = [((rho/dt)*I-0.5*mu*L)*u + grad p; -div u].
            bool cf_bdry_synch;
            cf_bdry_synch = true;
            d_hier_math_ops->grad(
                U_out_idx, U_out_sc_var,
                cf_bdry_synch,
                1.0, P_scratch_idx, P_scratch_cc_var, d_no_fill_op, d_new_time);
            cf_bdry_synch = false;
            d_hier_math_ops->div(
                P_out_idx, P_out_cc_var,
                -1.0, U_scratch_idx, U_scratch_sc_var, d_no_fill_op, d_new_time,
                cf_bdry_synch);
            d_hier_math_ops->laplace(
                U_out_idx, U_out_sc_var,
                d_helmholtz_spec,
                U_scratch_idx, U_scratch_sc_var,
                d_no_fill_op, d_new_time,
                1.0,
                U_out_idx, U_out_sc_var);

            // Deallocate the operator (if necessary).
            if (deallocate_at_completion) deallocateOperatorState();
            return;
        }// apply

    /*!
     * \brief Compute hierarchy dependent data required for computing y=Ax and
     * z=Ax+y.
     *
     * The vector arguments for apply(), applyAdjoint(), etc, need not match
     * those for initializeOperatorState().  However, there must be a certain
     * degree of similarity, including
     * - hierarchy configuration (hierarchy pointer and level range)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the input and output vectors
     *
     * \note It is generally necessary to reinitialize the operator state when
     * the hierarchy configuration changes.
     *
     * It is safe to call initializeOperatorState() when the state is already
     * initialized.  In this case, the operator state is first deallocated and
     * then reinitialized.
     *
     * Conditions on arguments:
     * - input and output vectors must have same hierarchy
     * - input and output vectors must have same structure, depth, etc.
     *
     * Call deallocateOperatorState() to remove any data allocated by this
     * method.
     *
     * \see deallocateOperatorState
     *
     * \param in input vector
     * \param out output vector
     */
    virtual void
    initializeOperatorState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& in,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out)
        {
            if (d_is_initialized) deallocateOperatorState();

            d_x_scratch = in.cloneVector("INSStaggeredStokesOperator::x_scratch");
            d_x_scratch->allocateVectorData();

            d_is_initialized = true;
            return;
        }// initializeOperatorState

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeOperatorState().
     *
     * Remove all hierarchy dependent data set by initializeOperatorState().  It
     * is safe to call deallocateOperatorState() when the operator state is
     * already deallocated.
     *
     * \see initializeOperatorState
     */
    virtual void
    deallocateOperatorState()
        {
            if (!d_is_initialized) return;

            d_x_scratch->deallocateVectorData();
            d_x_scratch->freeVectorComponents();
            d_x_scratch.setNull();

            d_is_initialized = false;
            return;
        }// deallocateOperatorState

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
    INSStaggeredStokesOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredStokesOperator(
        const INSStaggeredStokesOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredStokesOperator&
    operator=(
        const INSStaggeredStokesOperator& that);

    // Whether the operator is initialized.
    bool d_is_initialized;

    // The simulation time.
    double d_current_time, d_new_time, d_dt;

    // Problem coefficients.
    const double d_rho;
    const double d_mu;
    const double d_lambda;
    SAMRAI::solv::PoissonSpecifications d_helmholtz_spec;

    // Math objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;

    // Boundary condition objects.
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_U_P_bdry_fill_op, d_no_fill_op;

    // Scratch data objects.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_x_scratch;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSStaggeredStokesOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredStokesOperator
