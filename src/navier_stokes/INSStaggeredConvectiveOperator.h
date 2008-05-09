#ifndef included_INSStaggeredConvectiveOperator
#define included_INSStaggeredConvectiveOperator

// Filename: INSStaggeredConvectiveOperator.h
// Last modified: <08.May.2008 18:24:33 griffith@box230.cims.nyu.edu>
// Created on 08 May 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/GeneralOperator.h>

// SAMRAI INCLUDES
#include <RefineAlgorithm.h>
#include <RefineOperator.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredConvectiveOperator is a concrete
 * IBTK::GeneralOperator which implements a upwind convective differencing
 * operator based on the piecewise parabolic method (PPM).
 *
 * \see INSStaggeredHierarchyIntegrator
 */
class INSStaggeredConvectiveOperator
    : public IBTK::GeneralOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    INSStaggeredConvectiveOperator(
        const double rho,
        const double mu,
        const double lambda,
        const bool conservation_form,
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > refine_alg,
        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_op,
        std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > refine_scheds)
        : d_is_initialized(false),
          d_U_scratch_idx(),
          d_rho(rho),
          d_mu(mu),
          d_lambda(lambda),
          d_conservation_form(conservation_form),
          d_refine_alg(refine_alg),
          d_refine_op(refine_op),
          d_refine_scheds(refine_scheds),
          d_hierarchy(NULL),
          d_coarsest_ln(-1),
          d_finest_ln(-1)
        {
            // intentionally blank
            return;
        }// INSStaggeredConvectiveOperator

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSStaggeredConvectiveOperator()
        {
            deallocateOperatorState();
            return;
        }// ~INSStaggeredConvectiveOperator

    /*!
     * \brief Set the patch data index corresponding to the scratch velocity.
     */
    void
    setVelocityScratchPatchDataIndex(
        const int U_scratch_idx)
        {
            d_U_scratch_idx = U_scratch_idx;
            return;
        }// setVelocityScratchPatchDataIndex

    /*!
     * \name General operator functionality.
     */
    //\{

    /*!
     * \brief Compute \f$y=F[x]\f$.
     *
     * Before calling apply(), the form of the vectors \a x and \a y should be
     * set properly by the user on all patch interiors on the specified range of
     * levels in the patch hierarchy.  The user is responsible for all data
     * management for the quantities associated with the vectors.  In
     * particular, patch data in these vectors must be allocated prior to
     * calling this method.
     *
     * \param x input vector
     * \param y output vector, i.e., \f$y=F[x]\f$
     *
     * <b>Conditions on Parameters:</b>
     * - vectors \a x and \a y must have same hierarchy
     * - vectors \a x and \a y must have same structure, depth, etc.
     *
     * In general, the vectors \a x and \a y \em cannot be the same.
     *
     * \see initializeOperatorState
     */
    virtual void
    apply(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

    /*!
     * \brief Compute hierarchy dependent data required for computing y=F[x] and
     * z=F[x]+y.
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

            // Get the hierarchy configuration.
            d_hierarchy = in.getPatchHierarchy();
            d_coarsest_ln = in.getCoarsestLevelNumber();
            d_finest_ln = in.getFinestLevelNumber();
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(d_hierarchy == out.getPatchHierarchy());
            TBOX_ASSERT(d_coarsest_ln == out.getCoarsestLevelNumber());
            TBOX_ASSERT(d_finest_ln == out.getFinestLevelNumber());
#endif
            // Allocate scratch data.
            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                if (!level->checkAllocated(d_U_scratch_idx))
                {
                    level->allocatePatchData(d_U_scratch_idx);
                }
            }
            d_is_initialized = true;
            return;
        }// initializeOperatorState

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeOperatorState().
     *
     * \note It is safe to call deallocateOperatorState() when the operator
     * state is already deallocated.
     *
     * \see initializeOperatorState
     */
    virtual void
    deallocateOperatorState()
        {
            if (!d_is_initialized) return;

            // Deallocate scratch data.
            for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
            {
                SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                if (level->checkAllocated(d_U_scratch_idx))
                {
                    level->deallocatePatchData(d_U_scratch_idx);
                }
            }
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
    INSStaggeredConvectiveOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredConvectiveOperator(
        const INSStaggeredConvectiveOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredConvectiveOperator&
    operator=(
        const INSStaggeredConvectiveOperator& that);

    // Whether the operator is initialized.
    bool d_is_initialized;

    // Scratch data.
    int d_U_scratch_idx;

    // Problem coefficients.
    const double d_rho;
    const double d_mu;
    const double d_lambda;

    // Whether to use conservative or non-conservative differencing.
    const bool d_conservation_form;

    // Data communication algorithms, operators, and schedules.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_refine_alg;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_refine_op;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_refine_scheds;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSStaggeredConvectiveOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredConvectiveOperator
