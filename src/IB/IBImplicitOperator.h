// Filename: IBImplicitOperator.h
// Created on 26 Aug 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010 Boyce Griffith
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef included_IBImplicitOperator
#define included_IBImplicitOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petsc.h>

// IBAMR INCLUDES
#include <ibamr/INSStaggeredStokesOperator.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
class IBImplicitHierarchyIntegrator;

/*!
 * \brief Class IBImplicitOperator is a concrete IBTK::GeneralOperator which
 * implements an implicit staggered-grid (MAC) discretization of the IB method.
 *
 * This class is intended to be used with an iterative Newton-Krylov solver.
 *
 * \see IBImplicitHierarchyIntegrator
 */
class IBImplicitOperator
    : public IBTK::GeneralOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    IBImplicitOperator(
        SAMRAI::tbox::Pointer<INSStaggeredStokesOperator> stokes_op,
        IBImplicitHierarchyIntegrator* ib_implicit_integrator);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBImplicitOperator();

    /*!
     * \brief Set the current time interval.
     */
    void
    setTimeInterval(
        const double current_time,
        const double new_time);

    /*!
     * \name General operator functionality.
     */
    //\{

    /*!
     * \brief Compute \f$y=F[x]\f$.
     *
     * \note initializeOperatorState() must be called prior to any calls to
     * apply().
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=F[x]
     */
    virtual void
    apply(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

    /*!
     * \brief Compute hierarchy dependent data required for computing y=F[x] .
     *
     * \note Call deallocateOperatorState() to remove any data allocated by this
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
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeOperatorState().
     *
     * \see initializeOperatorState
     */
    virtual void
    deallocateOperatorState();

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
    IBImplicitOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBImplicitOperator(
        const IBImplicitOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBImplicitOperator&
    operator=(
        const IBImplicitOperator& that);

    // Whether the operator is initialized.
    bool d_is_initialized;

    // The simulation time.
    double d_current_time, d_new_time, d_dt;

    // The corresponding Stokes operator.
    SAMRAI::tbox::Pointer<INSStaggeredStokesOperator> d_stokes_op;

    // The IB implicit integrator provides most of the IB-related functionality.
    IBImplicitHierarchyIntegrator* d_ib_implicit_integrator;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBImplicitOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBImplicitOperator
