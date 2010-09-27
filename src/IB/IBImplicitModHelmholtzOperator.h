// Filename: IBImplicitModHelmholtzOperator.h
// Created on 30 Aug 2010 by Boyce Griffith
//
// Copyright (c) 2002-2010, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBImplicitModHelmholtzOperator
#define included_IBImplicitModHelmholtzOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petsc.h>

// IBTK INCLUDES
#include <ibtk/SCLaplaceOperator.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBImplicitModHelmholtzOperator is a concrete
 * IBTK::GeneralOperator which implements a linear operator arrising from an
 * implicit staggered-grid (MAC) discretization of the IB method.
 *
 * This class is intended to be used with an iterative Newton-Krylov solver.
 *
 * \see IBImplicitHierarchyIntegrator
 */
class IBImplicitModHelmholtzOperator
    : public IBTK::LinearOperator
{
public:
    /*!
     * \brief Class constructor.
     */
    IBImplicitModHelmholtzOperator(
        SAMRAI::tbox::Pointer<IBTK::SCLaplaceOperator> helmholtz_op,
        SAMRAI::tbox::Pointer<IBTK::LinearOperator> ib_SJSstar_op);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBImplicitModHelmholtzOperator();

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
     * The user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must have enough
     *   ghost cells for computation of Ax).
     *
     * In general, the vectors x and y \em cannot be the same.
     *
     * \note The operator MUST be initialized prior to calling apply.
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
     *
     * \note The default implementation is empty.
     */
    virtual void
    initializeOperatorState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& in,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeOperatorState().
     *
     * Remove all hierarchy dependent data set by initializeOperatorState().  It
     * is safe to call deallocateOperatorState() when the operator state is
     * already deallocated.
     *
     * \see initializeOperatorState
     *
     * \note The default implementation is empty.
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
    IBImplicitModHelmholtzOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBImplicitModHelmholtzOperator(
        const IBImplicitModHelmholtzOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBImplicitModHelmholtzOperator&
    operator=(
        const IBImplicitModHelmholtzOperator& that);

    // Whether the operator is initialized.
    bool d_is_initialized;

    // The Stokes operator.
    SAMRAI::tbox::Pointer<IBTK::SCLaplaceOperator> d_helmholtz_op;

    // The IB force operator.
    SAMRAI::tbox::Pointer<IBTK::LinearOperator> d_ib_SJSstar_op;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBImplicitModHelmholtzOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBImplicitModHelmholtzOperator
