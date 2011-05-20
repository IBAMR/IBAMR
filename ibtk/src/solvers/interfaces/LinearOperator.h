// Filename: LinearOperator.h
// Created on 14 Sep 2003 by Boyce Griffith
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

#ifndef included_LinearOperator
#define included_LinearOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/GeneralOperator.h>

// SAMRAI INCLUDES
#include <SAMRAIVectorReal.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class LinearOperator provides an abstract interface for the
 * specification of linear operators to compute \f$ y=Ax \f$ and \f$ z=Ax+y \f$
 * and, optionally, \f$ y=A^{T} x \f$ and \f$ z=A^{T}x+y \f$.
 */
class LinearOperator
    : public GeneralOperator
{
public:
    /*!
     * \brief Default constructor.
     *
     * \param is_symmetric  Boolean that indicates whether the operator is symmetric.
     */
    LinearOperator(
        bool is_symmetric=false);

    /*!
     * \brief Empty virtual destructor.
     */
    virtual
    ~LinearOperator();

    /*!
     * \name Linear operator functionality.
     */
    //\{

    /*!
     * \brief Indicates whether the linear operator is symmetric.
     */
    bool
    isSymmetric() const;

    /*!
     * \brief Modify y to account for inhomogeneous boundary conditions.
     *
     * Before calling this function, the form of the vector y should be set
     * properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in this vector should be allocated.
     * The user is responsible for managing the storage for the vectors.
     *
     * \note The operator MUST be initialized prior to calling
     * modifyRhsForInhomogeneousBc.
     *
     * \see initializeOperatorState
     *
     * \param y output: y=Ax
     *
     * \note A default implementation is provided which does nothing but warns
     * that inhomogeneous boundary conditions are not properly supported.
     */
    virtual void
    modifyRhsForInhomogeneousBc(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

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
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y) = 0;

    /*!
     * \brief Compute z=Ax+y.
     *
     * Before calling this function, the form of the vectors x, y, and z should
     * be set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * The user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must have enough
     *   ghost cells for computation of Ax).
     *
     * In general, the vectors x and z \em cannot be the same.
     *
     * \note The operator MUST be initialized prior to calling apply.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y input
     * \param z output: z=Ax+y
     *
     * \note A default implementation is provided which employs apply() and
     * SAMRAI::solv::SAMRAIVectorReal::add().
     */
    virtual void
    applyAdd(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& z);

    /*!
     * \brief Compute y=A'x.
     *
     * Before calling this function, the form of the vectors x and y should be
     * set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * The user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must have enough
     *   ghost cells for computation of A'x).
     *
     * In general, the vectors x and y \em cannot be the same.
     *
     * \note The operator MUST be initialized prior to calling apply.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y output: y=A'x
     *
     * \note A default implementation is provided.  If the operator is
     * symmetric, apply() is used.  Otherwise, the default implementation causes
     * an unrecoverable error to occur.
     */
    virtual void
    applyAdjoint(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

    /*!
     * \brief Compute z=A'x+y.
     *
     * Before calling this function, the form of the vectors x, y, and z should
     * be set properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in these vectors should be allocated.
     * The user is responsible for managing the storage for the vectors.
     *
     * Conditions on arguments:
     * - vectors must have same hierarchy
     * - vectors must have same variables (except that x \em must have enough
     *   ghost cells for computation of A'x).
     *
     * In general, the vectors x and z \em cannot be the same.
     *
     * \note The operator MUST be initialized prior to calling apply.
     *
     * \see initializeOperatorState
     *
     * \param x input
     * \param y input
     * \param z output: z=A'x+y
     *
     * \note A default implementation is provided which employs applyAdjoint()
     * and SAMRAI::solv::SAMRAIVectorReal::add().
     */
    virtual void
    applyAdjointAdd(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& z);

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
        bool enabled=true) = 0;

    //\}

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    LinearOperator(
        const LinearOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    LinearOperator&
    operator=(
        const LinearOperator& that);

    /*!
     * Indicates whether the linear operator is symmetric.
     */
    bool d_is_symmetric;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/LinearOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_LinearOperator
