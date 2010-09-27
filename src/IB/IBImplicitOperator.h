// Filename: IBImplicitOperator.h
// Created on 26 Aug 2010 by Boyce Griffith
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

#ifndef included_IBImplicitOperator
#define included_IBImplicitOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petsc.h>

// IBAMR INCLUDES
#include <ibamr/IBImplicitSFROperator.h>
#include <ibamr/INSStaggeredStokesOperator.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
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
        SAMRAI::tbox::Pointer<IBImplicitSFROperator> ib_SFR_op);

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

    // The Stokes operator.
    SAMRAI::tbox::Pointer<INSStaggeredStokesOperator> d_stokes_op;

    // The IB force operator.
    SAMRAI::tbox::Pointer<IBImplicitSFROperator> d_ib_SFR_op;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/IBImplicitOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBImplicitOperator
