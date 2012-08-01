// Filename: INSStaggeredStokesOperator.h
// Created on 29 Mar 2008 by Boyce Griffith
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

#ifndef included_INSStaggeredStokesOperator
#define included_INSStaggeredStokesOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBAMR INCLUDES
#include <ibamr/INSProblemCoefs.h>
#include <ibamr/INSStaggeredPhysicalBoundaryHelper.h>
#include <ibamr/ibamr_enums.h>

// IBTK INCLUDES
#include <ibtk/LinearOperator.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredStokesOperator is a concrete IBTK::LinearOperator
 * which implements a staggered-grid (MAC) discretization of the incompressible
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
        const std::string& object_name,
        const INSProblemCoefs* problem_coefs,
        TimeSteppingType viscous_time_stepping_type,
        const blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM>& U_bc_coefs,
        SAMRAI::tbox::Pointer<INSStaggeredPhysicalBoundaryHelper> U_bc_helper,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef,
        SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops=SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps>(),
        bool homogeneous_bc=true);

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredStokesOperator();

    /*!
     * \brief Implementation of the apply method which supports either
     * homogeneous or inhomogeneous boundary conditions.
     */
    void
    apply(
        bool homogeneous_bc,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

    /*!
     * \name Linear operator functionality.
     */
    //\{

    /*!
     * \brief Modify y to account for inhomogeneous boundary conditions.
     *
     * Before calling this function, the form of the vector y should be set
     * properly by the user on all patch interiors on the range of levels
     * covered by the operator.  All data in this vector should be allocated.
     * The user is responsible for managing the storage for the vectors.
     *
     * \see initializeOperatorState
     *
     * \param y output: y=Ax
     */
    void
    modifyRhsForInhomogeneousBc(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

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
    void
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
     */
    void
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
     */
    void
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
    void
    enableLogging(
        bool enabled=true);

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

    // Housekeeping.
    std::string d_object_name;

    // Operator parameters.
    bool d_is_initialized;

    // Cached communications operators.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_U_fill_pattern, d_P_fill_pattern;
    std::vector<IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_transaction_comps;
    SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> d_hier_bdry_fill, d_no_fill;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_x, d_b;

    // Problem specification and mathematical operators.
    const INSProblemCoefs* d_problem_coefs;
    const TimeSteppingType d_viscous_time_stepping_type;
    blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM> d_U_bc_coefs;
    SAMRAI::tbox::Pointer<INSStaggeredPhysicalBoundaryHelper> d_U_bc_helper;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const d_P_bc_coef;
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> d_hier_math_ops;
    const bool d_hier_math_ops_external;
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibamr/INSStaggeredStokesOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredStokesOperator
