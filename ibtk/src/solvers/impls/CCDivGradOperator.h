// Filename: CCDivGradOperator.h
// Created on 18 Feb 2010 by Boyce Griffith
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

#ifndef included_CCDivGradOperator
#define included_CCDivGradOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LinearOperator.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsReal.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CCDivGradOperator is a concrete LinearOperator which implements
 * a cell-centered finite difference operator which corresponds to
 * \f$-\nabla^{\mathrm cc} \cdot \nabla^{\mathrm cc}\f$, i.e., minus the
 * composition of the cell-centered divergence and cell-centered gradient
 * operators.
 *
 * \note The class implementation is intended to be used only for uniform grids
 * with periodic boundaries.  This class may also work in some cases for
 * locally-refined discretizations, but if it does, its success is purely
 * accidental.
 */
class CCDivGradOperator
    : public LinearOperator
{
public:
    /*!
     * \brief Constructor for class CCDivGradOperator.
     *
     * \param object_name  String used to register internal variables and for error reporting purposes.
     */
    CCDivGradOperator(
        const std::string& object_name);

    /*!
     * \brief Destructor.
     */
    ~CCDivGradOperator();

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
    void
    apply(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& y);

    /*!
     * \brief Compute hierarchy-dependent data required for computing y=Ax (and
     * y=A'x).
     *
     * \param in input vector
     * \param out output vector
     *
     * \see KrylovLinearSolver::initializeSolverState
     */
    void
    initializeOperatorState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& in,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& out);

    /*!
     * \brief Remove all hierarchy-dependent data computed by
     * initializeOperatorState().
     *
     * Remove all hierarchy-dependent data set by initializeOperatorState().  It
     * is safe to call deallocateOperatorState() even if the state is already
     * deallocated.
     *
     * \see initializeOperatorState
     * \see KrylovLinearSolver::deallocateSolverState
     */
    void
    deallocateOperatorState();

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    CCDivGradOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CCDivGradOperator(
        const CCDivGradOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CCDivGradOperator&
    operator=(
        const CCDivGradOperator& that);

    // Cached communications operators.
    SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> d_scalar_hier_bdry_fill, d_vector_hier_bdry_fill;

    // Scratch data.
    int d_scalar_scratch_idx, d_vector_scratch_idx;
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > d_scalar_scratch_var, d_vector_scratch_var;

    // Mathematical operators.
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM,double> > d_hier_cc_data_ops;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/CCDivGradOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CCDivGradOperator
