// Filename: SCLaplaceOperator.h
// Created on 11 Apr 2008 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_SCLaplaceOperator
#define included_SCLaplaceOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "IntVector.h"
#include "PatchHierarchy.h"
#include "SAMRAIVectorReal.h"
#include "VariableFillPattern.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/LaplaceOperator.h"
#include "tbox/Pointer.h"

namespace IBTK
{
class StaggeredPhysicalBoundaryHelper;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class SCLaplaceOperator is a concrete LaplaceOperator which implements
 * a globally second-order accurate side-centered finite difference
 * discretization of a scalar elliptic operator of the form \f$ L = C I + \nabla
 * \cdot D \nabla\f$.
 */
class SCLaplaceOperator : public LaplaceOperator
{
public:
    /*!
     * \brief Constructor for class SCLaplaceOperator initializes the operator
     * coefficients and boundary conditions to default values.
     */
    SCLaplaceOperator(const std::string& object_name, bool homogeneous_bc = true);

    /*!
     * \brief Destructor.
     */
    ~SCLaplaceOperator();

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
    void apply(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& y);

    /*!
     * \brief Compute hierarchy-dependent data required for computing y=Ax (and
     * y=A'x).
     *
     * \param in input vector
     * \param out output vector
     *
     * \see KrylovLinearSolver::initializeSolverState
     */
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& in,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& out);

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
    void deallocateOperatorState();

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SCLaplaceOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SCLaplaceOperator(const SCLaplaceOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SCLaplaceOperator& operator=(const SCLaplaceOperator& that);

    // Operator parameters.
    int d_ncomp;

    // Cached communications operators.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_fill_pattern;
    std::vector<HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_transaction_comps;
    SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> d_hier_bdry_fill, d_no_fill;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_x, d_b;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    // Dirichlet boundary condition utilities.
    std::vector<SAMRAI::tbox::Pointer<StaggeredPhysicalBoundaryHelper> > d_bc_helpers;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_SCLaplaceOperator
