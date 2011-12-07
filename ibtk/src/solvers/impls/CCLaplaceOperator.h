// Filename: CCLaplaceOperator.h
// Created on 19 Sep 2003 by Boyce Griffith
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

#ifndef included_CCLaplaceOperator
#define included_CCLaplaceOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/LinearOperator.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsReal.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PatchHierarchy.h>
#include <PoissonSpecifications.h>
#include <RobinBcCoefStrategy.h>
#include <SAMRAIVectorReal.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class CCLaplaceOperator is a concrete LinearOperator which implements
 * a globally second-order accurate cell-centered finite difference
 * discretization of a scalar elliptic operator of the form \f$ L = C I + \nabla
 * \cdot D \nabla\f$.
 */
class CCLaplaceOperator
    : public LinearOperator
{
public:
    /*!
     * \brief Constructor for class CCLaplaceOperator initializes the operator
     * coefficients and boundary conditions for a scalar-valued operator.
     *
     * \param object_name     String used to register internal variables and for error reporting purposes.
     * \param poisson_spec    Laplace operator coefficients.
     * \param bc_coef         Robin boundary conditions to use with this class.
     * \param homogeneous_bc  Whether to employ the homogeneous form of the boundary conditions.
     */
    CCLaplaceOperator(
        const std::string& object_name,
        const SAMRAI::solv::PoissonSpecifications& poisson_spec,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef,
        bool homogeneous_bc=true);

    /*!
     * \brief Constructor for class CCLaplaceOperator initializes the operator
     * coefficients and boundary conditions for a vector-valued operator.
     *
     * \param object_name     String used to register internal variables and for error reporting purposes.
     * \param poisson_spec    Laplace operator coefficients.
     * \param bc_coefs        Robin boundary conditions to use with this class.
     * \param homogeneous_bc  Whether to employ the homogeneous form of the boundary conditions.
     */
    CCLaplaceOperator(
        const std::string& object_name,
        const SAMRAI::solv::PoissonSpecifications& poisson_spec,
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs,
        bool homogeneous_bc=true);

    /*!
     * \brief Constructor for class CCLaplaceOperator initializes the operator
     * coefficients and boundary conditions for a vector-valued operator.
     *
     * \param object_name     String used to register internal variables and for error reporting purposes.
     * \param poisson_spec    Laplace operator coefficients.
     * \param bc_coefs        Robin boundary conditions to use with this class.
     * \param homogeneous_bc  Whether to employ the homogeneous form of the boundary conditions.
     */
    CCLaplaceOperator(
        const std::string& object_name,
        const SAMRAI::solv::PoissonSpecifications& poisson_spec,
        const blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs,
        bool homogeneous_bc=true);

    /*!
     * \brief Destructor.
     */
    ~CCLaplaceOperator();

    /*!
     * \brief Set the SAMRAI::solv::PoissonSpecifications object used to specify
     * the coefficients for the scalar-valued or vector-valued Laplace operator.
     */
    void
    setPoissonSpecifications(
        const SAMRAI::solv::PoissonSpecifications& poisson_spec);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy object used to specify
     * physical boundary conditions.
     *
     * \note \a bc_coef may be NULL.  In this case, homogeneous Dirichlet
     * boundary conditions are employed.
     *
     * \param bc_coef  Pointer to an object that can set the Robin boundary condition coefficients
     */
    void
    setPhysicalBcCoef(
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.
     *
     * \param bc_coefs  Vector of pointers to objects that can set the Robin boundary condition coefficients
     */
    void
    setPhysicalBcCoefs(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.
     *
     * \param bc_coefs  Vector of pointers to objects that can set the Robin boundary condition coefficients
     */
    void
    setPhysicalBcCoefs(
        const blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs);

    /*!
     * \brief Specify whether the boundary conditions are homogeneous.
     */
    void
    setHomogeneousBc(
        bool homogeneous_bc);

    /*!
     * \brief Set the hierarchy time, for use with the refinement schedules and
     * boundary condition routines employed by the object.
     */
    void
    setTime(
        double time);

    /*!
     * \brief Set the HierarchyMathOps object used by the operator.
     */
    void
    setHierarchyMathOps(
        SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops);

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
     * \note The operator MUST be initialized prior to calling
     * modifyRhsForInhomogeneousBc.
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
    CCLaplaceOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    CCLaplaceOperator(
        const CCLaplaceOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    CCLaplaceOperator&
    operator=(
        const CCLaplaceOperator& that);

    // Housekeeping.
    std::string d_object_name;

    // Operator parameters.
    bool d_is_initialized;
    int d_ncomp;
    double d_apply_time;

    // Cached communications operators.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_fill_pattern;
    std::vector<HierarchyGhostCellInterpolation::InterpolationTransactionComponent> d_transaction_comps;
    SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> d_hier_bdry_fill, d_no_fill;

    // Scratch data.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_x, d_b;

    // Problem specification.
    SAMRAI::solv::PoissonSpecifications d_poisson_spec;
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* const d_default_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
    bool d_homogeneous_bc, d_correcting_rhs;

    // Mathematical operators.
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM,double> > d_hier_cc_data_ops;
    SAMRAI::tbox::Pointer<HierarchyMathOps> d_hier_math_ops;
    bool d_hier_math_ops_external;

    // Hierarchy configuration.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/CCLaplaceOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_CCLaplaceOperator
