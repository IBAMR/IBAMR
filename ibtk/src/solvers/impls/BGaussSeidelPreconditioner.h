// Filename: BGaussSeidelPreconditioner.h
// Created on 11 Apr 2005 by Boyce Griffith
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

#ifndef included_BGaussSeidelPreconditioner
#define included_BGaussSeidelPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LinearOperator.h>
#include <ibtk/LinearSolver.h>

// SAMRAI INCLUDES
#include <SAMRAIVectorReal.h>
#include <tbox/ConstPointer.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <map>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class BGaussSeidelPreconditioner is a block Gauss-Seidel
 * preconditioner which implements the abstract LinearSolver interface.
 *
 * This solver class performs a single block Gauss-Seidel sweep, applying
 * specified component LinearSolver and LinearOperator objects to the components
 * of a supplied SAMRAI::solv::SAMRAIVectorReal vector, and doing so in a
 * multiplicative fashion.  Note that the block Gauss-Seidel algorithm is not
 * generally convergent, but can be used as a preconditioner for a
 * KrylovLinearSolver.
 *
 * Note that the default block Gauss-Seidel algorithm is not a symmetric linear
 * operator, even if the individual component linear operators and solvers are
 * symmetric.  Instead, the algorithm applies the component preconditioners to
 * the vector components starting with the first vector component and ending
 * with the last vector component.  The algorithm can be symmetrized via the
 * setSymmetricPreconditioner() member function, and the order in which vector
 * components are visited can be reversed via the setReversedOrder() member
 * function.
 *
 * \note Class BJacobiPreconditioner implements the additive (i.e., block
 * Jacobi) version of this algorithm.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 symmetric_preconditioner = FALSE   // see setSymmetricPreconditioner()
 reverse_order = FALSE              // see setReversedOrder()
 initial_guess_nonzero = FALSE      // see setInitialGuessNonzero()
 rel_residual_tol = 1.0e-6          // see setRelativeTolerance()
 abs_residual_tol = 1.0e-30         // see setAbsoluteTolerance()
 max_iterations = 1                 // see setMaxIterations()
 \endverbatim
 *
 */
class BGaussSeidelPreconditioner
    : public virtual LinearSolver
{
public:
    /*!
     * \brief Default constructor.
     *
     * \param input_db optional SAMRAI::tbox::Database for input.
     */
    BGaussSeidelPreconditioner(
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db=NULL);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~BGaussSeidelPreconditioner();

    /*!
     * \brief Set the preconditioner to be employed on the specified vector
     * component.
     */
    void
    setComponentPreconditioner(
        SAMRAI::tbox::Pointer<LinearSolver> preconditioner,
        const unsigned int component);

    /*!
     * \brief Set the linear operators to be employed on the specified vector
     * component.
     */
    void
    setComponentOperators(
        std::vector<SAMRAI::tbox::Pointer<LinearOperator> > linear_ops,
        const unsigned int component);

    /*!
     * \brief Indicate whether to apply the component preconditioners
     * symmetrically.
     */
    void
    setSymmetricPreconditioner(
        const bool symmetric_preconditioner);

    /*!
     * \brief Indicate whether to apply the component preconditioners in
     * reversed order (i.e., starting with the last component and ending with
     * the first component).
     */
    void
    setReversedOrder(
        const bool reverse_order);

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Solve the linear system of equations \f$Ax=b\f$ for \f$x\f$.
     *
     * Before calling solveSystem(), the form of the solution \a x and
     * right-hand-side \a b vectors must be set properly by the user on all
     * patch interiors on the specified range of levels in the patch hierarchy.
     * The user is responsible for all data management for the quantities
     * associated with the solution and right-hand-side vectors.  In particular,
     * patch data in these vectors must be allocated prior to calling this
     * method.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     *
     * <b>Conditions on Parameters:</b>
     * - vectors \a x and \a b must have same patch hierarchy
     * - vectors \a x and \a b must have same structure, depth, etc.
     *
     * \note The vector arguments for solveSystem() need not match those for
     * initializeSolverState().  However, there must be a certain degree of
     * similarity, including:\par
     * - hierarchy configuration (hierarchy pointer and range of levels)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the solution \a x and right-hand-side \a b
     *   vectors
     *
     * \note The solver need not be initialized prior to calling solveSystem();
     * however, see initializeSolverState() and deallocateSolverState() for
     * opportunities to save overhead when performing multiple consecutive
     * solves.
     *
     * \see initializeSolverState
     * \see deallocateSolverState
     *
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise
     */
    virtual bool
    solveSystem(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     *
     * By default, the solveSystem() method computes some required hierarchy
     * dependent data before solving and removes that data after the solve.  For
     * multiple solves that use the same hierarchy configuration, it is more
     * efficient to:
     *
     * -# initialize the hierarchy-dependent data required by the solver via
     *    initializeSolverState(),
     * -# solve the system one or more times via solveSystem(), and
     * -# remove the hierarchy-dependent data via deallocateSolverState().
     *
     * Note that it is generally necessary to reinitialize the solver state when
     * the hierarchy configuration changes.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     *
     * <b>Conditions on Parameters:</b>
     * - vectors \a x and \a b must have same patch hierarchy
     * - vectors \a x and \a b must have same structure, depth, etc.
     *
     * \note The vector arguments for solveSystem() need not match those for
     * initializeSolverState().  However, there must be a certain degree of
     * similarity, including:\par
     * - hierarchy configuration (hierarchy pointer and range of levels)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the solution \a x and right-hand-side \a b
     *   vectors
     *
     * \note It is safe to call initializeSolverState() when the state is
     * already initialized.  In this case, the solver state is first deallocated
     * and then reinitialized.
     *
     * \see deallocateSolverState
     */
    virtual void
    initializeSolverState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * \note It is safe to call deallocateSolverState() when the solver state is
     * already deallocated.
     *
     * \see initializeSolverState
     */
    virtual void
    deallocateSolverState();

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    virtual void
    setInitialGuessNonzero(
        bool initial_guess_nonzero=true);

    /*!
     * \brief Get whether the initial guess is non-zero.
     */
    virtual bool
    getInitialGuessNonzero() const;

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    virtual void
    setMaxIterations(
        int max_iterations);

    /*!
     * \brief Get the maximum number of iterations to use per solve.
     */
    virtual int
    getMaxIterations() const;

    /*!
     * \brief Set the absolute residual tolerance for convergence.
     */
    virtual void
    setAbsoluteTolerance(
        double abs_residual_tol);

    /*!
     * \brief Get the absolute residual tolerance for convergence.
     */
    virtual double
    getAbsoluteTolerance() const;

    /*!
     * \brief Set the relative residual tolerance for convergence.
     */
    virtual void
    setRelativeTolerance(
        double rel_residual_tol);

    /*!
     * \brief Get the relative residual tolerance for convergence.
     */
    virtual double
    getRelativeTolerance() const;

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent linear solve.
     */
    virtual int
    getNumIterations() const;

    /*!
     * \brief Return the residual norm from the most recent iteration.
     */
    virtual double
    getResidualNorm() const;

    //\}

    /*!
     * \name Logging functions.
     */
    //\{

    /*!
     * \brief Enable or disable logging.
     */
    virtual void
    enableLogging(
        bool enabled=true);

    //\}

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    BGaussSeidelPreconditioner(
        const BGaussSeidelPreconditioner& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    BGaussSeidelPreconditioner&
    operator=(
        const BGaussSeidelPreconditioner& that);

    /*!
     * \brief Extract the individual components of a
     * SAMRAI::solv::SAMRAIVectorReal object, and create individual
     * SAMRAI::solv::SAMRAIVectorReal objects to correspond to each of the
     * components.
     */
    static std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >
    getComponentVectors(
        const SAMRAI::tbox::ConstPointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > x);

    /*!
     * Boolean value to indicate whether the preconditioner is presently
     * initialized.
     */
    bool d_is_initialized;

    /*!
     * The component preconditioners.
     */
    std::map<unsigned int,SAMRAI::tbox::Pointer<LinearSolver> > d_pc_map;

    /*!
     * The component operators.
     */
    std::map<unsigned int,std::vector<SAMRAI::tbox::Pointer<LinearOperator> > > d_linear_ops_map;

    /*!
     * Parameters to specify the ordering of the application of the component
     * preconditioners.
     */
    bool d_symmetric_preconditioner, d_reverse_order;

    /*!
     * Solver configuration parameters.
     */
    bool d_initial_guess_nonzero;
    double d_rel_residual_tol;
    double d_abs_residual_tol;
    int d_max_iterations;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/BGaussSeidelPreconditioner.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_BGaussSeidelPreconditioner
