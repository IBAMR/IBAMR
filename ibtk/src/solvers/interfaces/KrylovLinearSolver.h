// Filename: KrylovLinearSolver.h
// Created on 08 Sep 2003 by Boyce Griffith
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

#ifndef included_KrylovLinearSolver
#define included_KrylovLinearSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/LinearOperator.h>
#include <ibtk/LinearSolver.h>

// SAMRAI INCLUDES
#include <SAMRAIVectorReal.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <ostream>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class KrylovLinearSolver provides an abstract interface for the
 * implementation of Krylov subspace solvers for linear problems of the form
 * \f$Ax=b\f$.
 */
class KrylovLinearSolver
    : public LinearSolver
{
public:
    /*!
     * \brief Empty constructor.
     */
    KrylovLinearSolver();

    /*!
     * \brief Empty virtual destructor.
     */
    virtual
    ~KrylovLinearSolver();

    /*!
     * \name Krylov solver functionality.
     */
    //\{

    /*!
     * \brief Set the linear operator used when solving \f$Ax=b\f$.
     */
    virtual void
    setOperator(
        SAMRAI::tbox::Pointer<LinearOperator> A) = 0;

    /*!
     * \brief Retrieve the linear operator used when solving \f$Ax=b\f$.
     */
    virtual SAMRAI::tbox::Pointer<LinearOperator>
    getOperator() const = 0;

    /*!
     * \brief Set the preconditioner used by the Krylov subspace method when
     * solving \f$Ax=b\f$.
     *
     * \note If the preconditioner is NULL, no preconditioning is performed.
     */
    virtual void
    setPreconditioner(
        SAMRAI::tbox::Pointer<LinearSolver> pc_solver=NULL) = 0;

    /*!
     * \brief Retrieve the preconditioner used by the Krylov subspace method
     * when solving \f$Ax=b\f$.
     */
    virtual SAMRAI::tbox::Pointer<LinearSolver>
    getPreconditioner() const = 0;

    /*!
     * \brief Set the nullspace of the linear system.
     *
     * Should not assume the basis vector, if any, to be normalized.  If the
     * basis vector is not normalized, the solver may normalize it.
     *
     * \note A default implementation is provided which calls the vector version
     * of setNullspace().
     */
    virtual void
    setNullspace(
        const bool contains_constant_vector,
        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > nullspace_basis_vec);

    /*!
     * \brief Set the nullspace of the linear system.
     *
     * Implementations can require the nullspace basis vectors to be orthogonal
     * but should not assume the basis vectors to be orthonormal.  If the basis
     * vectors are not orthonormal, the solver may normalize them.
     *
     * \note A default implementation is provided which does nothing.
     */
    virtual void
    setNullspace(
        const bool contains_constant_vector,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >& nullspace_basis_vecs);

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
     * \note Subclasses must be implemented so that the vector arguments for
     * solveSystem() need not match those for initializeSolverState().  However,
     * they are allowed to require a certain degree of similarity,
     * including:\par
     * - hierarchy configuration (hierarchy pointer and range of levels)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the solution \a x and right-hand-side \a b
     *   vectors
     *
     * \note Subclasses are required to be implemented so that the solver does
     * not need to be initialized prior to calling solveSystem(); however, see
     * initializeSolverState() and deallocateSolverState() for opportunities to
     * save overhead when performing multiple consecutive solves.
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
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b) = 0;

    /*!
     * \brief Compute hierarchy dependent data required for solving \f$Ax=b\f$.
     *
     * In a typical implementation, the solveSystem() method will compute some
     * required hierarchy dependent data before the solve, and then remove that
     * data after the solve.  For multiple solves that use the same hierarchy
     * configuration, it is more efficient to:
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
     * \note Subclasses must be implemented so that the vector arguments for
     * solveSystem() need not match those for initializeSolverState().  However,
     * they are allowed to require a certain degree of similarity,
     * including:\par
     * - hierarchy configuration (hierarchy pointer and range of levels)
     * - number, type and alignment of vector component data
     * - ghost cell widths of data in the solution \a x and right-hand-side \a b
     *   vectors
     *
     * \note Subclasses are required to be implemented so that it is safe to
     * call initializeSolverState() when the solver state is already
     * initialized.  In this case, the solver state should be first deallocated
     * and then reinitialized.
     *
     * \note Subclasses are required to be implemented so that when linear
     * operator or preconditioner objects have been registered with the solver
     * via setOperator() and setPreconditioner(), they are also initialized by
     * initializeSolverState().
     *
     * \see deallocateSolverState
     *
     * \note A default implementation is provided which does nothing.
     */
    virtual void
    initializeSolverState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * \note Subclasses are required to be implemented so that it is safe to
     * call deallocateSolverState() when the solver state is already
     * deallocated.
     *
     * \note Subclasses are required to be implemented so that when linear
     * operator or preconditioner objects have been registered with the solver
     * via setOperator() and setPreconditioner(), they are also deallocated by
     * deallocateSolverState().
     *
     * \see initializeSolverState
     *
     * \note A default implementation is provided which does nothing.
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
        bool initial_guess_nonzero=true) = 0;

    /*!
     * \brief Get whether the initial guess is non-zero.
     */
    virtual bool
    getInitialGuessNonzero() const = 0;

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    virtual void
    setMaxIterations(
        int max_iterations) = 0;

    /*!
     * \brief Get the maximum number of iterations to use per solve.
     */
    virtual int
    getMaxIterations() const = 0;

    /*!
     * \brief Set the absolute residual tolerance for convergence.
     */
    virtual void
    setAbsoluteTolerance(
        double abs_residual_tol) = 0;

    /*!
     * \brief Get the absolute residual tolerance for convergence.
     */
    virtual double
    getAbsoluteTolerance() const = 0;

    /*!
     * \brief Set the relative residual tolerance for convergence.
     */
    virtual void
    setRelativeTolerance(
        double rel_residual_tol) = 0;

    /*!
     * \brief Get the relative residual tolerance for convergence.
     */
    virtual double
    getRelativeTolerance() const = 0;

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent linear solve.
     */
    virtual int
    getNumIterations() const = 0;

    /*!
     * \brief Return the residual norm from the most recent iteration.
     */
    virtual double
    getResidualNorm() const = 0;

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
    KrylovLinearSolver(
        const KrylovLinearSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    KrylovLinearSolver&
    operator=(
        const KrylovLinearSolver& that);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/KrylovLinearSolver.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_KrylovLinearSolver
