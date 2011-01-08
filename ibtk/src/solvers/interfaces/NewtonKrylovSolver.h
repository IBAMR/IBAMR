// Filename: NewtonKrylovSolver.h
// Created on 18 Nov 2003 by Boyce Griffith
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

#ifndef included_NewtonKrylovSolver
#define included_NewtonKrylovSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/GeneralOperator.h>
#include <ibtk/JacobianOperator.h>
#include <ibtk/KrylovLinearSolver.h>

// SAMRAI INCLUDES
#include <SAMRAIVectorReal.h>
#include <tbox/DescribedClass.h>
#include <tbox/Pointer.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class NewtonKrylovSolver provides an abstract interface for the
 * implementation of inexact Newton-Krylov solvers for nonlinear problems of the
 * form \f$ F[x]=b \f$.
 */
class NewtonKrylovSolver
    : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Empty constructor.
     */
    NewtonKrylovSolver();

    /*!
     * \brief Empty virtual destructor.
     */
    virtual
    ~NewtonKrylovSolver();

    /*!
     * \name Newton-Krylov solver functionality.
     */
    //\{

    /*!
     * \brief Set the nonlinear operator \f$F[x]\f$ used by the solver.
     */
    virtual void
    setOperator(
        SAMRAI::tbox::Pointer<GeneralOperator> op) = 0;

    /*!
     * \brief Retrieve the nonlinear operator \f$F[x]\f$ used by the solver.
     */
    virtual SAMRAI::tbox::Pointer<GeneralOperator>
    getOperator() const = 0;

    /*!
     * \brief Return the vector in which the approximate solution is stored.
     *
     * \note Implementations of this member function are permitted to return a
     * NULL pointer if the solver is not initialized.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >
    getSolutionVector() const = 0;

    /*!
     * \brief Return the vector in which the nonlinear function evaluation is
     * stored.
     *
     * \note Implementations of this member function are permitted to return a
     * NULL pointer if the solver is not initialized.
     */
    virtual SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >
    getFunctionVector() const = 0;

    /*!
     * \brief Set the Jacobian operator \f$J[x] = F'[x]\f$ used by the solver.
     *
     * \note Subclasses should be implemented so that if a Jacobian object is
     * not explicitly provided to the solver, a Jacobian-free inexact
     * Newton-Krylov method is employed to approximate the action of the
     * Jacobian.
     */
    virtual void
    setJacobian(
        SAMRAI::tbox::Pointer<JacobianOperator> J) = 0;

    /*!
     * \brief Retrieve the Jacobian operator \f$J[x] = F'[x]\f$ used by the
     * solver.
     */
    virtual SAMRAI::tbox::Pointer<JacobianOperator>
    getJacobian() const = 0;

    /*!
     * \brief Retrieve the Krylov linear solver used in computing Newton step
     * directions.
     */
    virtual SAMRAI::tbox::Pointer<KrylovLinearSolver>
    getLinearSolver() const = 0;

    /*!
     * \brief Solve the system \f$F[x]=b\f$ for \f$x\f$.
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
     * \brief Compute hierarchy dependent data required for solving
     * \f$F[x]=b\f$.
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
     * \note Subclasses are required to be implemented so that when any operator
     * objects have been registered with the solver via setOperator() or
     * setJacobian(), they are also initialized by initializeSolverState().
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
     * \note Subclasses are required to be implemented so that when any operator
     * objects have been registered with the solver via setOperator() or
     * setJacobian(), they are also deallocated by deallocateSolverState().
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
     * \brief Set the maximum number of nonlinear iterations to use per solve.
     */
    virtual void
    setMaxIterations(
        int max_iterations) = 0;

    /*!
     * \brief Get the maximum number of nonlinear iterations to use per solve.
     */
    virtual int
    getMaxIterations() const = 0;

    /*!
     * \brief Set the maximum number of function evaluations to use per solve.
     */
    virtual void
    setMaxEvaluations(
        int max_evaluations) = 0;

    /*!
     * \brief Get the maximum number of function evaluations to use per solve.
     */
    virtual int
    getMaxEvaluations() const = 0;

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

    /*!
     * \brief Set the tolerance in terms of the norm of the change in the
     * solution between steps.
     */
    virtual void
    setSolutionTolerance(
        double solution_tol) = 0;

    /*!
     * \brief Get the tolerance in terms of the norm of the change in the
     * solution between steps.
     */
    virtual double
    getSolutionTolerance() const = 0;

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent nonlinear solve.
     */
    virtual int
    getNumIterations() const = 0;

    /*!
     * \brief Return the number of linear iterations from the most recent
     * nonlinear solve.
     */
    virtual int
    getNumLinearIterations() const = 0;

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
    NewtonKrylovSolver(
        const NewtonKrylovSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    NewtonKrylovSolver&
    operator=(
        const NewtonKrylovSolver& that);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/NewtonKrylovSolver.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_NewtonKrylovSolver
