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

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class NewtonKrylovSolver provides an abstract interface for the
 * implementation of inexact Newton-Krylov solvers for nonlinear problems of the
 * form \f$ F[x]=b \f$.
 */
class NewtonKrylovSolver
    : public GeneralSolver
{
public:
    /*!
     * \brief Constructor.
     */
    NewtonKrylovSolver(
        const std::string& object_name,
        bool homogeneous_bc=false);

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
        int max_iterations);

    /*!
     * \brief Get the maximum number of nonlinear iterations to use per solve.
     */
    virtual int
    getMaxIterations() const;

    /*!
     * \brief Set the maximum number of function evaluations to use per solve.
     */
    virtual void
    setMaxEvaluations(
        int max_evaluations);

    /*!
     * \brief Get the maximum number of function evaluations to use per solve.
     */
    virtual int
    getMaxEvaluations() const;

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

    /*!
     * \brief Set the tolerance in terms of the norm of the change in the
     * solution between steps.
     */
    virtual void
    setSolutionTolerance(
        double solution_tol);

    /*!
     * \brief Get the tolerance in terms of the norm of the change in the
     * solution between steps.
     */
    virtual double
    getSolutionTolerance() const;

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent nonlinear solve.
     */
    virtual int
    getNumIterations() const;

    /*!
     * \brief Return the number of linear iterations from the most recent
     * nonlinear solve.
     */
    virtual int
    getNumLinearIterations() const;

    /*!
     * \brief Return the residual norm from the most recent iteration.
     */
    virtual double
    getResidualNorm() const;

    //\}

protected:
    // Solver parameters.
    double d_rel_residual_tol;
    double d_abs_residual_tol;
    double d_solution_tol;
    int d_max_iterations;
    int d_max_evaluations;
    int d_current_its, d_current_lits;
    double d_current_residual_norm;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    NewtonKrylovSolver();

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
