// Filename: PoissonKrylovLinearSolverWrapper.h
// Created on 13 Aug 2012 by Boyce Griffith
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

#ifndef included_PoissonKrylovLinearSolverWrapper
#define included_PoissonKrylovLinearSolverWrapper

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/KrylovLinearSolver.h>
#include <ibtk/PoissonSolver.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PoissonKrylovLinearSolverWrapper provides a wrapper for
 * KrylovLinearSolvers that are to be used as Poisson solvers.
 */
class PoissonKrylovLinearSolverWrapper
    : public KrylovLinearSolver,
      public PoissonSolver
{
public:
    /*!
     * Constructor.
     */
    PoissonKrylovLinearSolverWrapper(
        SAMRAI::tbox::Pointer<KrylovLinearSolver> krylov_solver);

    /*!
     * Destructor.
     */
    ~PoissonKrylovLinearSolverWrapper();

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
     * \note \a bc_coef may be NULL.  In this case, default boundary conditions
     * (as supplied to the class constructor) are employed.
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
     * default boundary conditions (as supplied to the class constructor) are
     * employed for that data depth.
     *
     * \param bc_coefs  Vector of pointers to objects that can set the Robin boundary condition coefficients
     */
    void
    setPhysicalBcCoefs(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * \name General-purpose solver functionality.
     */
    //\{

    /*!
     * \brief Set whether the solver should use homogeneous boundary conditions.
     */
    void
    setHomogeneousBc(
        bool homogeneous_bc);

    /*!
     * \brief Return whether the solver is using homogeneous boundary
     * conditions.
     */
    bool
    getHomogeneousBc() const;

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    void
    setSolutionTime(
        double solution_time);

    /*!
     * \brief Get the time at which the solution is being evaluated.
     */
    double
    getSolutionTime() const;

    /*!
     * \brief Set the current time interval.
     */
    void
    setTimeInterval(
        double current_time,
        double new_time);

    /*!
     * \brief Get the current time interval.
     */
    std::pair<double,double>
    getTimeInterval() const;

    /*!
     * \brief Solve the system of equations.
     */
    bool
    solveSystem(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Compute hierarchy dependent data required for solving
     * \f$F[x]=b\f$.
     */
    void
    initializeSolverState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     */
    void
    deallocateSolverState();

    //\}

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Set the nullspace of the linear system.
     */
    void
    setNullspace(
        bool contains_constant_vector,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >& nullspace_basis_vecs=std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >());

    //\}

    /*!
     * \name Krylov solver functionality.
     */
    //\{

    /*!
     * \brief Set the linear operator used when solving \f$Ax=b\f$.
     */
    void
    setOperator(
        SAMRAI::tbox::Pointer<LinearOperator> A);

    /*!
     * \brief Retrieve the linear operator used when solving \f$Ax=b\f$.
     */
    SAMRAI::tbox::Pointer<LinearOperator>
    getOperator() const;

    /*!
     * \brief Set the preconditioner used by the Krylov subspace method when
     * solving \f$Ax=b\f$.
     *
     * \note If the preconditioner is NULL, no preconditioning is performed.
     */
    void
    setPreconditioner(
        SAMRAI::tbox::Pointer<LinearSolver> pc_solver=NULL);

    /*!
     * \brief Retrieve the preconditioner used by the Krylov subspace method
     * when solving \f$Ax=b\f$.
     */
    SAMRAI::tbox::Pointer<LinearSolver>
    getPreconditioner() const;

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    void
    setInitialGuessNonzero(
        bool initial_guess_nonzero=true);

    /*!
     * \brief Get whether the initial guess is non-zero.
     */
    bool
    getInitialGuessNonzero() const;

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    void
    setMaxIterations(
        int max_iterations);

    /*!
     * \brief Get the maximum number of iterations to use per solve.
     */
    int
    getMaxIterations() const;

    /*!
     * \brief Set the absolute residual tolerance for convergence.
     */
    void
    setAbsoluteTolerance(
        double abs_residual_tol);

    /*!
     * \brief Get the absolute residual tolerance for convergence.
     */
    double
    getAbsoluteTolerance() const;

    /*!
     * \brief Set the relative residual tolerance for convergence.
     */
    void
    setRelativeTolerance(
        double rel_residual_tol);

    /*!
     * \brief Get the relative residual tolerance for convergence.
     */
    double
    getRelativeTolerance() const;

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent linear solve.
     */
    int
    getNumIterations() const;

    /*!
     * \brief Return the residual norm from the most recent iteration.
     */
    double
    getResidualNorm() const;

    //\}

    /*!
     * \name Logging functions.
     */
    //\{

    /*!
     * \brief Enable or disable logging.
     */
    void
    enableLogging(
        bool enabled=true);

    //\}

protected:
    /*!
     * The wrapped KrylovLinearSolver.
     */
    SAMRAI::tbox::Pointer<KrylovLinearSolver> d_krylov_solver;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PoissonKrylovLinearSolverWrapper();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PoissonKrylovLinearSolverWrapper(
        const PoissonKrylovLinearSolverWrapper& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PoissonKrylovLinearSolverWrapper&
    operator=(
        const PoissonKrylovLinearSolverWrapper& that);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/PoissonKrylovLinearSolverWrapper.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PoissonKrylovLinearSolverWrapper
