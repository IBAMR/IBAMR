// Filename: KrylovLinearSolverWrapper.h
// Created on 17 Aug 2012 by Boyce Griffith
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

#ifndef included_KrylovLinearSolverWrapper
#define included_KrylovLinearSolverWrapper

/////////////////////////////// INCLUDES /////////////////////////////////////

// IBTK INCLUDES
#include <ibtk/KrylovLinearSolver.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class KrylovLinearSolverWrapper provides a wrapper for
 * KrylovLinearSolvers.  This class is not useful itself, but subclasses can be
 * used to extend the functionality of an existing Krylov solver class.
 *
 * \note To ensure consistency of all solver class data members, it is
 * recommended that the wrapped solver be accessed \em only though the wrapper
 * class.
 */
class KrylovLinearSolverWrapper
    : public KrylovLinearSolver
{
public:
    /*!
     * Constructor.
     */
    KrylovLinearSolverWrapper(
        SAMRAI::tbox::Pointer<KrylovLinearSolver> krylov_solver);

    /*!
     * Destructor.
     */
    ~KrylovLinearSolverWrapper();

    /*!
     * \brief Get the wrapped solver.
     */
    SAMRAI::tbox::Pointer<KrylovLinearSolver>
    getWrappedSolver() const;

    /*!
     * \name General-purpose solver functionality.
     */
    //\{

    /*!
     * \brief Return whether the operator is initialized.
     */
    bool
    getIsInitialized() const;

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
     * \brief Set the HierarchyMathOps object used by the solver.
     */
    void
    setHierarchyMathOps(
        SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops);

    /*!
     * \brief Get the HierarchyMathOps object used by the solver.
     */
    SAMRAI::tbox::Pointer<HierarchyMathOps>
    getHierarchyMathOps() const;

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
        bool contains_constant_vec,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >& nullspace_basis_vecs=std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >());

    /*!
     * \brief Get whether the nullspace of the linear system contains th
     * constant vector.
     */
    bool
    getNullspaceContainsConstantVector() const;

    /*!
     * \brief Get the basis vectors for the nullspace of the linear system.
     */
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > >&
    getNullspaceBasisVectors() const;

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
    setLoggingEnabled(
        bool enable_logging);

    /*!
     * \brief Determine whether logging is enabled or disabled.
     */
    bool
    getLoggingEnabled() const;

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
    KrylovLinearSolverWrapper();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    KrylovLinearSolverWrapper(
        const KrylovLinearSolverWrapper& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    KrylovLinearSolverWrapper&
    operator=(
        const KrylovLinearSolverWrapper& that);
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/KrylovLinearSolverWrapper.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_KrylovLinearSolverWrapper
