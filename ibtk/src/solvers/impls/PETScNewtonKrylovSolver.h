// Filename: PETScNewtonKrylovSolver.h
// Created on 26 Nov 2003 by Boyce Griffith
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

#ifndef included_PETScNewtonKrylovSolver
#define included_PETScNewtonKrylovSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscsnes.h>

// IBTK INCLUDES
#include <ibtk/GeneralOperator.h>
#include <ibtk/JacobianOperator.h>
#include <ibtk/NewtonKrylovSolver.h>
#include <ibtk/PETScKrylovLinearSolver.h>

// SAMRAI INCLUDES
#include <SAMRAIVectorReal.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <ostream>
#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScNewtonKrylovSolver provides a NewtonKrylovSolver interface
 * for a <A HREF="http://www-unix.mcs.anl.gov/petsc">PETSc</A> inexact
 * Newton-Krylov iterative nonlinear solver (SNES).
 *
 * This solver class provides access to inexact Newton-Krylov methods, including
 * line search and trust region Newton methods, using the PETSc SNES nonlinear
 * solver interface.  Note that solver configuration is typically done at
 * runtime via command line options.
 *
 * \note
 * - Users have the \em option of supplying the PETScNewtonKrylovSolver class
 *   constructor with a SNES object.  It is important to emphasize that doing so
 *   is optional.  When a SNES object is provided to the class constructor, the
 *   PETScNewtonKrylovSolver class acts as a \em wrapper for the supplied SNES
 *   object.  \par
 * - The functionality of the PETScNewtonKrylovSolver class is essentially
 *   identical regardless as to whether a PETSc SNES object is provided to the
 *   class constructor.  The main exception is that when an external SNES object
 *   is provided to the class constructor, memory management of that object is
 *   \em NOT handled by the PETScNewtonKrylovSolver.  In particular, it is the
 *   caller's responsibility to ensure that the supplied SNES object is properly
 *   destroyed via SNESDestroy().
 *
 * PETSc is developed in the Mathematics and Computer Science (MCS) Division at
 * Argonne National Laboratory (ANL).  For more information about PETSc, see <A
 * HREF="http://www-unix.mcs.anl.gov/petsc">
 * http://www-unix.mcs.anl.gov/petsc</A>.
 */
class PETScNewtonKrylovSolver
    : public NewtonKrylovSolver
{
public:
    /*!
     * \brief Constructor for a concrete NewtonKrylovSolver that employs the
     * PETSc SNES solver framework.
     *
     * \param object_name     Name of the solver
     * \param options_prefix  Prefix for accessing options set through the PETSc options database (optional)
     * \param petsc_comm      MPI communicator
     *
     * \note The value of \a petsc_comm is used to specify the MPI communicator
     * used when initializing any PETSc objects required by this class.
     */
    PETScNewtonKrylovSolver(
        const std::string& object_name,
        const std::string& options_prefix="",
        MPI_Comm petsc_comm=PETSC_COMM_WORLD);

    /*!
     * \brief Constructor for a concrete NewtonKrylovSolver that acts as a
     * "wrapper" for a supplied PETSc SNES object.
     *
     * \param object_name     Name of the solver
     * \param petsc_snes      PETSc SNES object
     * \param options_prefix  Prefix for accessing options set through the PETSc options database (optional)
     *
     * \note This constructor initializes a PETScNewtonKrylovSolver object that
     * acts as a "wrapper" for the provided SNES object.  Note that memory
     * management of the provided SNES object is \em NOT handled by this class.
     */
    PETScNewtonKrylovSolver(
        const std::string& object_name,
        const SNES& petsc_snes,
        const std::string& options_prefix="");

    /*!
     * \brief Destructor.
     */
    ~PETScNewtonKrylovSolver();

    /*!
     * \name Functions to access the underlying PETSc objects.
     */
    //\{

    /*!
     * \brief Get the PETSc SNES object.
     */
    const SNES&
    getPETScSNES() const;

    //\}

    /*!
     * \name Newton-Krylov solver functionality.
     */
    //\{

    /*!
     * \brief Set the current time interval (for a time-dependent solver).
     */
    void
    setTimeInterval(
        double current_time,
        double new_time);

    /*!
     * \brief Set the nonlinear operator \f$F[x]\f$ used by the solver.
     */
    void
    setOperator(
        SAMRAI::tbox::Pointer<GeneralOperator> op);

    /*!
     * \brief Retrieve the nonlinear operator \f$F[x]\f$ used by the solver.
     */
    SAMRAI::tbox::Pointer<GeneralOperator>
    getOperator() const;

    /*!
     * \brief Return the vector in which the approximate solution is stored.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >
    getSolutionVector() const;

    /*!
     * \brief Return the vector in which the nonlinear function evaluation is
     * stored.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> >
    getFunctionVector() const;

    /*!
     * \brief Set the Jacobian operator \f$J[x] = F'[x]\f$ used by the solver.
     *
     * \note If a Jacobian object is not explicitly provided to the solver, a
     * Jacobian-free inexact Newton-Krylov method is employed to approximate the
     * action of the Jacobian.
     */
    void
    setJacobian(
        SAMRAI::tbox::Pointer<JacobianOperator> J);

    /*!
     * \brief Retrieve the Jacobian operator \f$J[x] = F'[x]\f$ used by the
     * solver.
     */
    SAMRAI::tbox::Pointer<JacobianOperator>
    getJacobian() const;

    /*!
     * \brief Retrieve the Krylov linear solver used in computing Newton step
     * directions.
     */
    SAMRAI::tbox::Pointer<KrylovLinearSolver>
    getLinearSolver() const;

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
    bool
    solveSystem(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Compute hierarchy dependent data required for solving
     * \f$F[x]=b\f$.
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
     * When any operator objects have been registered with this class via
     * setOperator() or setJacobian(), they are also initialized by this member
     * function.
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
    void
    initializeSolverState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * When any operator objects have been registered with this class via
     * setOperator() or setJacobian(), they are also deallocated by this member
     * function.
     *
     * \note It is safe to call deallocateSolverState() when the solver state is
     * already deallocated.
     *
     * \see initializeSolverState
     */
    void
    deallocateSolverState();

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set the maximum number of nonlinear iterations to use per solve.
     */
    void
    setMaxIterations(
        int max_iterations);

    /*!
     * \brief Get the maximum number of nonlinear iterations to use per solve.
     */
    int
    getMaxIterations() const;

    /*!
     * \brief Set the maximum number of function evaluations to use per solve.
     */
    void
    setMaxEvaluations(
        int max_evaluations);

    /*!
     * \brief Get the maximum number of function evaluations to use per solve.
     */
    int
    getMaxEvaluations() const;

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

    /*!
     * \brief Set the tolerance in terms of the norm of the change in the
     * solution between steps.
     */
    void
    setSolutionTolerance(
        double solution_tol);

    /*!
     * \brief Get the tolerance in terms of the norm of the change in the
     * solution between steps.
     */
    double
    getSolutionTolerance() const;

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent nonlinear solve.
     */
    int
    getNumIterations() const;

    /*!
     * \brief Return the number of linear iterations from the most recent
     * nonlinear solve.
     */
    int
    getNumLinearIterations() const;

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

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScNewtonKrylovSolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScNewtonKrylovSolver(
        const PETScNewtonKrylovSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScNewtonKrylovSolver&
    operator=(
        const PETScNewtonKrylovSolver& that);

    /*!
     * \brief Common routine used by all class constructors.
     */
    void
    common_ctor();

    /*!
     * \brief Report the SNESConvergedReason.
     */
    void
    reportSNESConvergedReason(
        const SNESConvergedReason& reason,
        std::ostream& os) const;

    /*!
     * \brief Reset the values of the convergence tolerances for the PETSc SNES
     * object.
     */
    void
    resetSNESOptions();

    /*!
     * \brief Reset the function for the PETSc SNES object.
     */
    void
    resetSNESFunction();

    /*!
     * \brief Reset the Jacobian for the PETSc SNES object.
     */
    void
    resetSNESJacobian();

    /*!
     * \name Static functions for use by PETSc SNES and MatShell objects.
     */
    //\{

    /*!
     * \brief Evaluate f = F[x].
     */
    static PetscErrorCode
    FormFunction_SAMRAI(
        SNES snes,
        Vec x,
        Vec f,
        void* p_ctx);

    /*!
     * \brief Setup F'[x].
     */
    static PetscErrorCode
    FormJacobian_SAMRAI(
        SNES snes,
        Vec x,
        Mat* A,
        Mat* B,
        MatStructure* mat_structure,
        void* p_ctx);

    /*!
     * \brief Compute the matrix vector product y = Ax.
     */
    static PetscErrorCode
    MatVecMult_SAMRAI(
        Mat A,
        Vec x,
        Vec y);

    /*!
     * \brief Compute the matrix vector product y = Ax + z.
     */
    static PetscErrorCode
    MatVecMultAdd_SAMRAI(
        Mat A,
        Vec x,
        Vec y,
        Vec z);

    /*!
     * \brief Compute the matrix-transpose vector product y = A'x.
     */
    static PetscErrorCode
    MatVecMultTranspose_SAMRAI(
        Mat A,
        Vec x,
        Vec y);

    /*!
     * \brief Compute the matrix-transpose vector product y = A'x + z.
     */
    static PetscErrorCode
    MatVecMultTransposeAdd_SAMRAI(
        Mat A,
        Vec x,
        Vec y,
        Vec z);

    /*!
     * \brief Get vector(s) compatible with the matrix, i.e., with the same
     * parallel layout.
     */
    static PetscErrorCode
    MatGetVecs_SAMRAI(
        Mat mat,
        Vec* right,
        Vec* left);

    /*!
     * \brief Apply the preconditioner to x.
     */
    static PetscErrorCode
    PCApply_SAMRAI(
        void* ctx,
        Vec x,
        Vec y);

    //\}

    std::string d_object_name;

    bool d_is_initialized, d_reinitializing_solver;
    bool d_do_log;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_solver_x, d_solver_b, d_solver_r;
    Vec d_petsc_x, d_petsc_b, d_petsc_r;

    std::string d_options_prefix;

    MPI_Comm d_petsc_comm;
    SNES     d_petsc_snes;
    Mat      d_petsc_jac;
    bool d_managing_petsc_snes;
    bool d_user_provided_function;
    bool d_user_provided_jacobian;

    SAMRAI::tbox::Pointer<GeneralOperator>         d_F;
    SAMRAI::tbox::Pointer<JacobianOperator>        d_J;
    SAMRAI::tbox::Pointer<PETScKrylovLinearSolver> d_krylov_solver;

    double d_abs_residual_tol;
    double d_rel_residual_tol;
    double d_solution_tol;
    int d_max_iterations;
    int d_max_evaluations;

    int d_current_its, d_current_lits;
    double d_current_residual_norm;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/PETScNewtonKrylovSolver.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PETScNewtonKrylovSolver
