// Filename: PETScKrylovLinearSolver.h
// Created on 09 Sep 2003 by Boyce Griffith
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

#ifndef included_PETScKrylovLinearSolver
#define included_PETScKrylovLinearSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscksp.h>

// IBTK INCLUDES
#include <ibtk/KrylovLinearSolver.h>
#include <ibtk/LinearOperator.h>
#include <ibtk/LinearSolver.h>

// SAMRAI INCLUDES
#include <SAMRAIVectorReal.h>
#include <tbox/Pointer.h>

// C++ STDLIB INCLUDES
#include <ostream>
#include <string>
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScKrylovLinearSolver provides a KrylovLinearSolver interface
 * for a <A HREF="http://www-unix.mcs.anl.gov/petsc">PETSc</A> Krylov subspace
 * iterative linear solver (KSP).
 *
 * This solver class provides access to a large number of Krylov subspace
 * solvers for linear problems of the form \f$Ax=b\f$ using the PETSc KSP linear
 * solver interface.  See <A
 * HREF="http://www-unix.mcs.anl.gov/petsc/petsc-as/documentation/linearsolvertable.html">
 * http://www-unix.mcs.anl.gov/petsc/petsc-as/documentation/linearsolvertable.html</A>
 * for a complete list of the Krylov solvers provided by this class.  Note that
 * solver configuration is typically done at runtime via command line options.
 *
 * \note
 * - Preconditioners and direct solvers provided by PETSc \em cannot be used
 *   within the present IBTK solver framework.  However, this <em>does not
 *   mean</em> that preconditioners cannot be used with the
 *   PETScKrylovLinearSolver class; see, for instance, classes FACPreconditioner
 *   and CCPoissonFACOperator.  \par
 * - Users have the \em option of supplying the PETScKrylovLinearSolver class
 *   constructor with a KSP object.  It is important to emphasize that doing so
 *   is optional.  When a KSP object is provided to the class constructor, the
 *   PETScKrylovLinearSolver class acts as a \em wrapper for the supplied KSP
 *   object.  \par
 * - The functionality of the PETScKrylovLinearSolver class is essentially
 *   identical regardless as to whether a PETSc KSP object is provided to the
 *   class constructor.  The main exception is that when an external KSP object
 *   is provided to the class constructor, memory management of that object is
 *   \em NOT handled by the PETScKrylovLinearSolver.  In particular, it is the
 *   caller's responsibility to ensure that the supplied KSP object is properly
 *   destroyed via KSPDestroy().
 *
 * PETSc is developed in the Mathematics and Computer Science (MCS) Division at
 * Argonne National Laboratory (ANL).  For more information about PETSc, see <A
 * HREF="http://www-unix.mcs.anl.gov/petsc">
 * http://www-unix.mcs.anl.gov/petsc</A>.
 */
class PETScKrylovLinearSolver
    : public KrylovLinearSolver
{
public:
    /*!
     * \brief Constructor for a concrete KrylovLinearSolver that employs the
     * PETSc KSP solver framework.
     *
     * \param object_name     Name of the solver
     * \param options_prefix  Prefix for accessing options set through the PETSc options database (optional)
     * \param petsc_comm      MPI communicator
     *
     * \note The value of \a petsc_comm is used to specify the MPI communicator
     * used when initializing any PETSc objects required by this class.
     */
    PETScKrylovLinearSolver(
        const std::string& object_name,
        const std::string& options_prefix="",
        MPI_Comm petsc_comm=PETSC_COMM_WORLD);

    /*!
     * \brief Constructor for a concrete KrylovLinearSolver that acts as a
     * "wrapper" for a supplied PETSc KSP object.
     *
     * \param object_name     Name of the solver
     * \param petsc_ksp       PETSc KSP object
     * \param options_prefix  Prefix for accessing options set through the PETSc options database (optional)
     *
     * \note This constructor initializes a PETScKrylovLinearSolver object that
     * acts as a "wrapper" for the provided KSP object.  Note that memory
     * management of the provided KSP object is \em NOT handled by this class.
     */
    PETScKrylovLinearSolver(
        const std::string& object_name,
        const KSP& petsc_ksp,
        const std::string& options_prefix="");

    /*!
     * \brief Destructor.
     */
    ~PETScKrylovLinearSolver();

    /*!
     * \brief Set the KSP type.
     */
    void
    setKSPType(
        const std::string& ksp_type);

    /*!
     * \brief Set a list of valid preconditioner types, set via the
     * -pc_shell_type options database flag.
     */
    void
    setValidPCShellTypes(
        const std::vector<std::string>& pc_shell_types);

    /*!
     * \brief Get the preconditioner type, set via the -pc_shell_type options
     * database flag.
     *
     * \note If the PC type is not PCSHELL, or if no preconditioner types have
     * been registered with the solver via setValidPCShellTypes, this method
     * will return the value "none".
     */
    const std::string&
    getPCShellType() const;

    /*!
     * \name Functions to access the underlying PETSc objects.
     */
    //\{

    /*!
     * \brief Get the PETSc KSP object.
     */
    const KSP&
    getPETScKSP() const;

    //\}

    /*!
     * \name Krylov solver functionality.
     */
    //\{

    /*!
     * \brief Set the current time interval (for a time-dependent solver).
     */
    void
    setTimeInterval(
        const double current_time,
        const double new_time);

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

    /*!
     * \brief Set the nullspace of the linear system.
     *
     * The basis vector, if any, will be normalized by the solver.
     */
    void
    setNullspace(
        const bool contains_constant_vector,
        SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > nullspace_basis_vec);

    /*!
     * \brief Set the nullspace of the linear system.
     *
     * Basis vectors must be orthogonal but are not required to be orthonormal.
     * Basis vectors will be normalized automatically.
     */
    void
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
     * When linear operator or preconditioner objects have been registered with
     * this class via setOperator() and setPreconditioner(), they are also
     * initialized by this member function.
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
     * When linear operator or preconditioner objects have been registered with
     * this class via setOperator() and setPreconditioner(), they are also
     * deallocated by this member function.
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

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScKrylovLinearSolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScKrylovLinearSolver(
        const PETScKrylovLinearSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScKrylovLinearSolver&
    operator=(
        const PETScKrylovLinearSolver& that);

    /*!
     * \brief Common routine used by all class constructors.
     */
    void
    common_ctor();

    /*!
     * \brief Report the KSPConvergedReason.
     */
    void
    reportKSPConvergedReason(
        const KSPConvergedReason& reason,
        std::ostream& os) const;

    /*!
     * \brief Reset the values of the convergence tolerances for the PETSc KSP
     * object.
     */
    void
    resetKSPOptions();

    /*!
     * \brief Reset the KSP operators to correspond to the supplied
     * LinearOperator.
     */
    void
    resetKSPOperators();

    /*!
     * \brief Reset the KSP PC to correspond to the supplied preconditioner.
     */
    void
    resetKSPPC();

    /*!
     * \brief Reset the KSP nullspace object to correspond to the supplied
     * nullspace basis vectors.
     */
    void
    resetKSPNullspace();

    /*!
     * \brief Destroy data allocated to describe nullspace.
     */
    void
    deallocateNullspaceData();

    /*!
     * \name Static functions for use by PETSc KSP and MatShell objects.
     */
    //\{

    /*!
     * \brief Set PC options from the PETSc options database.
     */
    static PetscErrorCode
    PCShellSetFromOptions_SAMRAI(
        PC pc);

    /*!
     * \brief Compute the matrix vector product \f$y=Ax\f$.
     */
    static PetscErrorCode
    MatVecMult_SAMRAI(
        Mat A,
        Vec x,
        Vec y);

    /*!
     * \brief Compute the matrix vector product \f$z=Ax+y\f$.
     */
    static PetscErrorCode
    MatVecMultAdd_SAMRAI(
        Mat A,
        Vec x,
        Vec y,
        Vec z);

    /*!
     * \brief Compute the matrix-transpose vector product \f$y=A^{T}x\f$.
     */
    static PetscErrorCode
    MatVecMultTranspose_SAMRAI(
        Mat A,
        Vec x,
        Vec y);

    /*!
     * \brief Compute the matrix-transpose vector product \f$y=A^{T}x+z\f$.
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
     * \brief Apply the preconditioner to \a x and store the result in \a y.
     */
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 0)
    static PetscErrorCode
    PCApply_SAMRAI(
        void* ctx,
        Vec x,
        Vec y);
#endif
#if (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1)
    static PetscErrorCode
    PCApply_SAMRAI(
        PC pc,
        Vec x,
        Vec y);
#endif

    //\}

    std::string d_object_name;

    std::string d_ksp_type;

    std::vector<std::string> d_pc_shell_types;
    std::string d_pc_shell_type;

    bool d_is_initialized, d_reinitializing_solver;
    bool d_do_log;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_solver_x, d_solver_b;
    Vec d_petsc_x, d_petsc_b;

    std::string d_options_prefix;

    MPI_Comm     d_petsc_comm;
    KSP          d_petsc_ksp;
    Mat          d_petsc_mat;
    MatNullSpace d_petsc_nullsp;
    bool d_managing_petsc_ksp;
    bool d_user_provided_mat;
    bool d_user_provided_pc;

    SAMRAI::tbox::Pointer<LinearOperator> d_A;
    SAMRAI::tbox::Pointer<LinearSolver> d_pc_solver;

    bool d_nullsp_contains_constant_vector;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_solver_nullsp_constant;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > > d_solver_nullsp_basis;
    Vec d_petsc_nullsp_constant;
    std::vector<Vec> d_petsc_nullsp_basis;
    bool d_solver_has_attached_nullsp;

    bool d_initial_guess_nonzero;

    double d_rel_residual_tol;
    double d_abs_residual_tol;
    double d_divergence_tol;
    int d_max_iterations;

    int d_current_its;
    double d_current_residual_norm;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibtk/PETScKrylovLinearSolver.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PETScKrylovLinearSolver
