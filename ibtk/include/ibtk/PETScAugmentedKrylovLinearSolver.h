// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDE GUARD ////////////////////////////////

#ifndef included_IBTK_PETScAugmentedKrylovLinearSolver
#define included_IBTK_PETScAugmentedKrylovLinearSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/PETScLinearAugmentedOperator.h"

#include "IntVector.h"
#include "MultiblockDataTranslator.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscsys.h"
#include "petscvec.h"

#include <mpi.h>

#include <iosfwd>
#include <string>
#include <vector>

namespace IBTK
{
class LinearOperator;
} // namespace IBTK

/////////////////////////////// FORWARD DECLARATION //////////////////////////

namespace IBTK
{
}

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScAugmentedKrylovLinearSolver provides a KrylovLinearSolver interface
 * for a <A HREF="http://www.mcs.anl.gov/petsc">PETSc</A> Krylov subspace
 * iterative linear solver (KSP).
 *
 * This solver class provides access to a large number of Krylov subspace
 * solvers for linear problems of the form \f$[A; A_{aug}][x; x_{aug}]=[b; b_{aug}]\f$ using the PETSc KSP linear
 * solver interface.  See <A
 *
 HREF="http://www.mcs.anl.gov/petsc/documentation/linearsolvertable.html">http://www.mcs.anl.gov/petsc/documentation/linearsolvertable.html</A>
 * for a complete list of the Krylov solvers provided by this class.  Note that
 * solver configuration is typically done at runtime via command line options.
 *
 * \note
 * - Currently, preconditioners can not be used with this class.
 * - Users must register the following with this class
 *    -- a PETScLinearAugmentedOperator via the setOperator class.
 *    -- a Vec object that contains the data storage for the augmented vector
 *    -- a Mat object computes the action of the augmented operator (A_{aug} * [x; x_{aug}]).
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 options_prefix = ""           // see setOptionsPrefix()
 ksp_type = "gmres"            // see setKSPType()
 initial_guess_nonzero = TRUE  // see setInitialGuessNonzero()
 rel_residual_tol = 1.0e-5     // see setRelativeTolerance()
 abs_residual_tol = 1.0e-50    // see setAbsoluteTolerance()
 max_iterations = 10000        // see setMaxIterations()
 enable_logging = FALSE        // see setLoggingEnabled()
 \endverbatim
 *
 * PETSc is developed in the Mathematics and Computer Science (MCS) Division at
 * Argonne National Laboratory (ANL).  For more information about PETSc, see <A
 * HREF="http://www.mcs.anl.gov/petsc">http://www.mcs.anl.gov/petsc</A>.
 */
class PETScAugmentedKrylovLinearSolver : public KrylovLinearSolver
{
public:
    /*!
     * \brief Constructor for a concrete KrylovLinearSolver that employs the
     * PETSc KSP solver framework.
     */
    PETScAugmentedKrylovLinearSolver(std::string object_name,
                                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                     std::string default_options_prefix,
                                     MPI_Comm petsc_comm = PETSC_COMM_WORLD);

    /*!
     * \brief Destructor.
     */
    ~PETScAugmentedKrylovLinearSolver();

    /*!
     * \brief Set the KSP type.
     */
    void setKSPType(const std::string& ksp_type);

    /*!
     * \brief Set the options prefix used by this PETSc solver object.
     */
    void setOptionsPrefix(const std::string& options_prefix);

    /*!
     * \name Functions to set up augmented components of the system.
     */
    //\{
    /*!
     * \brief Set the augmented matrix. Note that the matrix should be set up to take a nested Vec containing Eulerian
     * dofs and augmented dofs, and return the action on the augmented DOFs.
     */
    void setAugmentedRHS(const Vec& vec);

    /*!
     * \brief Set the augmented vector. Note that the Vec sets up the internal data structures.
     */
    void setAugmentedVec(const Vec& vec);
    //\}

    /*!
     * \name Functions to access the underlying PETSc objects.
     */
    //\{

    /*!
     * \brief Get the PETSc KSP object.
     */
    const KSP& getPETScKSP() const;

    const Mat& getAugmentedMat() const;

    const Vec& getAugmentedVec() const;
    //\}

    /*!
     * \name Krylov solver functionality.
     */
    //\{

    /*!
     * \brief Set the linear operator used when solving \f$Ax=b\f$.
     */
    void setOperator(SAMRAI::tbox::Pointer<LinearOperator> A) override;

    /*!
     * \brief Set the preconditioner used by the Krylov subspace method when
     * solving \f$Ax=b\f$.
     *
     * \note If the preconditioner is NULL, no preconditioning is performed.
     */
    void setPreconditioner(SAMRAI::tbox::Pointer<LinearSolver> pc_solver = NULL) override;

    /*!
     * \brief Set the nullspace of the linear system.
     *
     * Basis vectors must be orthogonal but are not required to be orthonormal.
     * Basis vectors will be normalized automatically.
     */
    void setNullspace(
        bool contains_constant_vec,
        const std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > >& nullspace_basis_vecs =
            std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > >()) override;

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
    bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                     SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

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
    void initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                               const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) override;

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
    void deallocateSolverState() override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScAugmentedKrylovLinearSolver() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScAugmentedKrylovLinearSolver(const PETScAugmentedKrylovLinearSolver& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScAugmentedKrylovLinearSolver& operator=(const PETScAugmentedKrylovLinearSolver& that) = delete;

    /*!
     * \brief Common routine used by all class constructors.
     */
    void common_ctor();

    /*!
     * \brief Reset the KSP wrapped by this solver class.
     */
    void resetWrappedKSP(KSP& petsc_ksp);

    /*!
     * \brief Reset the values of the convergence tolerances for the PETSc KSP
     * object.
     */
    void resetKSPOptions();

    /*!
     * \brief Reset the KSP operators to correspond to the supplied
     * LinearOperator.
     */
    void resetKSPOperators();

    /*!
     * \brief Reset the KSP PC to correspond to the supplied preconditioner.
     */
    void resetKSPPC();

    /*!
     * \brief Reset the Mat nullspace object to correspond to the supplied
     * nullspace basis vectors.
     */
    void resetMatNullspace();

    /*!
     * \brief Destroy data allocated to describe nullspace.
     */
    void deallocateNullspaceData();

    /*!
     * \name Static functions for use by PETSc KSP and MatShell objects.
     */
    //\{

    /*!
     * \brief Compute the matrix vector product \f$y=Ax\f$.
     *
     * \note The x and y Vec's are block matrices containing the Eulerian and augmented DOFs.
     */
    static PetscErrorCode MatVecMult_SAMRAI(Mat A, Vec x, Vec y);
    //\}

    std::string d_ksp_type;

    bool d_reinitializing_solver = false;

    // Data structures for entire system, Eulerian + augmented dofs.
    // Note these are Vec and Mat nests.
    KSP d_petsc_ksp = nullptr;
    Mat d_petsc_mat = nullptr;
    Vec d_petsc_x = nullptr, d_petsc_b = nullptr;

    // Data structures for Eulerian dofs.
    Vec d_eul_x = nullptr, d_eul_b = nullptr;

    // Data structures for augmented DOFs
    // Note the Mat takes in a nested vec and returns an Augmented vector.
    Vec d_aug_vec, d_aug_x, d_aug_b;

    std::string d_options_prefix;

    MPI_Comm d_petsc_comm;
    MatNullSpace d_petsc_nullsp = nullptr;
    // bool d_managing_petsc_ksp = true;
    // bool d_user_provided_mat = false;
    bool d_user_provided_pc = false;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_nullspace_constant_vec;
    Vec d_petsc_nullspace_constant_vec = nullptr;
    std::vector<Vec> d_petsc_nullspace_basis_vecs;
    bool d_solver_has_attached_nullspace = false;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PETScAugmentedKrylovLinearSolver
