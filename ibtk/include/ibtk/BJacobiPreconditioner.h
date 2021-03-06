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

#ifndef included_IBTK_BJacobiPreconditioner
#define included_IBTK_BJacobiPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LinearSolver.h"

#include "tbox/Pointer.h"

#include <map>
#include <string>

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class BJacobiPreconditioner is a block Jacobi preconditioner which
 * implements the abstract LinearSolver interface.
 *
 * This solver class performs a single block Jacobi sweep, applying specified
 * component LinearSolver objects to the components of a supplied
 * SAMRAI::solv::SAMRAIVectorReal vector, and doing so in an additive fashion.
 * Note that the block Jacobi algorithm is not generally convergent, but can be
 * used as a preconditioner for a KrylovLinearSolver.
 *
 * \note Class BGaussSeidelPreconditioner implements the multiplicative (i.e.,
 * block Gauss-Seidel) version of this algorithm.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 initial_guess_nonzero = FALSE      // see setInitialGuessNonzero()
 rel_residual_tol = 1.0e-6          // see setRelativeTolerance()
 abs_residual_tol = 1.0e-30         // see setAbsoluteTolerance()
 max_iterations = 1                 // see setMaxIterations()
 \endverbatim
 *
 */
class BJacobiPreconditioner : public LinearSolver
{
public:
    /*!
     * \brief Constructor.
     */
    BJacobiPreconditioner(std::string object_name,
                          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                          const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~BJacobiPreconditioner();

    /*!
     * \brief Set the preconditioner to be employed on the specified vector
     * component.
     */
    void setComponentPreconditioner(SAMRAI::tbox::Pointer<LinearSolver> preconditioner, unsigned int component);

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
     * \note It is safe to call deallocateSolverState() when the solver state is
     * already deallocated.
     *
     * \see initializeSolverState
     */
    void deallocateSolverState() override;

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    void setInitialGuessNonzero(bool initial_guess_nonzero = true) override;

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    void setMaxIterations(int max_iterations) override;

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent linear solve.
     */
    int getNumIterations() const override;

    /*!
     * \brief Return the residual norm from the most recent iteration.
     */
    double getResidualNorm() const override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    BJacobiPreconditioner() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    BJacobiPreconditioner(const BJacobiPreconditioner& from) = delete;

    /*!
     * ]brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    BJacobiPreconditioner& operator=(const BJacobiPreconditioner& that) = delete;

    /*!
     * The component preconditioners.
     */
    std::map<unsigned int, SAMRAI::tbox::Pointer<LinearSolver> > d_pc_map;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_BJacobiPreconditioner
