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

#ifndef included_IBTK_PETScPCLSWrapper
#define included_IBTK_PETScPCLSWrapper

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/LinearSolver.h"

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "tbox/Pointer.h"

#include "petscpc.h"
#include "petscvec.h"

#include <string>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PETScPCLSWrapper provides a LinearSolver interface for a <A
 * HREF="http://www.mcs.anl.gov/petsc">PETSc</A> PC object.
 */
class PETScPCLSWrapper : public LinearSolver
{
public:
    /*!
     * \brief Constructor.
     *
     * \param object_name name of the solver
     * \param petsc_pc PETSc PC object "wrapped" by this object
     *
     */
    PETScPCLSWrapper(std::string object_name, PC petsc_pc);

    /*!
     * \brief Destructor.
     *
     * \note The argument \a petsc_pc provided to the class constructor <em>is
     * not</em> deallocated by the class destructor.
     */
    ~PETScPCLSWrapper();

    /*!
     * \name Functions to access the underlying PETSc PC object.
     */
    //\{

    /*!
     * \brief Get the PETSc PC object "wrapped" by this object.
     */
    const PC& getPETScPC() const;

    //\}

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
     * -# initialize the hierarchy-dependent data required by the
     *    solver via initializeSolverState(),
     * -# solve the system one or more times via solveSystem(), and
     * -# remove the hierarchy-dependent data via
     *    deallocateSolverState().
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
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    void setInitialGuessNonzero(bool initial_guess_nonzero = true) override;

    /*!
     * \brief Get whether the initial guess is non-zero.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    bool getInitialGuessNonzero() const override;

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    void setMaxIterations(int max_iterations) override;

    /*!
     * \brief Get the maximum number of iterations to use per solve.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    int getMaxIterations() const override;

    /*!
     * \brief Set the absolute residual tolerance for stopping.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    void setAbsoluteTolerance(double abs_residual_tol) override;

    /*!
     * \brief Get the absolute residual tolerance for stopping.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    double getAbsoluteTolerance() const override;

    /*!
     * \brief Set the relative residual tolerance for stopping.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    void setRelativeTolerance(double rel_residual_tol) override;

    /*!
     * \brief Get the relative residual tolerance for stopping.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    double getRelativeTolerance() const override;

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent linear solve.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    int getNumIterations() const override;

    /*!
     * \brief Return the residual norm from the most recent iteration.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    double getResidualNorm() const override;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScPCLSWrapper() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScPCLSWrapper(const PETScPCLSWrapper& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScPCLSWrapper& operator=(const PETScPCLSWrapper& that) = delete;

    const PC d_petsc_pc;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_x, d_b;
    Vec d_petsc_x = nullptr, d_petsc_b = nullptr;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PETScPCLSWrapper
