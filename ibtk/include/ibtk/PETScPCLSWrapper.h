// Filename: PETScPCLSWrapper.h
// Created on 19 Oct 2003 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

#ifndef included_PETScPCLSWrapper
#define included_PETScPCLSWrapper

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>

#include "IntVector.h"
#include "SAMRAIVectorReal.h"
#include "ibtk/LinearSolver.h"
#include "petscpc.h"
#include "petscvec.h"
#include "tbox/Pointer.h"

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
    PETScPCLSWrapper(const std::string& object_name, const PC& petsc_pc);

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
    bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x, SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b);

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
                               const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * \note It is safe to call deallocateSolverState() when the solver state is
     * already deallocated.
     *
     * \see initializeSolverState
     */
    void deallocateSolverState();

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
    void setInitialGuessNonzero(bool initial_guess_nonzero = true);

    /*!
     * \brief Get whether the initial guess is non-zero.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    bool getInitialGuessNonzero() const;

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    void setMaxIterations(int max_iterations);

    /*!
     * \brief Get the maximum number of iterations to use per solve.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    int getMaxIterations() const;

    /*!
     * \brief Set the absolute residual tolerance for stopping.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    void setAbsoluteTolerance(double abs_residual_tol);

    /*!
     * \brief Get the absolute residual tolerance for stopping.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    double getAbsoluteTolerance() const;

    /*!
     * \brief Set the relative residual tolerance for stopping.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    void setRelativeTolerance(double rel_residual_tol);

    /*!
     * \brief Get the relative residual tolerance for stopping.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    double getRelativeTolerance() const;

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
    int getNumIterations() const;

    /*!
     * \brief Return the residual norm from the most recent iteration.
     *
     * \warning Not supported.  If called, a warning will be printed via
     * SAMRAI::TBOX_WARNING.
     */
    double getResidualNorm() const;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PETScPCLSWrapper();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PETScPCLSWrapper(const PETScPCLSWrapper& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PETScPCLSWrapper& operator=(const PETScPCLSWrapper& that);

    const PC d_petsc_pc;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_x, d_b;
    Vec d_petsc_x, d_petsc_b;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PETScPCLSWrapper
