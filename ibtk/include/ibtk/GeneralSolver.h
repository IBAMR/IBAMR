// Filename: GeneralSolver.h
// Created on 07 Apr 2012 by Boyce Griffith
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

#ifndef included_GeneralSolver
#define included_GeneralSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <iosfwd>
#include <string>
#include <utility>

#include "ibtk/HierarchyMathOps.h"
#include "tbox/DescribedClass.h"
#include "tbox/Pointer.h"

namespace SAMRAI
{
namespace solv
{
template <int DIM, class TYPE>
class SAMRAIVectorReal;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class GeneralSolver provides an abstract interface for the
 * implementation of linear or nonlinear solvers for systems of equations
 * defined on an AMR patch hierarchy.
 */
class GeneralSolver : public virtual SAMRAI::tbox::DescribedClass
{
public:
    /*!
     * \brief Constructor.
     */
    GeneralSolver();

    /*!
     * \brief Empty virtual destructor.
     */
    virtual ~GeneralSolver();

    /*!
     * \name General-purpose solver functionality.
     */
    //\{

    /*!
     * \brief Return the object name.
     */
    const std::string& getName() const;

    /*!
     * \brief Return whether the operator is initialized.
     */
    virtual bool getIsInitialized() const;

    /*!
     * \brief Set whether the solver should use homogeneous boundary conditions.
     */
    virtual void setHomogeneousBc(bool homogeneous_bc);

    /*!
     * \brief Return whether the solver is using homogeneous boundary
     * conditions.
     */
    virtual bool getHomogeneousBc() const;

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    virtual void setSolutionTime(double solution_time);

    /*!
     * \brief Get the time at which the solution is being evaluated.
     */
    virtual double getSolutionTime() const;

    /*!
     * \brief Set the current time interval.
     */
    virtual void setTimeInterval(double current_time, double new_time);

    /*!
     * \brief Get the current time interval.
     */
    virtual std::pair<double, double> getTimeInterval() const;

    /*!
     * \brief Get the current time step size.
     */
    virtual double getDt() const;

    /*!
     * \brief Set the HierarchyMathOps object used by the solver.
     */
    virtual void setHierarchyMathOps(SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops);

    /*!
     * \brief Get the HierarchyMathOps object used by the solver.
     */
    virtual SAMRAI::tbox::Pointer<HierarchyMathOps> getHierarchyMathOps() const;

    /*!
     * \brief Solve the system of equations.
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
    virtual bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                             SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b) = 0;

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
     */
    virtual void initializeSolverState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& x,
                                       const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& b);

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
     */
    virtual void deallocateSolverState();

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set the maximum number of nonlinear iterations to use per solve.
     */
    virtual void setMaxIterations(int max_iterations);

    /*!
     * \brief Get the maximum number of nonlinear iterations to use per solve.
     */
    virtual int getMaxIterations() const;

    /*!
     * \brief Set the absolute residual tolerance for convergence.
     */
    virtual void setAbsoluteTolerance(double abs_residual_tol);

    /*!
     * \brief Get the absolute residual tolerance for convergence.
     */
    virtual double getAbsoluteTolerance() const;

    /*!
     * \brief Set the relative residual tolerance for convergence.
     */
    virtual void setRelativeTolerance(double rel_residual_tol);

    /*!
     * \brief Get the relative residual tolerance for convergence.
     */
    virtual double getRelativeTolerance() const;

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent solve.
     */
    virtual int getNumIterations() const;

    /*!
     * \brief Return the residual norm from the most recent iteration.
     */
    virtual double getResidualNorm() const;

    //\}

    /*!
     * \name Logging functions.
     */
    //\{

    /*!
     * \brief Enable or disable logging.
     */
    virtual void setLoggingEnabled(bool enable_logging = true);

    /*!
     * \brief Determine whether logging is enabled or disabled.
     */
    virtual bool getLoggingEnabled() const;

    /*!
     * \brief Print class data to stream.
     */
    virtual void printClassData(std::ostream& stream);

    //\}

protected:
    // Basic initialization.
    void init(const std::string& object_name, bool homogeneous_bc);

    // Object name.
    std::string d_object_name;

    // Boolean value to indicate whether the preconditioner is presently
    // initialized.
    bool d_is_initialized;

    // Solver configuration.
    bool d_homogeneous_bc;
    double d_solution_time, d_current_time, d_new_time;
    double d_rel_residual_tol;
    double d_abs_residual_tol;
    int d_max_iterations;
    int d_current_iterations;
    double d_current_residual_norm;

    // Mathematical operators.
    SAMRAI::tbox::Pointer<HierarchyMathOps> d_hier_math_ops;
    bool d_hier_math_ops_external;

    // Logging configuration.
    bool d_enable_logging;

private:
    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    GeneralSolver(const GeneralSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    GeneralSolver& operator=(const GeneralSolver& that);
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_GeneralSolver
