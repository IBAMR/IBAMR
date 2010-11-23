// Filename: IBImplicitModHelmholtzPETScLevelSolver.h
// Created on 27 Sep 2010 by Boyce Griffith
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

#ifndef included_IBImplicitModHelmholtzPETScLevelSolver
#define included_IBImplicitModHelmholtzPETScLevelSolver

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscksp.h>

// IBTK INCLUDES
#include <ibtk/LinearSolver.h>

// SAMRAI INCLUDES
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <RefineSchedule.h>

// C++ STDLIB INCLUDES
#include <vector>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBImplicitModHelmholtzPETScLevelSolver is a concrete
 * IBTK::LinearSolver for solving equations related to an implicit
 * staggered-grid (MAC) discretization of the IB method.
 */
class IBImplicitModHelmholtzPETScLevelSolver
    : public IBTK::LinearSolver
{
public:
    /*!
     * \brief Constructor.
     *
     * \param object_name  Name of object.
     * \param input_db     Optional SAMRAI::tbox::Database for input.
     */
    IBImplicitModHelmholtzPETScLevelSolver(
        const std::string& object_name,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db=NULL);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~IBImplicitModHelmholtzPETScLevelSolver();

    /*!
     * \brief Set the options prefix used by this PETSc solver object.
     */
    void
    setOptionsPrefix(
        const std::string options_prefix);

    /*!
     * \name Functions for specifying the linear system.
     */
    //\{

    /*!
     * \brief Set the S dF/dX R = S dF/dX S^{*} matrix.
     */
    void
    setSJRMat(
        Mat& SJR_mat);

    /*!
     * \brief Set the scalar Poisson equation specifications.
     */
    void
    setPoissonSpecifications(
        const SAMRAI::solv::PoissonSpecifications& poisson_spec);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy object used to specify
     * physical boundary conditions.
     *
     * \note \a bc_coef may be NULL.  In this case, homogeneous Dirichlet
     * boundary conditions are employed.
     *
     * \param bc_coef  Pointer to an object that can set the Robin boundary condition coefficients
     */
    virtual void
    setPhysicalBcCoef(
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const bc_coef);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy object used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.
     *
     * \param bc_coefs  Vector of pointers to objects that can set the Robin boundary condition coefficients
     */
    virtual void
    setPhysicalBcCoefs(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

    /*!
     * \brief Specify whether the boundary conditions are homogeneous.
     */
    virtual void
    setHomogeneousBc(
        const bool homogeneous_bc);

    /*!
     * \brief Set the hierarchy time, for use with the refinement schedules and
     * boundary condition routines employed by the object.
     */
    void
    setTime(
        const double time);

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
    virtual bool
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
    virtual void
    initializeSolverState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b);

    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * \note It is safe to call deallocateSolverState() when the solver state is
     * already deallocated.
     *
     * \see initializeSolverState
     */
    virtual void
    deallocateSolverState();

    //\}

    /*!
     * \name Functions to access solver parameters.
     */
    //\{

    /*!
     * \brief Set whether the initial guess is non-zero.
     */
    virtual void
    setInitialGuessNonzero(
        bool initial_guess_nonzero=true);

    /*!
     * \brief Get whether the initial guess is non-zero.
     */
    virtual bool
    getInitialGuessNonzero() const;

    /*!
     * \brief Set the maximum number of iterations to use per solve.
     */
    virtual void
    setMaxIterations(
        int max_iterations);

    /*!
     * \brief Get the maximum number of iterations to use per solve.
     */
    virtual int
    getMaxIterations() const;

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

    //\}

    /*!
     * \name Functions to access data on the most recent solve.
     */
    //\{

    /*!
     * \brief Return the iteration count from the most recent linear solve.
     */
    virtual int
    getNumIterations() const;

    /*!
     * \brief Return the residual norm from the most recent iteration.
     */
    virtual double
    getResidualNorm() const;

    //\}

    /*!
     * \name Logging functions.
     */
    //\{

    /*!
     * \brief Enable or disable logging.
     */
    virtual void
    enableLogging(
        bool enabled=true);

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBImplicitModHelmholtzPETScLevelSolver();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBImplicitModHelmholtzPETScLevelSolver(
        const IBImplicitModHelmholtzPETScLevelSolver& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBImplicitModHelmholtzPETScLevelSolver&
    operator=(
        const IBImplicitModHelmholtzPETScLevelSolver& that);

    /*!
     * \brief Object name.
     */
    std::string d_object_name;

    /*!
     * \brief Solver initialization status.
     */
    bool d_is_initialized;

    /*!
     * \brief Associated hierarchy.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;

    /*!
     * \brief Associated level number.
     */
    int d_level_num;

    /*!
     * \name Problem specification and boundary condition handling.
     */
    //\{
    Mat d_SJR_mat;
    SAMRAI::solv::PoissonSpecifications d_poisson_spec;

    /*!
     * \brief Robin boundary coefficient object for physical boundaries and
     * related data.
     */
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* const d_default_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;
    bool d_homogeneous_bc;
    double d_apply_time;

    //\}

    /*!
     * \name PETSc objects.
     */
    //\{
    std::string d_options_prefix;
    KSP d_petsc_ksp;
    Mat d_petsc_mat;
    Vec d_petsc_x, d_petsc_b;
    int d_max_iterations;
    double d_abs_residual_tol;
    double d_rel_residual_tol;
    bool d_initial_guess_nonzero;

    int d_current_its;
    double d_current_residual_norm;

    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;
    int d_dof_index_idx;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,int> > d_dof_index_var;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > d_dof_index_fill;

    //\}

    /*!
     * \name Variables for debugging and analysis.
     */
    //\{

    /*!
     * \brief Flag to print solver info.
     */
    bool d_enable_logging;

    //\}
};
}// namespace IBAMR

/////////////////////////////// INLINE ///////////////////////////////////////

#include <ibamr/IBImplicitModHelmholtzPETScLevelSolver.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBImplicitModHelmholtzPETScLevelSolver
