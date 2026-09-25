// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2026 by the IBAMR developers
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

#ifndef included_IBTK_FACPreconditioner
#define included_IBTK_FACPreconditioner

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <ibtk/LinearSolver.h>
#include <ibtk/ibtk_enums.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <HierarchyDataOpsReal.h>
#include <IntVector.h>
#include <PatchHierarchy.h>
#include <SAMRAIVectorReal.h>

#include <string>
#include <vector>

namespace SAMRAI
{
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// FORWARD DECLARATIONS /////////////////////////

namespace IBTK
{
class FACPreconditionerStrategy;
} // namespace IBTK

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class FACPreconditioner is a concrete LinearSolver for implementing
 * FAC (multilevel multigrid) preconditioners.
 *
 * This class is similar to the SAMRAI class SAMRAI::solv::FACPreconditioner,
 * except that this class has been optimized for the case in which the solver is
 * to be used as a single-pass preconditioner, especially for the case in which
 * pre-smoothing sweeps are not needed.  This class is not suitable for use as a
 * stand-alone solver; rather, it is intended to be used in conjunction with an
 * iterative Krylov method.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 cycle_type = "V_CYCLE"  // see setMGCycleType()
 cycle_multiplicity = 2  // see setMGCycleMultiplicity()
 num_pre_sweeps = 0      // see setNumPreSmoothingSweeps()
 num_post_sweeps = 2     // see setNumPostSmoothingSweeps()
 enable_logging = FALSE  // see setLoggingEnabled()
 \endverbatim
*/
class FACPreconditioner : public LinearSolver
{
public:
    /*!
     * Constructor.
     */
    FACPreconditioner(std::string object_name,
                      SAMRAI::tbox::Pointer<FACPreconditionerStrategy> fac_strategy,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      const std::string& default_options_prefix);

    /*!
     * Destructor.
     */
    ~FACPreconditioner();

    /*!
     * \name Linear solver functionality.
     */
    //\{

    /*!
     * \brief Set whether the solver should use homogeneous boundary conditions.
     *
     * \note FACPreconditioner implements only the FAC correction scheme, which is
     * always subject to homogeneous boundary conditions; see
     * FACPreconditionerStrategy::ALWAYS_HOMOGENEOUS_BC. Calling this with \c
     * false triggers a TBOX_ERROR.
     */
    void setHomogeneousBc(bool homogeneous_bc) override;

    /*!
     * \brief Set the time at which the solution is to be evaluated.
     */
    void setSolutionTime(double solution_time) override;

    /*!
     * \brief Set the current time interval.
     */
    void setTimeInterval(double current_time, double new_time) override;

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
     * \note V-cycles without presmoothing restrict covered coarse RHS data in
     * place. Other cycle configurations preserve the RHS.
     *
     * \param x solution vector
     * \param b right-hand-side vector
     *
     * <b>Conditions on Parameters:</b>
     * - vectors \a x and \a b must have same patch hierarchy
     * - vectors \a x and \a b must have same structure, depth, etc.
     * - vectors \a x and \a b must not share component patch-data indices
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

    /*!
     * \brief Set the multigrid algorithm cycle type.
     *
     * V and W cycles make one and two recursive child visits, respectively.
     * MU_CYCLE uses setMGCycleMultiplicity(). F_CYCLE recursively makes an
     * F-cycle child visit followed by a V-cycle child visit. FMG_CYCLE performs
     * nested iteration with one V-cycle on each successively finer hierarchy.
     *
     * Except for V-cycles without presmoothing, this preconditioner currently
     * requires the initialized hierarchy range to begin at level zero, since
     * strategy residual ghost fills may access data below that range. Nonzero
     * coarsest levels are rejected for these cycles.
     */
    void setMGCycleType(MGCycleType cycle_type);

    /*!
     * \brief Get the multigrid algorithm cycle type.
     */
    MGCycleType getMGCycleType() const;

    /*!
     * \brief Set the positive number of child visits used by MU_CYCLE.
     *
     * The default is two. A multiplicity of one gives a V-cycle, and two gives
     * a W-cycle. This parameter does not change the other cycle types.
     */
    void setMGCycleMultiplicity(int cycle_multiplicity);

    /*!
     * \brief Get the number of child visits used by MU_CYCLE.
     */
    int getMGCycleMultiplicity() const;

    /*!
     * \brief Set the number of pre-smoothing sweeps to employ.
     */
    void setNumPreSmoothingSweeps(int num_pre_sweeps);

    /*!
     * \brief Get the number of pre-smoothing sweeps employed by the
     * preconditioner.
     */
    int getNumPreSmoothingSweeps() const;

    /*!
     * \brief Set the number of post-smoothing sweeps to employ.
     */
    void setNumPostSmoothingSweeps(int num_post_sweeps);

    /*!
     * \brief Get the number of post-smoothing sweeps employed by the
     * preconditioner.
     */
    int getNumPostSmoothingSweeps() const;

    //\}

    /*!
     * \brief Get the FAC preconditioner strategy objects employed by the
     * preconditioner.
     */
    SAMRAI::tbox::Pointer<FACPreconditionerStrategy> getFACPreconditionerStrategy() const;

protected:
    void FACVCycleNoPreSmoothing(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u,
                                 SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& f,
                                 int level_num);

    // Apply a cycle to an owned correction that is zero, including ghosts, on
    // the complete truncated hierarchy. The RHS is borrowed and preserved.
    void zeroStartCycle(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u,
                        SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& f,
                        int level_num,
                        MGCycleType cycle_type);

    // Improve a retained iterate by u += B(f - A*u), preserving its auxiliary
    // covered-cell values while evaluating the composite residual.
    void improveCycle(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& f,
                      int level_num,
                      MGCycleType cycle_type);

    void FMGCycle(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& u,
                  SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& f,
                  int level_num);

    SAMRAI::tbox::Pointer<FACPreconditionerStrategy> d_fac_strategy;
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM>> d_hierarchy;
    int d_coarsest_ln = 0;
    int d_finest_ln = 0;
    MGCycleType d_cycle_type = V_CYCLE;
    int d_cycle_multiplicity = 2;
    int d_num_pre_sweeps = 0, d_num_post_sweeps = 2;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    FACPreconditioner() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    FACPreconditioner(const FACPreconditioner& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    FACPreconditioner& operator=(const FACPreconditioner& that) = delete;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>
    getRangeVector(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& vector, int coarsest_ln, int finest_ln) const;

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>
    allocateRangeVector(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& vector,
                        const std::string& name,
                        int finest_ln) const;

    void allocateCycleScratchData(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                  const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs);

    void getFromInput(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

    // Vector arithmetic resets the operation object's level range. Keep these
    // objects independent of the caller, which may use its operations directly.
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM, double>>> d_vector_data_ops;

    // Each recursion depth owns its child equation and any warm-visit
    // increment. Storage can be reused only after that invocation returns.
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>> d_residual_vectors;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>> d_rhs_vectors;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>> d_correction_vectors;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>>> d_fmg_rhs_vectors;
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double>> d_evaluation_vector;
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_FACPreconditioner
