// Filename: INSStaggeredFACPreconditionerStrategy.h
// Created on 18 Apr 2012 by Boyce Griffith
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

#ifndef included_INSStaggeredFACPreconditionerStrategy
#define included_INSStaggeredFACPreconditionerStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscksp.h>

// IBTK INCLUDES
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/CoarseFineBoundaryRefinePatchStrategy.h>
#include <ibtk/FACPreconditionerStrategy.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>

// SAMRAI INCLUDES
#include <CoarsenAlgorithm.h>
#include <RefineAlgorithm.h>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredFACPreconditionerStrategy is an abstract
 * FACPreconditionerStrategy implementing many of the operations required by
 * smoothers for staggered-grid (MAC) discretizations of the incompressible
 * Stokes equations and related problems.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 smoother_choice = "additive"                   // see setSmootherChoice()

 U_prolongation_method = "CONSTANT_REFINE"      // see setProlongationMethods()
 P_prolongation_method = "LINEAR_REFINE"        // see setProlongationMethods()
 U_restriction_method = "CONSERVATIVE_COARSEN"  // see setRestrictionMethods()
 P_restriction_method = "CONSERVATIVE_COARSEN"  // see setRestrictionMethods()

 coarse_solver_choice = "block_jacobi"          // see setCoarsestLevelSolverChoice()
 coarse_solver_tolerance = 1.0e-6               // see setCoarsestLevelSolverTolerance()
 coarse_solver_max_iterations = 10              // see setCoarsestLevelSolverMaxIterations()
 \endverbatim
*/
class INSStaggeredFACPreconditionerStrategy
    : public IBTK::FACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    INSStaggeredFACPreconditionerStrategy(
        const std::string& object_name,
        int ghost_cell_width,
        SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db);

    /*!
     * \brief Destructor.
     */
    ~INSStaggeredFACPreconditionerStrategy();

    /*!
     * \brief Set the current time interval, for use with the refinement
     * schedules and boundary condition routines employed by the object.
     */
    void
    setTimeInterval(
        double current_time,
        double new_time);

    /*!
     * \name Functions for configuring the solver.
     */
    //\{

    /*!
     * \brief Specify the levels that need to be reset the next time the
     * operator is re-initialized.
     *
     * When the operator is initialized, then only the specified range of levels
     * are reset in the operator state the next time that the operator is
     * initialized.  If the operator is not initialized, this method has no
     * effect.
     *
     * To ensure the range of levels that is reset includes all levels in the
     * patch hierarchy, use \a coarsest_ln = \a finest_ln = \p -1.
     *
     * \note This function is used to save some unnecessary computations when
     * the hierarchy is regridded.  The range of levels specified must include
     * all levels which need to be reset by
     * SAMRAI::mesh::StandardTagAndInitStrategy::resetHierarchyConfiguration().
     * Any data residing outside of this range of levels will not be reset.
     * This \b is \b not what you want to have happen if, for instance, the
     * Poisson specifications changes.
     */
    void
    setResetLevels(
        int coarsest_ln,
        int finest_ln);

    /*!
     * \brief Specify the smoother type.
     */
    virtual void
    setSmootherChoice(
        const std::string& smoother_choice) = 0;

    /*!
     * \brief Specify the coarse level solver.
     */
    virtual void
    setCoarsestLevelSolverChoice(
        const std::string& coarse_solver_choice) = 0;

    /*!
     * \brief Set tolerance for coarse level solve.
     *
     * If the coarse level solver uses a convergence tolerance parameter, the
     * specified value is used.  If the coarse level solver does not use such a
     * stopping parameter, implementations are free to ignore this value.
     */
    void
    setCoarsestLevelSolverTolerance(
        double coarse_solver_tol);

    /*!
     * \brief Set the maximum number of iterations for the coarsest level solve.
     *
     * If the coarse level solver uses a maximum number of iterations parameter,
     * the specified value is used.  If the coarse level solver does not use
     * such a stopping parameter, implementations are free to ignore this value.
     */
    void
    setCoarsestLevelSolverMaxIterations(
        int coarse_solver_max_its);

    /*!
     * \brief Set the prolongation methods.
     */
    void
    setProlongationMethods(
        const std::string& U_prolongation_method,
        const std::string& P_prolongation_method);

    /*!
     * \brief Set the restriction methods.
     */
    void
    setRestrictionMethods(
        const std::string& U_restriction_method,
        const std::string& P_restriction_method);

    //\}

    /*!
     * \name Partial implementation of FACPreconditionerStrategy interface.
     */
    //\{

    /*!
     * \brief Restrict the residual quantity to the specified level from the
     * next finer level.
     *
     * \param src source residual
     * \param dst destination residual
     * \param dst_ln destination level number
     */
    void
    restrictResidual(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& src,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dst,
        int dst_ln);

    /*!
     * \brief Prolong the error quantity to the specified level from the next
     * coarser level.
     *
     * \param src source error vector
     * \param dst destination error vector
     * \param dst_ln destination level number of data transfer
     */
    void
    prolongError(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& src,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dst,
        int dst_ln);

    /*!
     * \brief Prolong the error quantity to the specified level from the next
     * coarser level and apply the correction to the fine-level error.
     *
     * \param src source error vector
     * \param dst destination error vector
     * \param dst_ln destination level number of data transfer
     */
    void
    prolongErrorAndCorrect(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& src,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dst,
        int dst_ln);

    /*!
     * \brief Compute hierarchy-dependent data.
     *
     * Note that although the vector arguments given to other methods in this
     * class may not necessarily be the same as those given to this method,
     * there will be similarities, including:
     *
     * - hierarchy configuration (hierarchy pointer and level range)
     * - number, type and alignment of vector component data
     * - ghost cell width of data in the solution (or solution-like) vector
     *
     * \param solution solution vector u
     * \param rhs right hand side vector f
     */
    void
    initializeOperatorState(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs);

    /*!
     * \brief Remove all hierarchy-dependent data.
     *
     * Remove all hierarchy-dependent data set by initializeOperatorState().
     *
     * \see initializeOperatorState
     */
    void
    deallocateOperatorState();

    //\}

protected:
    /*!
     * \brief Compute implementation-specific hierarchy-dependent data.
     */
    virtual void
    initializeOperatorStateSpecialized(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs,
        int coarsest_reset_ln,
        int finest_reset_ln) = 0;

    /*!
     * \brief Remove implementation-specific hierarchy-dependent data.
     */
    virtual void
    deallocateOperatorStateSpecialized(
        int coarsest_reset_ln,
        int finest_reset_ln) = 0;

    /*!
     * \name Methods for executing, caching, and resetting communication
     * schedules.
     */
    //\{

    /*!
     * \brief Execute a refinement schedule for prolonging data.
     */
    void
    xeqScheduleProlongation(
        const std::pair<int,int>& dst_idxs,
        const std::pair<int,int>& src_idxs,
        int dst_ln);

    /*!
     * \brief Execute schedule for restricting solution or residual to the
     * specified level.
     */
    void
    xeqScheduleRestriction(
        const std::pair<int,int>& dst_idxs,
        const std::pair<int,int>& src_idxs,
        int dst_ln);

    /*!
     * \brief Execute schedule for filling ghosts on the specified level.
     */
    void
    xeqScheduleGhostFillNoCoarse(
        const std::pair<int,int>& dst_idxs,
        int dst_ln);

    /*!
     * \brief Execute schedule for synchronizing data on the specified level.
     */
    void
    xeqScheduleDataSynch(
        int dst_idx,
        int dst_ln);

    //\}

    /*
     * The object name is used for error reporting purposes.
     *
     * The boolean indicates whether this object has been initialized.
     */
    std::string d_object_name;
    bool d_is_initialized;

    /*
     * The current time interval.
     */
    double d_current_time, d_new_time;

    /*
     * Ghost cell width.
     */
    const SAMRAI::hier::IntVector<NDIM> d_gcw;

    /*!
     * \name Hierarchy-dependent objects.
     */
    //\{

    /*
     * Solution and rhs vectors.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_solution, d_rhs;

    /*
     * Reference patch hierarchy and range of levels involved in the solve.
     *
     * This variable is non-null between the initializeOperatorState() and
     * deallocateOperatorState() calls.  It is not truly needed, because the
     * hierarchy is obtainable through variables in most function argument
     * lists.  We use it to enforce working on one hierarchy at a time.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    /*
     * Level operators, used to compute composite-grid residuals.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> > d_hier_bdry_fill_ops;
    std::vector<SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> > d_hier_math_ops;

    /*
     * Range of levels to be reset the next time the operator is initialized.
     */
    bool d_in_initialize_operator_state;
    int d_coarsest_reset_ln, d_finest_reset_ln;

    //\}

    /*!
     * \name Solver configuration variables.
     */
    //\{

    /*
     * The kind of smoothing to perform.
     */
    std::string d_smoother_choice;

    /*
     * The names of the refinement operators used to prolong the coarse grid
     * correction.
     */
    std::string d_U_prolongation_method, d_P_prolongation_method;

    /*
     * The names of the coarsening operators used to restrict the fine grid
     * error or residual.
     */
    std::string d_U_restriction_method, d_P_restriction_method;

    /*
     * Coarse level solver parameters.
     */
    std::string d_coarse_solver_choice;
    double d_coarse_solver_tol;
    int d_coarse_solver_max_its;

    //\}

    /*!
     * \name Internal context and scratch data.
     */
    //\{

    /*
     * Variable context for internally maintained hierarchy data.
     */
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> d_context;

    /*
     * Patch descriptor indices for scratch data.
     */
    int d_side_scratch_idx, d_cell_scratch_idx;

    //\}

    /*!
     * \name Various refine and coarsen objects.
     */
    //\{

    /*
     * Physical boundary operators.
     */
    SAMRAI::tbox::Pointer<IBTK::CartSideRobinPhysBdryOp> d_U_bc_op;
    SAMRAI::tbox::Pointer<IBTK::CartCellRobinPhysBdryOp> d_P_bc_op;

    /*
     * Coarse-fine interface interpolation objects.
     */
    SAMRAI::tbox::Pointer<IBTK::CoarseFineBoundaryRefinePatchStrategy> d_U_cf_bdry_op, d_P_cf_bdry_op;

    /*
     * Variable fill pattern object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_U_op_stencil_fill_pattern,  d_P_op_stencil_fill_pattern, d_U_synch_fill_pattern;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredFACPreconditionerStrategy();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredFACPreconditionerStrategy(
        const INSStaggeredFACPreconditionerStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredFACPreconditionerStrategy& operator=(
        const INSStaggeredFACPreconditionerStrategy& that);

    /*
     * Combined U & P physical boundary operator.
     */
    SAMRAI::xfer::RefinePatchStrategy<NDIM>* d_U_P_bc_op;

    /*
     * Error prolongation (refinement) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_U_prolongation_refine_operator, d_P_prolongation_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_prolongation_refine_patch_strategy;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_prolongation_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_prolongation_refine_schedules;

    /*
     * Residual restriction (coarsening) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_U_restriction_coarsen_operator, d_P_restriction_coarsen_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_restriction_coarsen_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_restriction_coarsen_schedules;

    /*
     * Refine operator for side and cell data from same level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_nocoarse_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_nocoarse_refine_schedules;

    /*
     * Operator for side data synchronization on same level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_synch_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_synch_refine_schedules;
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/INSStaggeredFACPreconditionerStrategy.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredFACPreconditionerStrategy
