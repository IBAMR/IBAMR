// Filename: INSStaggeredBoxRelaxationFACOperator.h
// Created on 11 Jun 2010 by Boyce Griffith
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

#ifndef included_INSStaggeredBoxRelaxationFACOperator
#define included_INSStaggeredBoxRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscksp.h>

// IBAMR INCLUDES
#include <ibamr/INSCoefs.h>

// IBTK INCLUDES
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/CoarseFineBoundaryRefinePatchStrategy.h>
#include <ibtk/FACPreconditionerStrategy.h>

// SAMRAI INCLUDES
#include <CoarsenAlgorithm.h>
#include <LocationIndexRobinBcCoefs.h>
#include <RefineAlgorithm.h>
#include <tbox/ConstPointer.h>

// C++ STDLIB INCLUDES
#include <map>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class INSStaggeredBoxRelaxationFACOperator is a concrete
 * FACPreconditionerStrategy implementing a box relaxation (Vanka-type) smoother
 * for use as a multigrid preconditioner.
 */
class INSStaggeredBoxRelaxationFACOperator
    : public virtual IBTK::FACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    INSStaggeredBoxRelaxationFACOperator(
        const std::string& object_name,
        const INSCoefs& problem_coefs,
        const double dt,
        const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db=NULL);

    /*!
     * \brief Virtual destructor.
     */
    virtual
    ~INSStaggeredBoxRelaxationFACOperator();

    /*!
     * \name Functions for specifying the problem coefficients.
     */
    //\{

    /*!
     * \brief Set the INSCoefs object and timestep size used to specify the
     * coefficients for the time-dependent incompressible Stokes operator.
     */
    void
    setProblemCoefficients(
        const INSCoefs& problem_coefs,
        const double dt);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a U_bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be NULL; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param U_bc_coefs  Vector of pointers to objects that can set the Robin boundary condition coefficients
     * \param P_bc_coef   Pointers to object that can set the Robin boundary condition coefficients
     */
    void
    setPhysicalBcCoefs(
        const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
        SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef);

    /*!
     * \brief Specify whether the boundary conditions are homogeneous.
     */
    void
    setHomogeneousBc(
        const bool homogeneous_bc);

    /*!
     * \brief Set the current time interval.
     */
    void
    setTimeInterval(
        const double current_time,
        const double new_time);

    //\}

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
        const int coarsest_ln,
        const int finest_ln);

    /*!
     * \brief Specify the ghost cell width for \em both the solution and the
     * right-hand side patch data.
     */
    void
    setGhostCellWidth(
        const SAMRAI::hier::IntVector<NDIM>& ghost_cell_width);

    /*!
     * \brief Specify the smoother type.
     *
     * Select from:
     * - \c "additive"
     *
     * \note The smoother is always additive between processors ("processor
     * block Gauss-Seidel").
     */
    void
    setSmootherChoice(
        const std::string& smoother_choice);

    /*!
     * \brief Specify the coarse level solver.
     *
     * Select from:
     * - \c "block_jacobi"
     */
    void
    setCoarsestLevelSolverChoice(
        const std::string& coarse_solver_choice);

    /*!
     * \brief Set tolerance for coarse level solve.
     *
     * If the coarse level solver requires a tolerance, the specified value is
     * used.
     *
     * \note Currently, this value is not used.
     */
    void
    setCoarsestLevelSolverTolerance(
        double coarse_solver_tol);

    /*!
     * \brief Set the maximum number of iterations for the coarsest level solve.
     */
    void
    setCoarsestLevelSolverMaxIterations(
        int coarse_solver_max_its);

    /*!
     * \brief Set the name of the prolongation methods.
     */
    void
    setProlongationMethods(
        const std::string& U_prolongation_method,
        const std::string& P_prolongation_method);

    /*!
     * \brief Set the name of the restriction methods.
     */
    void
    setRestrictionMethods(
        const std::string& U_restriction_method,
        const std::string& P_restriction_method);

    /*!
     * \brief Set the maximum number of FAC cycles used by the FACPreconditioner
     * that employs this concrete FACPreconditionerStrategy.
     */
    void
    setFACPreconditionerMaxCycles(
        int fac_max_cycles);

    /*!
     * \brief Set whether the FACPreconditioner that employs this concrete
     * FACPreconditionerStrategy is employing presmoothing.
     */
    void
    setFACPreconditionerUsesPresmoothing(
        bool fac_uses_presmoothing);

    /*!
     * \brief Set whether the FACPreconditioner that employs this concrete
     * FACPreconditionerStrategy uses a nonzero initial guess.
     */
    void
    setFACPreconditionerInitialGuessNonzero(
        bool fac_initial_guess_nonzero);

    /*!
     * \brief Set whether to skip restricting the solution.
     *
     * If this operator is being used as a preconditioner, it may be possible to
     * avoid restricting the solution.
     */
    void
    setSkipRestrictSolution(
        bool skip_restrict_sol);

    /*!
     * \brief Set whether to skip restricting the residual.
     *
     * If this operator is being used as a preconditioner, it may be possible to
     * avoid restricting the residual.
     */
    void
    setSkipRestrictResidual(
        bool skip_restrict_residual);

    //\}

    ///
    ///  The following routines:
    ///
    ///      setFACPreconditioner(),
    ///      restrictResidual(),
    ///      prolongError(),
    ///      prolongErrorAndCorrect(),
    ///      smoothError(),
    ///      solveCoarsestLevel(),
    ///      computeResidual(),
    ///      initializeOperatorState(),
    ///      deallocateOperatorState()
    ///
    ///  are concrete implementations of functions declared in the
    ///  FACPreconditionerStrategy abstract base class.
    ///

    /*!
     * \name Implementation of FACPreconditionerStrategy interface.
     */
    //\{

    /*!
     * \brief Set the FACPreconditioner object that is using this concrete
     * FACPreconditionerStrategy object.
     *
     * \param preconditioner  Pointer to the FAC preconditioner that is using this concrete FAC strategy
     */
    virtual void
    setFACPreconditioner(
        SAMRAI::tbox::ConstPointer<IBTK::FACPreconditioner> preconditioner);

    /*!
     * \brief Restrict the residual quantity to the specified level from the
     * next finer level.
     *
     * \param src source residual
     * \param dst destination residual
     * \param dst_ln destination level number
     */
    virtual void
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
    virtual void
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
    virtual void
    prolongErrorAndCorrect(
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& src,
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& dst,
        int dst_ln);

    /*!
     * \brief Perform a given number of relaxations on the error.
     *
     * \param error error vector
     * \param residual residual vector
     * \param level_num level number
     * \param num_sweeps number of sweeps to perform
     * \param performing_pre_sweeps boolean value that is true when pre-smoothing sweeps are being performed
     * \param performing_post_sweeps boolean value that is true when post-smoothing sweeps are being performed
     */
    virtual void
    smoothError(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& error,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        int level_num,
        int num_sweeps,
        bool performing_pre_sweeps,
        bool performing_post_sweeps);

    /*!
     * \brief Solve the residual equation Ae=r on the coarsest level of the
     * patch hierarchy.
     *
     * \param error error vector
     * \param residual residual vector
     * \param coarsest_ln coarsest level number
     */
    virtual bool
    solveCoarsestLevel(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& error,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        int coarsest_ln);

    /*!
     * \brief Compute composite grid residual on a single level.
     *
     * \param residual residual vector
     * \param solution solution vector
     * \param rhs source (right hand side) vector
     * \param level_num level number
     * \param error_equation_indicator flag stating whether u is an error vector or a solution vector
     */
    virtual void
    computeResidual(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs,
        int level_num);

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
    virtual void
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
    virtual void
    deallocateOperatorState();

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    INSStaggeredBoxRelaxationFACOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    INSStaggeredBoxRelaxationFACOperator(
        const INSStaggeredBoxRelaxationFACOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    INSStaggeredBoxRelaxationFACOperator& operator=(
        const INSStaggeredBoxRelaxationFACOperator& that);

    /*!
     * \name For executing, caching and resetting communication schedules.
     */
    //\{

    /*!
     * \brief Execute a refinement schedule for prolonging data.
     */
    void
    xeqScheduleProlongation(
        const std::pair<int,int>& dst_idxs,
        const std::pair<int,int>& src_idxs,
        const int dst_ln,
        const bool homogeneous_bc);

    /*!
     * \brief Execute schedule for restricting solution to the specified level.
     */
    void
    xeqScheduleURestriction(
        const std::pair<int,int>& dst_idxs,
        const std::pair<int,int>& src_idxs,
        const int dst_ln);

    /*!
     * \brief Execute schedule for restricting residual to the specified level.
     */
    void
    xeqScheduleRRestriction(
        const std::pair<int,int>& dst_idxs,
        const std::pair<int,int>& src_idxs,
        const int dst_ln);

    /*!
     * \brief Execute schedule for filling ghosts on the specified level.
     */
    void
    xeqScheduleGhostFill(
        const std::pair<int,int>& dst_idxs,
        const int dst_ln,
        const bool homogeneous_bc);

    /*!
     * \brief Execute schedule for filling ghosts on the specified level.
     */
    void
    xeqScheduleGhostFillNoCoarse(
        const std::pair<int,int>& dst_idxs,
        const int dst_ln,
        const bool homogeneous_bc);

    /*!
     * \brief Execute schedule for synchronizing data on the specified level.
     */
    void
    xeqScheduleSideDataSynch(
        const int U_dst_idx,
        const int dst_ln);

    //\}

    /*!
     * \brief Check to make sure that all of the options make sense.
     */
    void
    sanityCheck();

    /*
     * The object name is used for error reporting purposes.
     *
     * The boolean indicates whether this object has been initialized.
     */
    std::string d_object_name;
    bool d_is_initialized;

    /*!
     * \name Hierarchy-dependent objects.
     */
    //\{

    /*
     * Solution and rhs vectors.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > d_solution, d_rhs;

    /*
     * Ghost cell width.
     */
    SAMRAI::hier::IntVector<NDIM> d_gcw;

    /*
     * Box operator data.
     */
    std::vector<Mat> d_box_op;
    std::vector<Vec> d_box_e, d_box_r;
    std::vector<KSP> d_box_ksp;

    /*
     * Mappings from patch indices to patch operators.
     */
    std::vector<std::map<int,std::vector<SAMRAI::hier::BoxList<NDIM> > > > d_patch_side_bc_box_overlap;
    std::vector<std::map<int,SAMRAI::hier::BoxList<NDIM> > > d_patch_cell_bc_box_overlap;

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
     * Range of levels to be reset the next time the operator is initialized.
     */
    bool d_in_initialize_operator_state;
    int d_coarsest_reset_ln, d_finest_reset_ln;

    //\}

    /*!
     * \name Private state variables for solution process.
     */
    //\{

    /*
     * Problem coefficient specifications.
     */
    INSCoefs d_problem_coefs;
    double d_dt;

    /*
     * The kind of smoothing to perform.
     */
    std::string d_smoother_choice;

    /*
     * The name of the refinement operators used to prolong the coarse grid
     * correction.
     */
    std::string d_U_prolongation_method, d_P_prolongation_method;

    /*
     * The name of the coarsening operators used to restrict the fine grid error
     * or residual.
     */
    std::string d_U_restriction_method, d_P_restriction_method;

    /*
     * Pointer to the FACPreconditioner that is using this operator.
     *
     * Integer indicates the maximum number of FAC cycles.  When this number is
     * 1, some computations and parallel communication can be avoided.
     *
     * Booleans that indicate whether the FAC preconditioner associated with
     * this concrete FACPreconditionerStrategy is performing presmoothing or
     * using a nonzero initial guess.
     *
     * If the FACPreconditioner is:
     *
     *     1) being used as a preconditioner for a Krylov subspace linear
     *        solver, and
     *     2) not performing presmoothing (which doesn't seem to work anyway)
     *
     * then we can save a substantial amount of work.
     */
    SAMRAI::tbox::ConstPointer<IBTK::FACPreconditioner> d_preconditioner;
    int d_fac_max_cycles;
    bool d_fac_uses_presmoothing, d_fac_initial_guess_nonzero;

    /*
     * In some cases, it may be desirable to skip restricting the solution or
     * the residual.
     */
    bool d_skip_restrict_sol, d_skip_restrict_residual;

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
     * \name Boundary condition handling objects.
     */
    //\{

    SAMRAI::tbox::Pointer<IBTK::CartSideRobinPhysBdryOp> d_U_bc_op;
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* const d_default_U_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_U_bc_coefs;

    SAMRAI::tbox::Pointer<IBTK::CartCellRobinPhysBdryOp> d_P_bc_op;
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* const d_default_P_bc_coef;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef;

    SAMRAI::xfer::RefinePatchStrategy<NDIM>* d_U_P_bc_op;

    bool d_homogeneous_bc;
    double d_current_time, d_new_time;

    //\}

    /*!
     * \name Various refine and coarsen objects.
     */
    //\{

    /*
     * Variable fill pattern object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_U_op_stencil_fill_pattern,  d_P_op_stencil_fill_pattern, d_U_synch_fill_pattern;

    /*
     * Error prolongation (refinement) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_U_prolongation_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_P_prolongation_refine_operator;
    SAMRAI::tbox::Pointer<IBTK::CoarseFineBoundaryRefinePatchStrategy> d_U_cf_bdry_op;
    SAMRAI::tbox::Pointer<IBTK::CoarseFineBoundaryRefinePatchStrategy> d_P_cf_bdry_op;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_prolongation_refine_patch_strategy;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_prolongation_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_prolongation_refine_schedules;

    /*
     * Solution restriction (coarsening) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_U_urestriction_coarsen_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_P_urestriction_coarsen_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_urestriction_coarsen_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_urestriction_coarsen_schedules;

    /*
     * Residual restriction (coarsening) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_U_rrestriction_coarsen_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_P_rrestriction_coarsen_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_rrestriction_coarsen_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_rrestriction_coarsen_schedules;

    /*
     * Refine operator for side data from coarser level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_U_ghostfill_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_P_ghostfill_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_refine_schedules;

    /*
     * Refine operator for side data from same level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_U_ghostfill_nocoarse_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_P_ghostfill_nocoarse_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_nocoarse_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_nocoarse_refine_schedules;

    /*
     * Operator for side data synchronization on same level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_U_synch_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_U_synch_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_U_synch_refine_schedules;

    //\}
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/INSStaggeredBoxRelaxationFACOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_INSStaggeredBoxRelaxationFACOperator
