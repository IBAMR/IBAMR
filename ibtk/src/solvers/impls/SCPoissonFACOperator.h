// Filename: SCPoissonFACOperator.h
// Created on 13 Nov 2008 by Boyce Griffith
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

#ifndef included_SCPoissonFACOperator
#define included_SCPoissonFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

// PETSc INCLUDES
#include <petscmat.h>

// IBTK INCLUDES
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/CoarseFineBoundaryRefinePatchStrategy.h>
#include <ibtk/FACPreconditionerStrategy.h>
#include <ibtk/SCPoissonHypreLevelSolver.h>
#include <ibtk/SCPoissonPETScLevelSolver.h>
#include <ibtk/HierarchyMathOps.h>

// SAMRAI INCLUDES
#include <CoarsenAlgorithm.h>
#include <LocationIndexRobinBcCoefs.h>
#include <RefineAlgorithm.h>
#include <tbox/ConstPointer.h>

// C++ STDLIB INCLUDES
#include <map>

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class SCPoissonFACOperator is a concrete FACPreconditionerStrategy for
 * solving elliptic equations of the form \f$ \mbox{$L u$} = \mbox{$(C I +
 * \nabla \cdot D \nabla) u$} = f \f$ using a globally second-order accurate
 * side-centered finite-difference discretization, with support for Robin and
 * periodic boundary conditions.
 *
 * \warning This class was originally intended to be used with the SAMRAI class
 * SAMRAI::solv::FACPreconditioner but is now designed to be used with the IBTK
 * class IBTK::FACPreconditioner.
 *
 * This class provides operators that are used by class FACPreconditioner to
 * solve scalar Poisson-type equations of the form \f[ (C I + \nabla \cdot D
 * \nabla) u = f \f] using a side-centered, globally second-order accurate
 * finite-difference discretization, where
 *
 * - \f$ C \f$, \f$ D \f$ and \f$ f \f$ are independent of \f$ u \f$,
 * - \f$ C \f$ is a constant damping factor,
 * - \f$ D \f$ is a constant diffusion coefficient, and
 * - \f$ f \f$ is a side-centered scalar function.
 *
 * Robin boundary conditions may be specified at physical boundaries; see class
 * SAMRAI::solv::RobinBcCoefStrategy.
 *
 * By default, the class is configured to solve the Poisson problem \f$
 * -\nabla^2 u = f \f$, subject to homogeneous Dirichlet boundary conditions.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 smoother_choice = "additive"                 // see setSmootherChoice()

 prolongation_method = "CONSTANT_REFINE"      // see setProlongationMethod()
 restriction_method = "CONSERVATIVE_COARSEN"  // see setRestrictionMethod()

 coarse_solver_choice = "block_jacobi"        // see setCoarsestLevelSolverChoice()
 coarse_solver_tolerance = 1.0e-6             // see setCoarsestLevelSolverTolerance()
 coarse_solver_max_iterations = 10            // see setCoarsestLevelSolverMaxIterations()

 hypre_solver = { ... }                       // SAMRAI::tbox::Database for initializing class SCPoissonHypreLevelSolver

 petsc_solver = { ... }                       // SAMRAI::tbox::Database for initializing class SCPoissonPetscLevelSolver
 \endverbatim
*/
class SCPoissonFACOperator
    : public FACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    SCPoissonFACOperator(
        const std::string& object_name,
        const SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>& input_db=NULL);

    /*!
     * \brief Destructor.
     */
    ~SCPoissonFACOperator();

    /*!
     * \name Functions for specifying the Poisson problem.
     */
    //\{

    /*!
     * \brief Set the SAMRAI::solv::PoissonSpecifications object used to specify
     * the coefficients for the scalar Poisson equation.
     */
    void
    setPoissonSpecifications(
        const SAMRAI::solv::PoissonSpecifications& poisson_spec);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.
     *
     * \param bc_coefs  Vector of pointers to objects that can set the Robin boundary condition coefficients
     */
    void
    setPhysicalBcCoefs(
        const blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM>& bc_coefs);

    /*!
     * \brief Set the hierarchy time, for use with the refinement schedules and
     * boundary condition routines employed by the object.
     */
    void
    setTime(
        const double time);

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
     * - \c "multiplicative"
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
     * - \c "hypre"
     * - \c "petsc"
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
     * \note This value is ignored if the bottom solver choice is
     * "block_jacobi".
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
     * \brief Set the name of the prolongation method.
     */
    void
    setProlongationMethod(
        const std::string& prolongation_method);

    /*!
     * \brief Set the name of the restriction method.
     */
    void
    setRestrictionMethod(
        const std::string& restriction_method);

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
    void
    setFACPreconditioner(
        SAMRAI::tbox::ConstPointer<FACPreconditioner> preconditioner);

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
     * \brief Perform a given number of relaxations on the error.
     *
     * \param error error vector
     * \param residual residual vector
     * \param level_num level number
     * \param num_sweeps number of sweeps to perform
     * \param performing_pre_sweeps boolean value that is true when pre-smoothing sweeps are being performed
     * \param performing_post_sweeps boolean value that is true when post-smoothing sweeps are being performed
     */
    void
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
    bool
    solveCoarsestLevel(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& error,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        int coarsest_ln);

    /*!
     * \brief Compute composite grid residual on the specified range of levels.
     *
     * \param residual residual vector
     * \param solution solution vector
     * \param rhs source (right hand side) vector
     * \param coarsest_level_num coarsest level number
     * \param finest_level_num finest level number
     */
    void
    computeResidual(
        SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& residual,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& solution,
        const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& rhs,
        int coarsest_level_num,
        int finest_level_num);

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

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    SCPoissonFACOperator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    SCPoissonFACOperator(
        const SCPoissonFACOperator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    SCPoissonFACOperator& operator=(
        const SCPoissonFACOperator& that);

    /*!
     * \name For executing, caching and resetting communication schedules.
     */
    //\{

    /*!
     * \brief Execute a refinement schedule for prolonging data.
     */
    void
    xeqScheduleProlongation(
        const int dst_idx,
        const int src_idx,
        const int dst_ln);

    /*!
     * \brief Execute schedule for restricting solution or residual to the
     * specified level.
     */
    void
    xeqScheduleRestriction(
        const int dst_idx,
        const int src_idx,
        const int dst_ln);

    /*!
     * \brief Execute schedule for filling ghosts on the specified level.
     */
    void
    xeqScheduleGhostFillNoCoarse(
        const int dst_idx,
        const int dst_ln);

    /*!
     * \brief Execute schedule for synchronizing data on the specified level.
     */
    void
    xeqScheduleSideDataSynch(
        const int dst_idx,
        const int dst_ln);

    //\}

    /*!
     * \brief Initialize the hypre bottom solver.
     */
    void
    initializeHypreLevelSolver();

    /*!
     * \brief Initialize the PETSc bottom solver.
     */
    void
    initializePETScLevelSolver();

    /*!
     * \brief Construct a matrix corresponding to a Laplace operator restricted
     * to a single patch.
     */
    static void
    buildPatchLaplaceOperator(
        Mat& A,
        const SAMRAI::solv::PoissonSpecifications& poisson_spec,
        const SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> >& patch,
        const int component_axis,
        const SAMRAI::hier::IntVector<NDIM>& ghost_cell_width);

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
    int d_depth;

    /*
     * Mappings from patch indices to patch operators.
     */
    bool d_using_petsc_smoothers;
    SAMRAI::hier::IntVector<NDIM> d_gcw;
    std::vector<std::vector<blitz::TinyVector<Vec,NDIM> > > d_patch_vec_e, d_patch_vec_f;
    std::vector<std::vector<blitz::TinyVector<Mat,NDIM> > > d_patch_mat;
    std::vector<std::vector<blitz::TinyVector<SAMRAI::hier::BoxList<NDIM>,NDIM> > > d_patch_bc_box_overlap;
    std::vector<std::vector<blitz::TinyVector<std::map<int,SAMRAI::hier::Box<NDIM> >,NDIM> > > d_patch_smoother_bc_boxes;

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
    std::vector<SAMRAI::tbox::Pointer<HierarchyGhostCellInterpolation> > d_hier_bdry_fill_ops;
    std::vector<SAMRAI::tbox::Pointer<HierarchyMathOps> > d_hier_math_ops;

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
     * Scalar Poisson equations specifications.
     */
    SAMRAI::solv::PoissonSpecifications d_poisson_spec;

    /*
     * The kind of smoothing to perform.
     */
    std::string d_smoother_choice;

    /*
     * The name of the refinement operator used to prolong the coarse grid
     * correction and to set ghost cell values at coarse-fine interfaces.
     */
    std::string d_prolongation_method;

    /*
     * The name of the coarsening operator used to restrict the fine grid error
     * or residual.
     */
    std::string d_restriction_method;

    /*
     * Pointer to the FACPreconditioner that is using this operator.
     */
    SAMRAI::tbox::ConstPointer<FACPreconditioner> d_preconditioner;

    /*
     * Coarse level solver parameters.
     */
    std::string d_coarse_solver_choice;
    double d_coarse_solver_tol;
    int d_coarse_solver_max_its;
    bool d_using_hypre;
    SAMRAI::tbox::Pointer<SCPoissonHypreLevelSolver> d_hypre_solver;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_hypre_db;
    bool d_using_petsc;
    SAMRAI::tbox::Pointer<SCPoissonPETScLevelSolver> d_petsc_solver;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_petsc_db;

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
    int d_side_scratch_idx;

    //\}

    /*!
     * \name Boundary condition handling objects.
     */
    //\{

    SAMRAI::tbox::Pointer<CartSideRobinPhysBdryOp> d_bc_op;
    SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* const d_default_bc_coef;
    blitz::TinyVector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*,NDIM> d_bc_coefs;
    double d_apply_time;

    //\}

    /*!
     * \name Various refine and coarsen objects.
     */
    //\{

    /*
     * Coarse-fine interface interpolation object.
     */
    SAMRAI::tbox::Pointer<CoarseFineBoundaryRefinePatchStrategy> d_cf_bdry_op;

    /*
     * Variable fill pattern object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_op_stencil_fill_pattern, d_side_synch_fill_pattern;

    /*
     * Error prolongation (refinement) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_prolongation_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_prolongation_refine_patch_strategy;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_prolongation_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_prolongation_refine_schedules;

    /*
     * Residual restriction (coarsening) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_restriction_coarsen_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenAlgorithm<NDIM> > d_restriction_coarsen_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > > d_restriction_coarsen_schedules;

    /*
     * Refine operator for side data from same level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_ghostfill_nocoarse_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_nocoarse_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_nocoarse_refine_schedules;

    /*
     * Operator for side data synchronization on same level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_side_synch_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_side_synch_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_side_synch_refine_schedules;

    //\}
};
}// namespace IBTK

/////////////////////////////// INLINE ///////////////////////////////////////

//#include <ibtk/SCPoissonFACOperator.I>

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_SCPoissonFACOperator
