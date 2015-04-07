// Filename: StaggeredStokesIBLevelRelaxationFACOperator.h
// Created on 22 Mar 2015 by Amneet Bhalla
//
// Copyright (c) 2002-2015, Amneet Bhalla and Boyce Griffith
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

#ifndef included_StaggeredStokesIBLevelRelaxationFACOperator
#define included_StaggeredStokesIBLevelRelaxationFACOperator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <utility>
#include <vector>

#include "CoarsenAlgorithm.h"
#include "CoarsenOperator.h"
#include "IntVector.h"
#include "PatchHierarchy.h"
#include "PoissonSpecifications.h"
#include "RefineAlgorithm.h"
#include "RefineOperator.h"
#include "RefinePatchStrategy.h"
#include "SAMRAIVectorReal.h"
#include "VariableContext.h"
#include "VariableFillPattern.h"
#include "ibamr/StaggeredStokesFACPreconditioner.h"
#include "ibamr/StaggeredStokesPhysicalBoundaryHelper.h"
#include "ibtk/CartCellRobinPhysBdryOp.h"
#include "ibtk/CartSideRobinPhysBdryOp.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"
#include "petscao.h"
#include "petscmat.h"

namespace boost
{
template <class T, std::size_t N>
class array;
} // namespace boost
namespace IBAMR
{
class StaggeredStokesPETScLevelSolver;
}// namespace IBAMR
namespace IBTK
{
class HierarchyGhostCellInterpolation;
class HierarchyMathOps;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class BoxList;
} // namespace hier
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class StaggeredStokesIBLevelRelaxationFACOperator is a concrete
 * FACPreconditionerStrategy class that implements the operations required by
 * smoothers for staggered-grid (MAC) discretizations of the implicit 
 * incompressible Stokes-IB equations.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 smoother_type = "ADDITIVE"                     // see setSmootherType()
 U_prolongation_method = "CONSTANT_REFINE"      // see setProlongationMethods()
 P_prolongation_method = "LINEAR_REFINE"        // see setProlongationMethods()
 U_restriction_method = "CONSERVATIVE_COARSEN"  // see setRestrictionMethods()
 P_restriction_method = "CONSERVATIVE_COARSEN"  // see setRestrictionMethods()
 coarse_solver_type = "BLOCK_JACOBI"            // see setCoarseSolverType()
 coarse_solver_rel_residual_tol = 1.0e-5        // see setCoarseSolverRelativeTolerance()
 coarse_solver_abs_residual_tol = 1.0e-50       // see setCoarseSolverAbsoluteTolerance()
 coarse_solver_max_iterations = 10              // see setCoarseSolverMaxIterations()
 coarse_solver_db = { ... }                     // SAMRAI::tbox::Database for initializing
 coarse
 level solver
 \endverbatim
*/
class StaggeredStokesIBLevelRelaxationFACOperator : public IBTK::FACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    StaggeredStokesIBLevelRelaxationFACOperator(const std::string& object_name,
                                             SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                                             const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~StaggeredStokesIBLevelRelaxationFACOperator();

	/*!
	 * \brief Static function to construct a StaggeredStokesFACPreconditioner with a
	 * StaggeredStokesIBLevelRelaxationFACOperator FAC strategy.
	 */
	static SAMRAI::tbox::Pointer<StaggeredStokesSolver> allocate_solver(const std::string& object_name,
																		SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
																		const std::string& default_options_prefix)
	{
		SAMRAI::tbox::Pointer<IBTK::FACPreconditionerStrategy> fac_operator =
			new StaggeredStokesIBLevelRelaxationFACOperator(object_name + "::StaggeredStokesLevelRelaxationFACOperator",
															input_db,
															default_options_prefix);
		return new StaggeredStokesFACPreconditioner(object_name, fac_operator, input_db, default_options_prefix);
	} // allocate_solver

    /*!
     * \brief Set the PoissonSpecifications object used to specify the
     * coefficients for the momentum equation in the incompressible Stokes
     * operator.
     */
    virtual void setVelocityPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& U_problem_coefs);

	/*!
	 * \brief Set if velocity and pressure have nullspace.
	 */
	virtual void setComponentsHaveNullspace(const bool has_velocity_nullspace,
											const bool has_pressure_nullspace);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a U_bc_coefs may be NULL.  In this case,
     * homogeneous Dirichlet boundary conditions are employed for that data
     * depth.  \a P_bc_coef may also be NULL; in that case, homogeneous Neumann
     * boundary conditions are employed for the pressure.
     *
     * \param U_bc_coefs  IBTK::Vector of pointers to objects that can set the Robin boundary
     *condition coefficients for the velocity
     * \param P_bc_coef   Pointer to object that can set the Robin boundary condition
     *coefficients
     *for the pressure
     */
    virtual void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                                    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef);

    /*!
     * \brief Set the StokesSpecifications object and timestep size used to specify
     * the coefficients for the time-dependent incompressible Stokes operator.
     */
    virtual void setPhysicalBoundaryHelper(SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper);

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
    void setResetLevels(int coarsest_ln, int finest_ln);

    /*!
     * \brief Specify the smoother type.
     */
    void setSmootherType(const std::string& smoother_type);

    /*!
     * \brief Specify the coarse level solver.
     */
    void setCoarseSolverType(const std::string& coarse_solver_type);

    /*!
     * \brief Set the maximum number of iterations for the coarse level solve.
     *
     * If the coarse level solver uses a maximum number of iterations parameter,
     * the specified value is used.  If the coarse level solver does not use
     * such a stopping parameter, implementations are free to ignore this value.
     */
    void setCoarseSolverMaxIterations(int coarse_solver_max_iterations);

    /*!
     * \brief Set the absolute residual tolerance for convergence for coarse
     * level solve.
     *
     * If the coarse level solver uses a absolute convergence tolerance
     * parameter, the specified value is used.  If the coarse level solver does
     * not use such a stopping parameter, implementations are free to ignore
     * this value.
     */
    void setCoarseSolverAbsoluteTolerance(double coarse_solver_abs_residual_tol);

    /*!
     * \brief Set the relative residual tolerance for convergence for coarse
     * level solve.
     *
     * If the coarse level solver uses a relative convergence tolerance
     * parameter, the specified value is used.  If the coarse level solver does
     * not use such a stopping parameter, implementations are free to ignore
     * this value.
     */
    void setCoarseSolverRelativeTolerance(double coarse_solver_rel_residual_tol);

    /*!
     * \brief Set the prolongation methods.
     */
    void setProlongationMethods(const std::string& U_prolongation_method, const std::string& P_prolongation_method);

    /*!
     * \brief Set the restriction methods.
     */
    void setRestrictionMethods(const std::string& U_restriction_method, const std::string& P_restriction_method);

	/*!
	 * \brief Set the IB-force Jacobian at the finest patch level (where the
	 * structure resides).
	 */
	void setIBForceJacobian(Mat& A);

	/*!
	 * \brief Set the IB-interpolation operator at the finest patch level (where the
	 * structure resides).
	 */
	void setIBInterpOp(Mat& J);

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
    void restrictResidual(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& src,
                          SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dst,
                          int dst_ln);

    /*!
     * \brief Prolong the error quantity to the specified level from the next
     * coarser level.
     *
     * \param src source error vector
     * \param dst destination error vector
     * \param dst_ln destination level number of data transfer
     */
    void prolongError(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& src,
                      SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dst,
                      int dst_ln);

    /*!
     * \brief Prolong the error quantity to the specified level from the next
     * coarser level and apply the correction to the fine-level error.
     *
     * \param src source error vector
     * \param dst destination error vector
     * \param dst_ln destination level number of data transfer
     */
    void prolongErrorAndCorrect(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& src,
                                SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& dst,
                                int dst_ln);

    /*!
     * \brief Solve the residual equation Ae=r on the coarsest level of the
     * patch hierarchy.
     *
     * \param error error vector
     * \param residual residual vector
     * \param coarsest_ln coarsest level number
     */
    bool solveCoarsestLevel(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
                            const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                            int coarsest_ln);

    /*!
     * \brief Compute the composite-grid residual on the specified range of
     * levels of the patch hierarchy.
     */
    void computeResidual(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                         const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
                         int coarsest_level_num,
                         int finest_level_num);

	/*!
	 * \brief Perform a given number of relaxations on the error.
	 *
	 * \param error error vector
	 * \param residual residual vector
	 * \param level_num level number
	 * \param num_sweeps number of sweeps to perform
	 * \param performing_pre_sweeps boolean value that is true when pre-smoothing sweeps are
	 *being
	 *performed
	 * \param performing_post_sweeps boolean value that is true when post-smoothing sweeps are
	 *being
	 *performed
	 */
	void smoothError(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& error,
					 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& residual,
					 int level_num,
					 int num_sweeps,
					 bool performing_pre_sweeps,
					 bool performing_post_sweeps);

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
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs);

    /*!
     * \brief Remove all hierarchy-dependent data.
     *
     * Remove all hierarchy-dependent data set by initializeOperatorState().
     *
     * \see initializeOperatorState
     */
    void deallocateOperatorState();

    /*!
     * \brief Allocate scratch data.
     */
    void allocateScratchData();

    /*!
     * \brief Deallocate scratch data.
     */
    void deallocateScratchData();

    //\}

private:

	/*!
	 * \brief Default constructor.
	 *
	 * \note This constructor is not implemented and should not be used.
	 */
	StaggeredStokesIBLevelRelaxationFACOperator();

	/*!
	 * \brief Copy constructor.
	 *
	 * \note This constructor is not implemented and should not be used.
	 *
	 * \param from The value to copy to this object.
	 */
	StaggeredStokesIBLevelRelaxationFACOperator(const StaggeredStokesIBLevelRelaxationFACOperator& from);

	/*!
	 * \brief Assignment operator.
	 *
	 * \note This operator is not implemented and should not be used.
	 *
	 * \param that The value to assign to this object.
	 *
	 * \return A reference to this object.
	 */
	StaggeredStokesIBLevelRelaxationFACOperator& operator=(const StaggeredStokesIBLevelRelaxationFACOperator& that);

    /*!
     * \name Methods for executing, caching, and resetting communication
     * schedules.
     */
    //\{

    /*!
     * \brief Execute a refinement schedule for prolonging data.
     */
    void xeqScheduleProlongation(const std::pair<int, int>& dst_idxs, const std::pair<int, int>& src_idxs, int dst_ln);

    /*!
     * \brief Execute schedule for restricting solution or residual to the
     * specified level.
     */
    void xeqScheduleRestriction(const std::pair<int, int>& dst_idxs, const std::pair<int, int>& src_idxs, int dst_ln);

    /*!
     * \brief Execute schedule for filling ghosts on the specified level.
     */
    void xeqScheduleGhostFillNoCoarse(const std::pair<int, int>& dst_idxs, int dst_ln);

    /*!
     * \brief Execute schedule for synchronizing data on the specified level.
     */
    void xeqScheduleDataSynch(int dst_idx, int dst_ln);

    //\}

    /*
     * Problem specification.
     */
    SAMRAI::solv::PoissonSpecifications d_U_problem_coefs;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_default_U_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_U_bc_coefs;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_default_P_bc_coef;
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* d_P_bc_coef;

    /*
     * Boundary condition helper object.
     */
    SAMRAI::tbox::Pointer<StaggeredStokesPhysicalBoundaryHelper> d_bc_helper;

    /*!
     * \name Hierarchy-dependent objects.
     */
    //\{

    /*
     * Solution and rhs vectors.
     */
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM, double> > d_solution, d_rhs;

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
    std::vector<SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> > d_level_bdry_fill_ops;
    std::vector<SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> > d_level_math_ops;

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
    std::string d_smoother_type;

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
    std::string d_coarse_solver_type, d_coarse_solver_default_options_prefix;
    double d_coarse_solver_rel_residual_tol;
    double d_coarse_solver_abs_residual_tol;
    int d_coarse_solver_max_iterations;
	SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesPETScLevelSolver> d_coarse_solver;
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_coarse_solver_db;

	/*
	 * Level solvers and solver parameters.
	 */
	std::string d_level_solver_type, d_level_solver_default_options_prefix;
	double d_level_solver_abs_residual_tol, d_level_solver_rel_residual_tol;
	int d_level_solver_max_iterations;
	std::vector<SAMRAI::tbox::Pointer<IBAMR::StaggeredStokesPETScLevelSolver> >
		d_level_solvers;
	SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> d_level_solver_db;

	// Nullspace info.
	bool d_has_velocity_nullspace, d_has_pressure_nullspace;

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
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_U_op_stencil_fill_pattern,
        d_P_op_stencil_fill_pattern, d_U_synch_fill_pattern;

    //\}

    /*
     * Combined U & P physical boundary operator.
     */
    SAMRAI::xfer::RefinePatchStrategy<NDIM>* d_U_P_bc_op;

    /*
     * Error prolongation (refinement) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > d_U_prolongation_refine_operator,
        d_P_prolongation_refine_operator;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefinePatchStrategy<NDIM> > d_prolongation_refine_patch_strategy;
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_prolongation_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_prolongation_refine_schedules;

    /*
     * Residual restriction (coarsening) operator.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > d_U_restriction_coarsen_operator,
        d_P_restriction_coarsen_operator;
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

    /*
	 * Application ordering of u from MAC DOFs on various patch levels.
	 */
	std::vector<AO> d_u_app_ordering;

	/*
	 * Eulerian data for storing u and p DOFs indexing.
	 */
	std::vector<std::vector<int> > d_num_dofs_per_proc;
	int d_u_dof_index_idx, d_p_dof_index_idx;
	SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM, int> > d_u_dof_index_var;
	SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM, int> > d_p_dof_index_var;

	/*
	 * Jacobian of the elasticity force at the finest patch level.
	 */
	Mat d_mat_A;

	/*
	 * IB interpolation operator J for the finest patch level.
	 */
	Mat d_mat_J;

	/*
	 * Data structures for elasticity and prolongation operator respresentation 
	 * on various patch levels.
	 */
	std::vector<Mat> d_mat_SAJ, d_mat_prolongation;
	std::vector<Vec> d_mat_scale_restriction;

	/*
	 * Mappings from patch indices to patch operators.
	 */
	std::vector<std::vector<boost::array<SAMRAI::hier::BoxList<NDIM>, NDIM> > > d_patch_side_bc_box_overlap;
	std::vector<std::vector<SAMRAI::hier::BoxList<NDIM> > > d_patch_cell_bc_box_overlap;

};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_StaggeredStokesIBLevelRelaxationFACOperator
