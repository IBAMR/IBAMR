// Filename: PoissonFACPreconditionerStrategy.h
// Created on 27 Apr 2012 by Boyce Griffith
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

#ifndef included_PoissonFACPreconditionerStrategy
#define included_PoissonFACPreconditionerStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefinePatchStrategy.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/xfer/VariableFillPattern.h"
#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"


namespace IBTK
{
class HierarchyGhostCellInterpolation;
class HierarchyMathOps;
} // namespace IBTK
namespace SAMRAI
{
namespace hier
{

class Variable;
} // namespace hier
namespace math
{
template < class TYPE>
class HierarchyDataOpsReal;
} // namespace math
namespace solv
{

class RobinBcCoefStrategy;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
namespace xfer
{

class CoarsenSchedule;

class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBTK
{
/*!
 * \brief Class PoissonFACPreconditionerStrategy is an abstract
 * FACPreconditionerStrategy implementing many of the operations required by
 * smoothers for the Poisson equation and related problems.
 *
 * Sample parameters for initialization from database (and their default
 * values): \verbatim

 smoother_type = "DEFAULT"                    // see setSmootherType()
 prolongation_method = "DEFAULT"              // see setProlongationMethod()
 restriction_method = "DEFAULT"               // see setRestrictionMethod()
 coarse_solver_type = "DEFAULT"               // see setCoarseSolverType()
 coarse_solver_rel_residual_tol = 1.0e-5      // see setCoarseSolverRelativeTolerance()
 coarse_solver_abs_residual_tol = 1.0e-50     // see setCoarseSolverAbsoluteTolerance()
 coarse_solver_max_iterations = 10            // see setCoarseSolverMaxIterations()
 \endverbatim
*/
class PoissonFACPreconditionerStrategy : public FACPreconditionerStrategy
{
public:
    /*!
     * \brief Constructor.
     */
    PoissonFACPreconditionerStrategy(const std::string& object_name,
                                     boost::shared_ptr<SAMRAI::hier::Variable > scratch_var,
                                     int ghost_cell_width,
                                     boost::shared_ptr<SAMRAI::tbox::Database> input_db,
                                     const std::string& default_options_prefix);

    /*!
     * \brief Destructor.
     */
    ~PoissonFACPreconditionerStrategy();

    /*!
     * \brief Set the SAMRAI::solv::PoissonSpecifications object used to specify
     * the coefficients for the scalar-valued or vector-valued Laplace operator.
     */
    virtual void setPoissonSpecifications(const SAMRAI::solv::PoissonSpecifications& poisson_spec);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy object used to specify
     * physical boundary conditions.
     *
     * \note \a bc_coef may be NULL.  In this case, default boundary conditions
     * (as supplied to the class constructor) are employed.
     *
     * \param bc_coef  boost::shared_ptr to an object that can set the Robin boundary condition
     *coefficients
     */
    virtual void setPhysicalBcCoef(SAMRAI::solv::RobinBcCoefStrategy* bc_coef);

    /*!
     * \brief Set the SAMRAI::solv::RobinBcCoefStrategy objects used to specify
     * physical boundary conditions.
     *
     * \note Any of the elements of \a bc_coefs may be NULL.  In this case,
     * default boundary conditions (as supplied to the class constructor) are
     * employed for that data depth.
     *
     * \param bc_coefs  Vector of pointers to objects that can set the Robin boundary condition
     *coefficients
     */
    virtual void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy*>& bc_coefs);

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
    virtual void setSmootherType(const std::string& smoother_type) = 0;

    /*!
     * \brief Specify the coarse level solver.
     */
    virtual void setCoarseSolverType(const std::string& coarse_solver_type) = 0;

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
     * \brief Set the name of the prolongation method.
     */
    void setProlongationMethod(const std::string& prolongation_method);

    /*!
     * \brief Set the name of the restriction method.
     */
    void setRestrictionMethod(const std::string& restriction_method);

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
    void restrictResidual(const SAMRAI::solv::SAMRAIVectorReal<double>& src,
                          SAMRAI::solv::SAMRAIVectorReal<double>& dst,
                          int dst_ln);

    /*!
     * \brief Prolong the error quantity to the specified level from the next
     * coarser level.
     *
     * \param src source error vector
     * \param dst destination error vector
     * \param dst_ln destination level number of data transfer
     */
    void prolongError(const SAMRAI::solv::SAMRAIVectorReal<double>& src,
                      SAMRAI::solv::SAMRAIVectorReal<double>& dst,
                      int dst_ln);

    /*!
     * \brief Prolong the error quantity to the specified level from the next
     * coarser level and apply the correction to the fine-level error.
     *
     * \param src source error vector
     * \param dst destination error vector
     * \param dst_ln destination level number of data transfer
     */
    void prolongErrorAndCorrect(const SAMRAI::solv::SAMRAIVectorReal<double>& src,
                                SAMRAI::solv::SAMRAIVectorReal<double>& dst,
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
    void initializeOperatorState(const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
                                 const SAMRAI::solv::SAMRAIVectorReal<double>& rhs);

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

protected:
    /*!
     * \brief Compute implementation-specific hierarchy-dependent data.
     */
    virtual void initializeOperatorStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<double>& solution,
                                                    const SAMRAI::solv::SAMRAIVectorReal<double>& rhs,
                                                    int coarsest_reset_ln,
                                                    int finest_reset_ln) = 0;

    /*!
     * \brief Remove implementation-specific hierarchy-dependent data.
     */
    virtual void deallocateOperatorStateSpecialized(int coarsest_reset_ln, int finest_reset_ln) = 0;

    /*!
     * \name Methods for executing, caching, and resetting communication
     * schedules.
     */
    //\{

    /*!
     * \brief Execute a refinement schedule for prolonging data.
     */
    void xeqScheduleProlongation(int dst_idx, int src_idx, int dst_ln);

    /*!
     * \brief Execute schedule for restricting solution or residual to the
     * specified level.
     */
    void xeqScheduleRestriction(int dst_idx, int src_idx, int dst_ln);

    /*!
     * \brief Execute schedule for filling ghosts on the specified level.
     */
    void xeqScheduleGhostFillNoCoarse(int dst_idx, int dst_ln);

    /*!
     * \brief Execute schedule for synchronizing data on the specified level.
     */
    void xeqScheduleDataSynch(int dst_idx, int dst_ln);

    //\}

    /*
     * Problem specification.
     */
    SAMRAI::solv::PoissonSpecifications d_poisson_spec;
    boost::shared_ptr<SAMRAI::solv::RobinBcCoefStrategy> d_default_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy*> d_bc_coefs;

    /*
     * Ghost cell width.
     */
    const SAMRAI::hier::IntVector d_gcw;

    /*!
     * \name Hierarchy-dependent objects.
     */
    //\{

    /*
     * Solution and rhs vectors.
     */
    boost::shared_ptr<SAMRAI::solv::SAMRAIVectorReal<double> > d_solution, d_rhs;

    /*
     * Reference patch hierarchy and range of levels involved in the solve.
     *
     * This variable is non-null between the initializeOperatorState() and
     * deallocateOperatorState() calls.  It is not truly needed, because the
     * hierarchy is obtainable through variables in most function argument
     * lists.  We use it to enforce working on one hierarchy at a time.
     */
    boost::shared_ptr<SAMRAI::hier::PatchHierarchy > d_hierarchy;
    int d_coarsest_ln, d_finest_ln;

    /*
     * HierarchyDataOpsReal objects restricted to a single level of the patch
     * hierarchy.
     */
    std::vector<boost::shared_ptr<SAMRAI::math::HierarchyDataOpsReal<double> > > d_level_data_ops;

    /*
     * Level operators, used to compute composite-grid residuals.
     */
    std::vector<boost::shared_ptr<IBTK::HierarchyGhostCellInterpolation> > d_level_bdry_fill_ops;
    std::vector<boost::shared_ptr<IBTK::HierarchyMathOps> > d_level_math_ops;

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
    std::string d_prolongation_method;

    /*
     * The names of the coarsening operators used to restrict the fine grid
     * error or residual.
     */
    std::string d_restriction_method;

    /*
     * Coarse level solver parameters.
     */
    std::string d_coarse_solver_type, d_coarse_solver_default_options_prefix;
    double d_coarse_solver_rel_residual_tol;
    double d_coarse_solver_abs_residual_tol;
    int d_coarse_solver_max_iterations;

    //\}

    /*!
     * \name Internal context and scratch data.
     */
    //\{

    /*
     * Variable context for internally maintained hierarchy data.
     */
    boost::shared_ptr<SAMRAI::hier::VariableContext> d_context;

    /*
     * Patch descriptor index for scratch data.
     */
    int d_scratch_idx;

    //\}

    /*!
     * \name Various refine and coarsen objects.
     */
    //\{

    /*
     * Physical boundary operators.
     */
    boost::shared_ptr<RobinPhysBdryPatchStrategy> d_bc_op;

    /*
     * Coarse-fine interface interpolation objects.
     */
    boost::shared_ptr<CoarseFineBoundaryRefinePatchStrategy> d_cf_bdry_op;

    /*
     * Variable fill pattern object.
     */
    boost::shared_ptr<SAMRAI::xfer::VariableFillPattern > d_op_stencil_fill_pattern, d_synch_fill_pattern;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PoissonFACPreconditionerStrategy();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PoissonFACPreconditionerStrategy(const PoissonFACPreconditionerStrategy& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PoissonFACPreconditionerStrategy& operator=(const PoissonFACPreconditionerStrategy& that);

    /*!
     * \name Various refine and coarsen objects.
     */
    //\{

    /*
     * Error prolongation (refinement) operator.
     */
    boost::shared_ptr<SAMRAI::hier::RefineOperator > d_prolongation_refine_operator;
    boost::shared_ptr<SAMRAI::xfer::RefinePatchStrategy > d_prolongation_refine_patch_strategy;
    boost::shared_ptr<SAMRAI::xfer::RefineAlgorithm > d_prolongation_refine_algorithm;
    std::vector<boost::shared_ptr<SAMRAI::xfer::RefineSchedule > > d_prolongation_refine_schedules;

    /*
     * Residual restriction (coarsening) operator.
     */
    boost::shared_ptr<SAMRAI::hier::CoarsenOperator > d_restriction_coarsen_operator;
    boost::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm > d_restriction_coarsen_algorithm;
    std::vector<boost::shared_ptr<SAMRAI::xfer::CoarsenSchedule > > d_restriction_coarsen_schedules;

    /*
     * Refine operator for cell data from same level.
     */
    boost::shared_ptr<SAMRAI::xfer::RefineAlgorithm > d_ghostfill_nocoarse_refine_algorithm;
    std::vector<boost::shared_ptr<SAMRAI::xfer::RefineSchedule > > d_ghostfill_nocoarse_refine_schedules;

    /*
     * Operator for data synchronization on same level.
     */
    boost::shared_ptr<SAMRAI::xfer::RefineAlgorithm > d_synch_refine_algorithm;
    std::vector<boost::shared_ptr<SAMRAI::xfer::RefineSchedule > > d_synch_refine_schedules;

    //\}
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_PoissonFACPreconditionerStrategy
