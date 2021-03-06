// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2020 by the IBAMR developers
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

#ifndef included_IBTK_PoissonFACPreconditionerStrategy
#define included_IBTK_PoissonFACPreconditionerStrategy

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include "ibtk/CoarseFineBoundaryRefinePatchStrategy.h"
#include "ibtk/FACPreconditionerStrategy.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"

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
#include "tbox/Pointer.h"

#include <memory>
#include <string>
#include <vector>

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
class Variable;
} // namespace hier
namespace math
{
template <int DIM, class TYPE>
class HierarchyDataOpsReal;
} // namespace math
namespace solv
{
template <int DIM>
class RobinBcCoefStrategy;
} // namespace solv
namespace tbox
{
class Database;
} // namespace tbox
namespace xfer
{
template <int DIM>
class CoarsenSchedule;
template <int DIM>
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
    PoissonFACPreconditionerStrategy(std::string object_name,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > scratch_var,
                                     int ghost_cell_width,
                                     SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
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
     * \param bc_coef  Pointer to an object that can set the Robin boundary condition
     *coefficients
     */
    virtual void setPhysicalBcCoef(SAMRAI::solv::RobinBcCoefStrategy<NDIM>* bc_coef);

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
    virtual void setPhysicalBcCoefs(const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& bc_coefs);

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
     * \brief Zero the supplied vector.
     */
    void setToZero(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& vec, int level_num) override;

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
                          int dst_ln) override;

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
                      int dst_ln) override;

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
                                int dst_ln) override;

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
                                 const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs) override;

    /*!
     * \brief Remove all hierarchy-dependent data.
     *
     * Remove all hierarchy-dependent data set by initializeOperatorState().
     *
     * \see initializeOperatorState
     */
    void deallocateOperatorState() override;

    //\}

protected:
    /*!
     * \brief Compute implementation-specific hierarchy-dependent data.
     */
    virtual void initializeOperatorStateSpecialized(const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& solution,
                                                    const SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& rhs,
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
    std::unique_ptr<SAMRAI::solv::RobinBcCoefStrategy<NDIM> > d_default_bc_coef;
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> d_bc_coefs;

    /*
     * Ghost cell width.
     */
    SAMRAI::hier::IntVector<NDIM> d_gcw;

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
    int d_coarsest_ln = IBTK::invalid_level_number, d_finest_ln = IBTK::invalid_level_number;

    /*
     * HierarchyDataOpsReal objects restricted to a single level of the patch
     * hierarchy.
     */
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyDataOpsReal<NDIM, double> > > d_level_data_ops;

    /*
     * Level operators, used to compute composite-grid residuals.
     */
    std::vector<SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation> > d_level_bdry_fill_ops;
    std::vector<SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> > d_level_math_ops;

    /*
     * Range of levels to be reset the next time the operator is initialized.
     */
    bool d_in_initialize_operator_state = false;
    int d_coarsest_reset_ln = IBTK::invalid_level_number, d_finest_reset_ln = IBTK::invalid_level_number;

    //\}

    /*!
     * \name Solver configuration variables.
     */
    //\{

    /*
     * The kind of smoothing to perform.
     */
    std::string d_smoother_type = "DEFAULT";

    /*
     * The names of the refinement operators used to prolong the coarse grid
     * correction.
     */
    std::string d_prolongation_method = "DEFAULT";

    /*
     * The names of the coarsening operators used to restrict the fine grid
     * error or residual.
     */
    std::string d_restriction_method = "DEFAULT";

    /*
     * Coarse level solver parameters.
     */
    std::string d_coarse_solver_type = "DEFAULT", d_coarse_solver_default_options_prefix;
    double d_coarse_solver_rel_residual_tol = 1.0e-5;
    double d_coarse_solver_abs_residual_tol = 1.0e-50;
    int d_coarse_solver_max_iterations = 10;

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
     * Patch descriptor index for scratch data.
     */
    int d_scratch_idx = IBTK::invalid_index;

    //\}

    /*!
     * \name Various refine and coarsen objects.
     */
    //\{

    /*
     * Physical boundary operators.
     */
    SAMRAI::tbox::Pointer<RobinPhysBdryPatchStrategy> d_bc_op;

    /*
     * Coarse-fine interface interpolation objects.
     */
    SAMRAI::tbox::Pointer<CoarseFineBoundaryRefinePatchStrategy> d_cf_bdry_op;

    /*
     * Variable fill pattern object.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::VariableFillPattern<NDIM> > d_op_stencil_fill_pattern, d_synch_fill_pattern;

    //\}

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    PoissonFACPreconditionerStrategy() = delete;

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    PoissonFACPreconditionerStrategy(const PoissonFACPreconditionerStrategy& from) = delete;

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    PoissonFACPreconditionerStrategy& operator=(const PoissonFACPreconditionerStrategy& that) = delete;

    /*!
     * \name Various refine and coarsen objects.
     */
    //\{

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
     * Refine operator for cell data from same level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_ghostfill_nocoarse_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_ghostfill_nocoarse_refine_schedules;

    /*
     * Operator for data synchronization on same level.
     */
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > d_synch_refine_algorithm;
    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > d_synch_refine_schedules;

    //\}
};
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBTK_PoissonFACPreconditionerStrategy
