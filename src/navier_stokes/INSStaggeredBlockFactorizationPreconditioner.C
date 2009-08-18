// Filename: INSStaggeredBlockFactorizationPreconditioner.C
// Last modified: <17.Aug.2009 16:20:45 griffith@boyce-griffiths-mac-pro.local>
// Created on 22 Sep 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredBlockFactorizationPreconditioner.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// C++ STDLIB INCLUDES
#include <limits>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);
static const int SIDEG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSStaggeredBlockFactorizationPreconditioner::INSStaggeredBlockFactorizationPreconditioner(
    const INSCoefs& problem_coefs,
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* P_bc_coef,
    const bool normalize_pressure,
    SAMRAI::tbox::Pointer<IBTK::LinearSolver> velocity_helmholtz_solver,
    SAMRAI::tbox::Pointer<IBTK::LinearSolver> pressure_poisson_solver,
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchyCellDataOpsReal<NDIM,double> > hier_cc_data_ops,
    SAMRAI::tbox::Pointer<SAMRAI::math::HierarchySideDataOpsReal<NDIM,double> > hier_sc_data_ops,
    SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops)
    : d_do_log(false),
      d_is_initialized(false),
      d_current_time(std::numeric_limits<double>::quiet_NaN()),
      d_new_time(std::numeric_limits<double>::quiet_NaN()),
      d_dt(std::numeric_limits<double>::quiet_NaN()),
      d_problem_coefs(problem_coefs),
      d_pressure_helmholtz_spec("INSStaggeredBlockFactorizationPreconditioner::pressure_helmholtz_spec"),
      d_normalize_pressure(normalize_pressure),
      d_velocity_helmholtz_solver(velocity_helmholtz_solver),
      d_pressure_poisson_solver(pressure_poisson_solver),
      d_hier_cc_data_ops(hier_cc_data_ops),
      d_hier_sc_data_ops(hier_sc_data_ops),
      d_hier_math_ops(hier_math_ops),
      d_wgt_cc_var(NULL),
      d_wgt_sc_var(NULL),
      d_wgt_cc_idx(-1),
      d_wgt_sc_idx(-1),
      d_volume(std::numeric_limits<double>::quiet_NaN()),
      d_P_bc_coef(P_bc_coef),
      d_P_bdry_fill_op(SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>(NULL)),
      d_no_fill_op(NULL),
      d_hierarchy(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_U_var(NULL),
      d_U_scratch_idx(-1),
      d_P_var(NULL),
      d_P_scratch_idx(-1)
{
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> context = var_db->getContext("INSStaggeredBlockFactorizationPreconditionerOperator::CONTEXT");

    const std::string U_var_name = "INSStaggeredBlockFactorizationPreconditioner::U";
    d_U_var = var_db->getVariable(U_var_name);
    if (d_U_var.isNull())
    {
        d_U_var = new SAMRAI::pdat::SideVariable<NDIM,double>(U_var_name);
        d_U_scratch_idx = var_db->registerVariableAndContext(d_U_var, context, SAMRAI::hier::IntVector<NDIM>(SIDEG));
    }
    else
    {
        d_U_scratch_idx = var_db->mapVariableAndContextToIndex(d_U_var, context);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_U_scratch_idx >= 0);
#endif
    const std::string P_var_name = "INSStaggeredBlockFactorizationPreconditioner::P";
    d_P_var = var_db->getVariable(P_var_name);
    if (d_P_var.isNull())
    {
        d_P_var = new SAMRAI::pdat::CellVariable<NDIM,double>(P_var_name);
        d_P_scratch_idx = var_db->registerVariableAndContext(d_P_var, context, SAMRAI::hier::IntVector<NDIM>(CELLG));
    }
    else
    {
        d_P_scratch_idx = var_db->mapVariableAndContextToIndex(d_P_var, context);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_P_scratch_idx >= 0);
#endif
    return;
}// INSStaggeredBlockFactorizationPreconditioner

INSStaggeredBlockFactorizationPreconditioner::~INSStaggeredBlockFactorizationPreconditioner()
{
    deallocateSolverState();
    return;
}// ~INSStaggeredBlockFactorizationPreconditioner

void
INSStaggeredBlockFactorizationPreconditioner::setTimeInterval(
    const double current_time,
    const double new_time)
{
    const double rho    = d_problem_coefs.getRho();
    const double mu     = d_problem_coefs.getMu();
    const double lambda = d_problem_coefs.getLambda();
    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = d_new_time-d_current_time;
    d_pressure_helmholtz_spec.setCConstant(-(rho/d_dt+0.5*lambda));
    d_pressure_helmholtz_spec.setDConstant(-(        -0.5*mu    ));
    return;
}// setTimeInterval

bool
INSStaggeredBlockFactorizationPreconditioner::solveSystem(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b)
{
    // Initialize the solver (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeSolverState(x,b);

    // Get the vector components.
    const int U_in_idx = b.getComponentDescriptorIndex(0);
    const int P_in_idx = b.getComponentDescriptorIndex(1);

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_in_var = b.getComponentVariable(0);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_in_var = b.getComponentVariable(1);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_in_sc_var = U_in_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_in_cc_var = P_in_var;

    const int U_out_idx = x.getComponentDescriptorIndex(0);
    const int P_out_idx = x.getComponentDescriptorIndex(1);

    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& U_out_var = x.getComponentVariable(0);
    const SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> >& P_out_var = x.getComponentVariable(1);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > U_out_sc_var = U_out_var;
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > P_out_cc_var = P_out_var;

    // Setup the component solver vectors.
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_scratch_vec;
    U_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredBlockFactorizationPreconditioner::U_scratch", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_scratch_vec->addComponent(d_U_var, d_U_scratch_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > P_scratch_vec;
    P_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredBlockFactorizationPreconditioner::P_scratch", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_scratch_vec->addComponent(d_P_var, d_P_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_out_vec;
    U_out_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredBlockFactorizationPreconditioner::U_out", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_out_vec->addComponent(U_out_sc_var, U_out_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > P_in_vec;
    P_in_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredBlockFactorizationPreconditioner::P_in", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_in_vec->addComponent(P_in_cc_var, P_in_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > P_out_vec;
    P_out_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        "INSStaggeredBlockFactorizationPreconditioner::P_out", d_hierarchy, d_coarsest_ln, d_finest_ln);
    P_out_vec->addComponent(P_out_cc_var, P_out_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    // Setup the interpolation transaction information.
    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent P_out_transaction_comp(P_out_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef);
    InterpolationTransactionComponent P_scratch_transaction_comp(d_P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef);

    // Solve the pressure sub-problem.
    //
    // P_out := -((rho/dt)*I-0.5*mu*L) * (-L)^{-1} * P_in
    d_pressure_poisson_solver->solveSystem(*P_scratch_vec,*P_in_vec);
    d_hier_math_ops->laplace(
        P_out_idx, P_out_cc_var,   // dst
        d_pressure_helmholtz_spec, // Poisson spec
        d_P_scratch_idx, d_P_var,  // src
        d_P_bdry_fill_op,          // src_bdry_fill
        d_current_time+0.5*d_dt);  // src_bdry_fill_time
    if (d_normalize_pressure)
    {
        const double P_mean = (1.0/d_volume)*d_hier_cc_data_ops->integral(P_out_idx, d_wgt_cc_idx);
        d_hier_cc_data_ops->addScalar(P_out_idx, P_out_idx, -P_mean);
    }

    // Compute the right-hand-side for the Helmholtz solve.
    d_P_bdry_fill_op->resetTransactionComponent(P_out_transaction_comp);
    static const bool cf_bdry_synch = true;
    d_hier_math_ops->grad(
        d_U_scratch_idx, d_U_var, // dst
        cf_bdry_synch,            // dst_cf_bdry_synch
        -1.0,                     // alpha
        P_out_idx, P_out_cc_var,  // src1
        d_P_bdry_fill_op,         // src1_bdry_fill
        d_current_time+0.5*d_dt,  // src1_bdry_fill_time
        1.0,                      // beta
        U_in_idx, U_in_sc_var);   // src2
    d_P_bdry_fill_op->resetTransactionComponent(P_scratch_transaction_comp);

    // Solve the velocity sub-problem.
    //
    // U_out := (rho/dt)*I-0.5*mu*L)^{-1} * [U_in + Grad * ((rho/dt)*I-0.5*mu*L) * (-L)^{-1} * P_in]
    //        = (rho/dt)*I-0.5*mu*L)^{-1} * [U_in - Grad * P_out]
    d_velocity_helmholtz_solver->solveSystem(*U_out_vec,*U_scratch_vec);

    // Deallocate the solver (if necessary).
    if (deallocate_at_completion) deallocateSolverState();
    return true;
}// solveSystem

void
INSStaggeredBlockFactorizationPreconditioner::initializeSolverState(
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    const SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b)
{
    if (d_is_initialized) deallocateSolverState();

    // Get the hierarchy configuration.
    d_hierarchy = x.getPatchHierarchy();
    d_coarsest_ln = x.getCoarsestLevelNumber();
    d_finest_ln = x.getFinestLevelNumber();
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_hierarchy == b.getPatchHierarchy());
    TBOX_ASSERT(d_coarsest_ln == b.getCoarsestLevelNumber());
    TBOX_ASSERT(d_finest_ln == b.getFinestLevelNumber());
#endif
    d_wgt_cc_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_sc_var = d_hier_math_ops->getSideWeightVariable();
    d_wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    d_wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    typedef IBTK::HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
    InterpolationTransactionComponent P_scratch_component(d_P_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_P_bc_coef);
    d_P_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_P_bdry_fill_op->initializeOperatorState(P_scratch_component, d_hierarchy);

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_U_scratch_idx))
        {
            level->allocatePatchData(d_U_scratch_idx);
        }
        if (!level->checkAllocated(d_P_scratch_idx))
        {
            level->allocatePatchData(d_P_scratch_idx);
        }
    }
    d_is_initialized = true;
    return;
}// initializeSolverState

void
INSStaggeredBlockFactorizationPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_U_scratch_idx))
        {
            level->deallocatePatchData(d_U_scratch_idx);
        }
        if (level->checkAllocated(d_P_scratch_idx))
        {
            level->deallocatePatchData(d_P_scratch_idx);
        }
    }
    d_is_initialized = false;
    return;
}// deallocateSolverState

void
INSStaggeredBlockFactorizationPreconditioner::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    // intentionally blank
    return;
}// setInitialGuessNonzero

bool
INSStaggeredBlockFactorizationPreconditioner::getInitialGuessNonzero() const
{
    // intentionally blank
    return true;
}// getInitialGuessNonzero

void
INSStaggeredBlockFactorizationPreconditioner::setMaxIterations(
    int max_iterations)
{
    // intentionally blank
    return;
}// setMaxIterations

int
INSStaggeredBlockFactorizationPreconditioner::getMaxIterations() const
{
    // intentionally blank
    return 1;
}// getMaxIterations

void
INSStaggeredBlockFactorizationPreconditioner::setAbsoluteTolerance(
    double abs_residual_tol)
{
    // intentionally blank
    return;
}// setAbsoluteTolerance

double
INSStaggeredBlockFactorizationPreconditioner::getAbsoluteTolerance() const
{
    // intentionally blank
    return 0.0;
}// getAbsoluteTolerance

void
INSStaggeredBlockFactorizationPreconditioner::setRelativeTolerance(
    double rel_residual_tol)
{
    // intentionally blank
    return;
}// setRelativeTolerance

double
INSStaggeredBlockFactorizationPreconditioner::getRelativeTolerance() const
{
    // intentionally blank
    return 0.0;
}// getRelativeTolerance

int
INSStaggeredBlockFactorizationPreconditioner::getNumIterations() const
{
    // intentionally blank
    return 0;
}// getNumIterations

double
INSStaggeredBlockFactorizationPreconditioner::getResidualNorm() const
{
    return 0.0;
}// getResidualNorm

void
INSStaggeredBlockFactorizationPreconditioner::enableLogging(
    bool enabled)
{
    d_do_log = enabled;
    return;
}// enableLogging

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredBlockFactorizationPreconditioner>;

//////////////////////////////////////////////////////////////////////////////
