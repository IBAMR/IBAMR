// Filename: INSStaggeredProjectionPreconditioner.C
// Last modified: <17.Aug.2009 16:21:09 griffith@boyce-griffiths-mac-pro.local>
// Created on 29 Apr 2008 by Boyce Griffith (griffith@box230.cims.nyu.edu)

#include "INSStaggeredProjectionPreconditioner.h"

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

INSStaggeredProjectionPreconditioner::INSStaggeredProjectionPreconditioner(
    const INSCoefs& problem_coefs,
    SAMRAI::solv::RobinBcCoefStrategy<NDIM>* Phi_bc_coef,
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
      d_pressure_helmholtz_spec("INSStaggeredProjectionPreconditioner::pressure_helmholtz_spec"),
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
      d_Phi_bc_coef(Phi_bc_coef),
      d_Phi_bdry_fill_op(SAMRAI::tbox::Pointer<IBTK::HierarchyGhostCellInterpolation>(NULL)),
      d_no_fill_op(NULL),
      d_hierarchy(NULL),
      d_coarsest_ln(-1),
      d_finest_ln(-1),
      d_Phi_var(NULL),
      d_F_var(NULL),
      d_Phi_scratch_idx(-1),
      d_F_scratch_idx(-1)
{
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext> context = var_db->getContext("INSStaggeredProjectionPreconditionerOperator::CONTEXT");

    const std::string Phi_var_name = "INSStaggeredProjectionPreconditioner::Phi";
    d_Phi_var = var_db->getVariable(Phi_var_name);
    if (d_Phi_var.isNull())
    {
        d_Phi_var = new SAMRAI::pdat::CellVariable<NDIM,double>(Phi_var_name);
        d_Phi_scratch_idx = var_db->registerVariableAndContext(d_Phi_var, context, SAMRAI::hier::IntVector<NDIM>(CELLG));
    }
    else
    {
        d_Phi_scratch_idx = var_db->mapVariableAndContextToIndex(d_Phi_var, context);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Phi_scratch_idx >= 0);
#endif
    const std::string F_var_name = "INSStaggeredProjectionPreconditioner::F";
    d_F_var = var_db->getVariable(F_var_name);
    if (d_F_var.isNull())
    {
        d_F_var = new SAMRAI::pdat::CellVariable<NDIM,double>(F_var_name);
        d_F_scratch_idx = var_db->registerVariableAndContext(d_F_var, context, SAMRAI::hier::IntVector<NDIM>(CELLG));
    }
    else
    {
        d_F_scratch_idx = var_db->mapVariableAndContextToIndex(d_F_var, context);
    }
#ifdef DEBFG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_F_scratch_idx >= 0);
#endif
    return;
}// INSStaggeredProjectionPreconditioner

INSStaggeredProjectionPreconditioner::~INSStaggeredProjectionPreconditioner()
{
    deallocateSolverState();
    return;
}// ~INSStaggeredProjectionPreconditioner

void
INSStaggeredProjectionPreconditioner::setTimeInterval(
    const double current_time,
    const double new_time,
    const double dt)
{
    const double rho    = d_problem_coefs.getRho();
    const double mu     = d_problem_coefs.getMu();
    const double lambda = d_problem_coefs.getLambda();
    d_current_time = current_time;
    d_new_time = new_time;
    d_dt = dt;
    d_pressure_helmholtz_spec.setCConstant(1.0+0.5*d_dt*lambda/rho);
    d_pressure_helmholtz_spec.setDConstant(   -0.5*d_dt*mu    /rho);
    return;
}// setTimeInterval

bool
INSStaggeredProjectionPreconditioner::solveSystem(
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& x,
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double>& b)
{
    // Initialize the solver (if necessary).
    const bool deallocate_at_completion = !d_is_initialized;
    if (!d_is_initialized) initializeSolverState(x,b);

    // Problem coefficients.
    const double rho    = d_problem_coefs.getRho();
//  const double mu     = d_problem_coefs.getMu();
//  const double lambda = d_problem_coefs.getLambda();

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
    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_in_vec;
    U_in_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>("INSStaggeredProjectionPreconditioner::U_in", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_in_vec->addComponent(U_in_sc_var, U_in_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > U_out_vec;
    U_out_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>("INSStaggeredProjectionPreconditioner::U_out", d_hierarchy, d_coarsest_ln, d_finest_ln);
    U_out_vec->addComponent(U_out_sc_var, U_out_idx, d_wgt_sc_idx, d_hier_sc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > Phi_scratch_vec;
    Phi_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>("INSStaggeredProjectionPreconditioner::Phi_scratch", d_hierarchy, d_coarsest_ln, d_finest_ln);
    Phi_scratch_vec->addComponent(d_Phi_var, d_Phi_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::solv::SAMRAIVectorReal<NDIM,double> > F_scratch_vec;
    F_scratch_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>("INSStaggeredProjectionPreconditioner::F_scratch", d_hierarchy, d_coarsest_ln, d_finest_ln);
    F_scratch_vec->addComponent(d_F_var, d_F_scratch_idx, d_wgt_cc_idx, d_hier_cc_data_ops);

    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = x.getPatchHierarchy();
    TBOX_ASSERT(hierarchy == b.getPatchHierarchy());

    // Solve for u^{*}.
    d_velocity_helmholtz_solver->solveSystem(*U_out_vec,*U_in_vec);

    // Compute F = -(rho/dt)*(P_in + div u^{*}).
    const bool u_star_cf_bdry_synch = true;
    d_hier_math_ops->div(
        d_F_scratch_idx, d_F_var, // dst
        -rho/d_dt,                // alpha
        U_out_idx, U_out_sc_var,  // src1
        d_no_fill_op,             // src1_bdry_fill
        d_new_time,               // src1_bdry_fill_time
        u_star_cf_bdry_synch,     // src1_cf_bdry_synch
        -rho/d_dt,                // beta
        P_in_idx, P_in_cc_var);   // src2

    // Solve -div grad Phi = F = -(rho/dt)*(P_in + div u^{*}).
    d_pressure_poisson_solver->solveSystem(*Phi_scratch_vec,*F_scratch_vec);

    // Use Phi to project u^{*}.
    const bool u_new_cf_bdry_synch = true;
    d_hier_math_ops->grad(
        U_out_idx, U_out_sc_var,         // dst
        u_new_cf_bdry_synch,             // dst_cf_bdry_synch
        -d_dt/rho,                       // alpha
        d_Phi_scratch_idx, d_Phi_var,    // src1
        d_Phi_bdry_fill_op,              // src1_bdry_fill
        0.5*(d_current_time+d_new_time), // src1_bdry_fill_time
        1.0,                             // beta
        U_out_idx, U_out_sc_var);        // src2

    // Compute P_out.
    d_hier_math_ops->laplace(
        P_out_idx, P_out_cc_var,          // dst
        d_pressure_helmholtz_spec,        // Poisson spec
        d_Phi_scratch_idx, d_Phi_var,     // src
        d_no_fill_op,                     // src_bdry_fill
        0.5*(d_current_time+d_new_time)); // src_bdry_fill_time
    if (d_normalize_pressure)
    {
        const double P_mean = (1.0/d_volume)*d_hier_cc_data_ops->integral(P_out_idx, d_wgt_cc_idx);
        d_hier_cc_data_ops->addScalar(P_out_idx, P_out_idx, -P_mean);
    }

    // Deallocate the solver (if necessary).
    if (deallocate_at_completion) deallocateSolverState();
    return true;
}// solveSystem

void
INSStaggeredProjectionPreconditioner::initializeSolverState(
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
    InterpolationTransactionComponent Phi_scratch_component(d_Phi_scratch_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_Phi_bc_coef);
    d_Phi_bdry_fill_op = new IBTK::HierarchyGhostCellInterpolation();
    d_Phi_bdry_fill_op->initializeOperatorState(Phi_scratch_component, d_hierarchy);

    // Allocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_Phi_scratch_idx))
        {
            level->allocatePatchData(d_Phi_scratch_idx);
        }
        if (!level->checkAllocated(d_F_scratch_idx))
        {
            level->allocatePatchData(d_F_scratch_idx);
        }
    }
    d_is_initialized = true;
    return;
}// initializeSolverState

void
INSStaggeredProjectionPreconditioner::deallocateSolverState()
{
    if (!d_is_initialized) return;

    // Deallocate scratch data.
    for (int ln = d_coarsest_ln; ln <= d_finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_Phi_scratch_idx))
        {
            level->deallocatePatchData(d_Phi_scratch_idx);
        }
        if (level->checkAllocated(d_F_scratch_idx))
        {
            level->deallocatePatchData(d_F_scratch_idx);
        }
    }
    d_is_initialized = false;
    return;
}// deallocateSolverState

void
INSStaggeredProjectionPreconditioner::setInitialGuessNonzero(
    bool initial_guess_nonzero)
{
    // intentionally blank
    return;
}// setInitialGuessNonzero

bool
INSStaggeredProjectionPreconditioner::getInitialGuessNonzero() const
{
    // intentionally blank
    return true;
}// getInitialGuessNonzero

void
INSStaggeredProjectionPreconditioner::setMaxIterations(
    int max_iterations)
{
    // intentionally blank
    return;
}// setMaxIterations

int
INSStaggeredProjectionPreconditioner::getMaxIterations() const
{
    // intentionally blank
    return 1;
}// getMaxIterations

void
INSStaggeredProjectionPreconditioner::setAbsoluteTolerance(
    double abs_residual_tol)
{
    // intentionally blank
    return;
}// setAbsoluteTolerance

double
INSStaggeredProjectionPreconditioner::getAbsoluteTolerance() const
{
    // intentionally blank
    return 0.0;
}// getAbsoluteTolerance

void
INSStaggeredProjectionPreconditioner::setRelativeTolerance(
    double rel_residual_tol)
{
    // intentionally blank
    return;
}// setRelativeTolerance

double
INSStaggeredProjectionPreconditioner::getRelativeTolerance() const
{
    // intentionally blank
    return 0.0;
}// getRelativeTolerance

int
INSStaggeredProjectionPreconditioner::getNumIterations() const
{
    // intentionally blank
    return 0;
}// getNumIterations

double
INSStaggeredProjectionPreconditioner::getResidualNorm() const
{
    return 0.0;
}// getResidualNorm

void
INSStaggeredProjectionPreconditioner::enableLogging(
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
template class SAMRAI::tbox::Pointer<IBAMR::INSStaggeredProjectionPreconditioner>;

//////////////////////////////////////////////////////////////////////////////
