// Filename: HierarchyProjector.C
// Last modified: <25.Feb.2007 19:13:05 boyce@boyce-griffiths-powerbook-g4-15.local>
// Created on 30 Mar 2004 by Boyce Griffith (boyce@trasnaform.speakeasy.net)

#include "HierarchyProjector.h"

/////////////////////////////// INCLUDES /////////////////////////////////////

#ifndef included_IBAMR_config
#include <IBAMR_config.h>
#define included_IBAMR_config
#endif

#ifndef included_SAMRAI_config
#include <SAMRAI_config.h>
#define included_SAMRAI_config
#endif

// STOOLS INCLUDES
#include <stools/KrylovLinearSolver.h>
#include <stools/FACPreconditionerLSWrapper.h>
#include <stools/PETScKrylovLinearSolver.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>
#include <IntVector.h>
#include <Patch.h>
#include <RefineOperator.h>
#include <Variable.h>
#include <VariableDatabase.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_project_hierarchy_face;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_project_hierarchy_side;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_hierarchy_configuration;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;

// Version of HierarchyProjector restart file data
static const int HIERARCHY_PROJECTOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

HierarchyProjector::HierarchyProjector(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    bool register_for_restart)
    : d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_do_log(false),
      d_hierarchy(hierarchy),
      d_grid_geom(d_hierarchy->getGridGeometry()),
      d_hier_cc_data_ops(NULL),
      d_hier_fc_data_ops(NULL),
      d_hier_sc_data_ops(NULL),
      d_hier_math_ops(NULL),
      d_is_managing_hier_math_ops(false),
      d_wgt_var(NULL),
      d_wgt_idx(-1),
      d_volume(0.0),
      d_context(NULL),
      d_div_w_var(NULL),
      d_div_w_idx(-1),
      d_sol_vec(NULL),
      d_rhs_vec(NULL),
      d_max_iterations(50),
      d_abs_residual_tol(1.0e-12),
      d_rel_residual_tol(1.0e-8),
      d_poisson_spec(d_object_name+"::Poisson spec"),
      d_default_bc_coef(new SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>(
                            d_object_name+"::default_bc_coef",
                            SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL))),
      d_bc_coef(NULL),
      d_poisson_solver(NULL),
      d_laplace_op(NULL),
      d_poisson_fac_op(NULL),
      d_poisson_fac_pc(NULL),
      d_sol_var(NULL),
      d_rhs_var(NULL),
      d_sol_idx(-1),
      d_rhs_idx(-1)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
    assert(!hierarchy.isNull());
#endif

    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->
            registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the input and restart
    // databases.
    bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Setup a default bc strategy object that specifies homogeneous
    // Neumann boundary conditions.
    for (int d = 0; d < NDIM; ++d)
    {
        d_default_bc_coef->setBoundarySlope(2*d  ,0.0);
        d_default_bc_coef->setBoundarySlope(2*d+1,0.0);
    }

    // Initialize the boundary conditions objects.
    setPhysicalBcCoef(d_default_bc_coef);

    // Get initialization data for the FAC ops and FAC preconditioners
    // and initialize them.
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> fac_op_db, fac_pc_db;

    if (input_db->keyExists("FACOp"))
    {
        fac_op_db = input_db->getDatabase("FACOp");
    }
    if (input_db->keyExists("FACPreconditioner"))
    {
        fac_pc_db = input_db->getDatabase("FACPreconditioner");
    }

    d_poisson_fac_op = new STOOLS::CCPoissonFACOperator(
        d_object_name+"::FAC Op", fac_op_db);
    d_poisson_fac_op->setPoissonSpecifications(d_poisson_spec);
    d_poisson_fac_op->setPhysicalBcCoef(d_bc_coef);

    d_poisson_fac_pc = new SAMRAI::solv::FACPreconditioner<NDIM>(
        d_object_name+"::FAC Preconditioner", *d_poisson_fac_op, fac_pc_db);
    d_poisson_fac_op->setPreconditioner(d_poisson_fac_pc);

    // Initialize the Poisson solver.
    static const bool homogeneous_bc = false;
    d_laplace_op = new STOOLS::CCLaplaceOperator(
        d_object_name+"::Laplace Operator", &d_poisson_spec,
        d_bc_coef, homogeneous_bc);

    d_poisson_solver = new STOOLS::PETScKrylovLinearSolver(
        d_object_name+"::PETSc Krylov solver", "proj_");
    d_poisson_solver->setMaxIterations(d_max_iterations);
    d_poisson_solver->setAbsoluteTolerance(d_abs_residual_tol);
    d_poisson_solver->setRelativeTolerance(d_rel_residual_tol);
    d_poisson_solver->setInitialGuessNonzero(true);
    d_poisson_solver->setOperator(d_laplace_op);
    d_poisson_solver->setPreconditioner(
        new STOOLS::FACPreconditionerLSWrapper(d_poisson_fac_pc, fac_pc_db));

    // Initialize Variables and contexts.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    d_context = var_db->getContext(d_object_name+"::CONTEXT");

    d_sol_var   = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::sol",1);
    d_rhs_var   = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::rhs",1);
    d_div_w_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::div w",1);

    const SAMRAI::hier::IntVector<NDIM> ghosts = 1;
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    d_sol_idx = var_db->registerVariableAndContext(
        d_sol_var, d_context, ghosts);
    d_rhs_idx = var_db->registerVariableAndContext(
        d_rhs_var, d_context, ghosts);
    d_div_w_idx = var_db->registerVariableAndContext(
        d_div_w_var, d_context, ghosts);

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager =
        SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var = new SAMRAI::pdat::CellVariable<NDIM,double>(
        "cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(
        cc_var, hierarchy);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > fc_var = new SAMRAI::pdat::FaceVariable<NDIM,double>(
        "fc_var");
    d_hier_fc_data_ops = hier_ops_manager->getOperationsDouble(
        fc_var, hierarchy);

    SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> > sc_var = new SAMRAI::pdat::SideVariable<NDIM,double>(
        "sc_var");
    d_hier_sc_data_ops = hier_ops_manager->getOperationsDouble(
        sc_var, hierarchy);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_project_hierarchy_face = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::HierarchyProjector::projectHierarchy[Face->Face]");
        t_project_hierarchy_side = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::HierarchyProjector::projectHierarchy[Side->Side]");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::HierarchyProjector::initializeLevelData()");
        t_reset_hierarchy_configuration = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::HierarchyProjector::resetHierarchyConfiguration()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::HierarchyProjector::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// HierarchyProjector

HierarchyProjector::~HierarchyProjector()
{
    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    // Deallocate solver.
    d_poisson_solver->deallocateSolverState();
    d_poisson_solver.setNull();

    // Deallocate other components.
    delete d_default_bc_coef;

    return;
}// ~HierarchyProjector

const std::string&
HierarchyProjector::getName() const
{
    return d_object_name;
}// getName

///
///  The following routines:
///
///      getHierarchyMathOps(),
///      setHierarchyMathOps(),
///      isManagingHierarchyMathOps()
///
///  allow for the sharing of a single HierarchyMathOps object between
///  mutiple HierarchyIntegrator objects.
///

SAMRAI::tbox::Pointer<STOOLS::HierarchyMathOps>
HierarchyProjector::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!d_hier_math_ops.isNull());
#endif
    return d_hier_math_ops;
}// getHierarchyMathOps

void
HierarchyProjector::setHierarchyMathOps(
    SAMRAI::tbox::Pointer<STOOLS::HierarchyMathOps> hier_math_ops,
    const bool manage_ops)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hier_math_ops.isNull());
#endif
    d_hier_math_ops = hier_math_ops;
    d_is_managing_hier_math_ops = manage_ops;
    return;
}// setHierarchyMathOps

bool
HierarchyProjector::isManagingHierarchyMathOps() const
{
    return d_is_managing_hier_math_ops;
}// isManagingHierarchyMathOps

///
///  The following routines:
///
///      setPhysicalBcCoef(),
///      getPhysicalBcCoef(),
///      getPoissonSolver()
///
///  allow other objects to access the Poisson solver and related data
///  used by this integrator.
///

void
HierarchyProjector::setPhysicalBcCoef(
    const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const bc_coef)
{
    if (bc_coef != NULL)
    {
        d_bc_coef = bc_coef;
    }
    else
    {
        d_bc_coef = d_default_bc_coef;
    }

    if (!d_poisson_fac_op.isNull()) d_poisson_fac_op->setPhysicalBcCoef(d_bc_coef);
    if (!d_laplace_op.isNull()) d_laplace_op->setPhysicalBcCoef(d_bc_coef);
    return;
}// setPhysicalBcCoef

const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*
HierarchyProjector::getPhysicalBcCoef() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_bc_coef != NULL);
#endif
    return d_bc_coef;
}// getPhysicalBcCoef

SAMRAI::tbox::Pointer<STOOLS::KrylovLinearSolver>
HierarchyProjector::getPoissonSolver() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!d_poisson_solver.isNull());
#endif
    return d_poisson_solver;
}// getPoissonSolver

///
///  The following routines:
///
///      projectHierarchy()
///
///  provide the projection functionality.
///

void
HierarchyProjector::projectHierarchy(
    const int u_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& u_var,
    const int Phi_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Phi_var,
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& Phi_bdry_fill,
    const double Phi_bdry_fill_time,
    const int grad_Phi_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& grad_Phi_var,
    const int w_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& w_var,
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& w_bdry_fill,
    const double w_bdry_fill_time,
    const bool w_cf_bdry_synch,
    const int Q_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Q_var)
{
    t_project_hierarchy_face->start();

    d_laplace_op->setTime(Phi_bdry_fill_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Compute div w.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_div_w_idx, w_bdry_fill_time);
    }

    d_hier_math_ops->div(
        d_div_w_idx, d_div_w_var, // dst
        -1.0,                     // alpha
        w_idx, w_var,             // src1
        w_bdry_fill,              // src1_bdry_fill
        w_bdry_fill_time,         // src1_bdry_fill_time
        w_cf_bdry_synch,          // src1_cf_bdry_synch
        1.0,                      // beta
        Q_idx, Q_var);            // src2

    // Solve -div grad Phi = div u - div w.
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double> sol_vec(
        d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(Phi_var, Phi_idx, d_wgt_idx, d_hier_cc_data_ops);

    SAMRAI::solv::SAMRAIVectorReal<NDIM,double> rhs_vec(
        d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_div_w_var, d_div_w_idx, d_wgt_idx, d_hier_cc_data_ops);

    if (d_do_log) SAMRAI::tbox::plog << "HierarchyProjector::projectHierarchy(): about to solve projection equation . . . \n";
    d_poisson_solver->solveSystem(sol_vec,rhs_vec);
    if (d_do_log) SAMRAI::tbox::plog << "HierarchyProjector::projectHierarchy(): number of iterations = " << d_poisson_solver->getNumIterations() << "\n";
    if (d_do_log) SAMRAI::tbox::plog << "HierarchyProjector::projectHierarchy(): residual norm        = " << d_poisson_solver->getResidualNorm()  << "\n";

    // Deallocate div w.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_div_w_idx);
    }

    // Normalize Phi.
    const double Phi_mean = d_hier_cc_data_ops->integral(
        Phi_idx, d_wgt_idx)/d_volume;
    d_hier_cc_data_ops->addScalar(Phi_idx, Phi_idx, -Phi_mean);

    // Set u = w - grad Phi.
    const bool grad_Phi_cf_bdry_synch = true;

    d_hier_math_ops->grad(
        grad_Phi_idx, grad_Phi_var, // dst
        grad_Phi_cf_bdry_synch,     // dst_cf_bdry_synch
        1.0,                        // alpha
        Phi_idx, Phi_var,           // src
        Phi_bdry_fill,              // src_bdry_fill
        Phi_bdry_fill_time);        // src_bdry_fill_time

    d_hier_fc_data_ops->subtract(u_idx, w_idx, grad_Phi_idx);

    t_project_hierarchy_face->stop();
    return;
}// projectHierarchy

void
HierarchyProjector::projectHierarchy(
    const int u_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& u_var,
    const int Phi_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Phi_var,
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& Phi_bdry_fill,
    const double Phi_bdry_fill_time,
    const int grad_Phi_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& grad_Phi_var,
    const int w_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::SideVariable<NDIM,double> >& w_var,
    const std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >& w_bdry_fill,
    const double w_bdry_fill_time,
    const bool w_cf_bdry_synch,
    const int Q_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Q_var)
{
    t_project_hierarchy_side->start();

    d_laplace_op->setTime(Phi_bdry_fill_time);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Compute div w.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_div_w_idx, w_bdry_fill_time);
    }

    d_hier_math_ops->div(
        d_div_w_idx, d_div_w_var, // dst
        -1.0,                     // alpha
        w_idx, w_var,             // src1
        w_bdry_fill,              // src1_bdry_fill
        w_bdry_fill_time,         // src1_bdry_fill_time
        w_cf_bdry_synch,          // src1_cf_bdry_synch
        1.0,                      // beta
        Q_idx, Q_var);            // src2

    // Solve -div grad Phi = div u - div w.
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double> sol_vec(
        d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(Phi_var, Phi_idx, d_wgt_idx, d_hier_cc_data_ops);

    SAMRAI::solv::SAMRAIVectorReal<NDIM,double> rhs_vec(
        d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_div_w_var, d_div_w_idx, d_wgt_idx, d_hier_cc_data_ops);

    if (d_do_log) SAMRAI::tbox::plog << "HierarchyProjector::projectHierarchy(): about to solve projection equation . . . \n";
    d_poisson_solver->solveSystem(sol_vec,rhs_vec);
    if (d_do_log) SAMRAI::tbox::plog << "HierarchyProjector::projectHierarchy(): number of iterations = " << d_poisson_solver->getNumIterations() << "\n";
    if (d_do_log) SAMRAI::tbox::plog << "HierarchyProjector::projectHierarchy(): residual norm        = " << d_poisson_solver->getResidualNorm()  << "\n";

    // Deallocate div w.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_div_w_idx);
    }

    // Normalize Phi.
    const double Phi_mean = d_hier_cc_data_ops->integral(
        Phi_idx, d_wgt_idx)/d_volume;
    d_hier_cc_data_ops->addScalar(Phi_idx, Phi_idx, -Phi_mean);

    // Set u = w - grad Phi.
    const bool grad_Phi_cf_bdry_synch = true;

    d_hier_math_ops->grad(
        grad_Phi_idx, grad_Phi_var, // dst
        grad_Phi_cf_bdry_synch,     // dst_cf_bdry_synch
        1.0,                        // alpha
        Phi_idx, Phi_var,           // src
        Phi_bdry_fill,              // src_bdry_fill
        Phi_bdry_fill_time);        // src_bdry_fill_time

    d_hier_fc_data_ops->subtract(u_idx, w_idx, grad_Phi_idx);

    t_project_hierarchy_side->stop();
    return;
}// projectHierarchy

///
///  The following routines:
///
///      initializeLevelData(),
///      resetHierarchyConfiguration()
///
///  are concrete implementations of functions declared in the
///  mesh::StandardTagAndInitStrategy<NDIM> abstract base class.
///

void
HierarchyProjector::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((level_number >= 0)
           && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        assert(level_number == old_level->getLevelNumber());
    }
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif

    // intentionally blank

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
HierarchyProjector::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((coarsest_level >= 0)
           && (coarsest_level <= finest_level)
           && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        assert(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // Reset the Hierarchy data operations for the new hierarchy
    // configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_fc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_fc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_sc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_sc_data_ops->resetLevels(0, finest_hier_level);

    // Reset the Hierarchy math operations for the new configuration.
    if (d_is_managing_hier_math_ops)
    {
        d_hier_math_ops->setPatchHierarchy(hierarchy);
        d_hier_math_ops->resetLevels(0, finest_hier_level);
    }
    else if (d_hier_math_ops.isNull())
    {
        d_hier_math_ops = new STOOLS::HierarchyMathOps(
            d_object_name+"::HierarchyMathOps", hierarchy);
        d_is_managing_hier_math_ops = true;
    }

    // Get the cell weights data.
    d_wgt_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Get the volume of the physical domain.
    d_volume = d_hier_math_ops->getVolumeOfPhysicalDomain();

    // Reset the solution and rhs vectors.
    d_sol_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::sol_vec", d_hierarchy, 0, finest_hier_level);
    d_sol_vec->addComponent(d_sol_var,d_sol_idx,d_wgt_idx,d_hier_cc_data_ops);

    d_rhs_vec = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
        d_object_name+"::rhs_vec", d_hierarchy, 0, finest_hier_level);
    d_rhs_vec->addComponent(d_rhs_var,d_rhs_idx,d_wgt_idx,d_hier_cc_data_ops);

    // (Re)-initialize the Poisson solver.
    d_laplace_op->setHierarchyMathOps(d_hier_math_ops);
    d_poisson_fac_op->setResetLevels(coarsest_level, finest_level);
    d_poisson_solver->initializeSolverState(*d_sol_vec,*d_rhs_vec);

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::tbox::Serializable abstract base class.
///

void
HierarchyProjector::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif
    db->putInteger("HIERARCHY_PROJECTOR_VERSION", HIERARCHY_PROJECTOR_VERSION);

    t_put_to_database->stop();
    return;
}// putToDatabase

///
///  The following routines:
///
///      printClassData()
///
///  are provided for your viewing pleasure.
///

void
HierarchyProjector::printClassData(
    std::ostream& os) const
{
    os << "\nHierarchyProjector::printClassData..." << endl;
    os << "\nHierarchyProjector: this = " << const_cast<HierarchyProjector*>(this) << endl;
    os << "d_object_name = " << d_object_name << endl;
    return;
}// printClassData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
HierarchyProjector::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif
    (void) is_from_restart;
    d_max_iterations = db->getIntegerWithDefault("max_iterations", d_max_iterations);
    d_abs_residual_tol = db->getDoubleWithDefault("abs_residual_tol", d_abs_residual_tol);
    d_rel_residual_tol = db->getDoubleWithDefault("rel_residual_tol", d_rel_residual_tol);
    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);
    return;
}// getFromInput

void
HierarchyProjector::getFromRestart()
{
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> restart_db =
        SAMRAI::tbox::RestartManager::getManager()->getRootDatabase();

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  \n"
                   << "Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("HIERARCHY_PROJECTOR_VERSION");
    if (ver != HIERARCHY_PROJECTOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::HierarchyProjector>;

//////////////////////////////////////////////////////////////////////////////
