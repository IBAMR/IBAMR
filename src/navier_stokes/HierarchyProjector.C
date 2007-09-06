// Filename: HierarchyProjector.C
// Last modified: <06.Sep.2007 03:05:24 griffith@box221.cims.nyu.edu>
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
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_project_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_hierarchy_configuration;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;

// Version of HierarchyProjector restart file data.
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
      d_hier_math_ops(NULL),
      d_is_managing_hier_math_ops(false),
      d_wgt_var(NULL),
      d_wgt_idx(-1),
      d_volume(0.0),
      d_context(NULL),
      d_F_var(NULL),
      d_F_idx(-1),
      d_w_var(NULL),
      d_w_idx(-1),
      d_sol_vec(NULL),
      d_rhs_vec(NULL),
      d_max_iterations(50),
      d_abs_residual_tol(1.0e-12),
      d_rel_residual_tol(1.0e-8),
      d_poisson_spec(d_object_name+"::Poisson spec"),
      d_u_bc_coefs(NDIM,NULL),
      d_default_u_bc_coefs(),
      d_Phi_bc_coef(NULL),
      d_poisson_solver(NULL),
      d_laplace_op(NULL),
      d_poisson_fac_op(NULL),
      d_poisson_fac_pc(NULL),
      d_sol_var(NULL),
      d_rhs_var(NULL),
      d_sol_idx(-1),
      d_rhs_idx(-1),
      d_inhomogeneous_bc_fill_alg(NULL),
      d_inhomogeneous_bc_fill_scheds(),
      d_bc_refine_strategy(NULL),
      d_bc_op(NULL)
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

    // Initialize object with data read from the input and restart databases.
    bool from_restart = SAMRAI::tbox::RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Initialize Variables and contexts.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    d_context = var_db->getContext(d_object_name+"::CONTEXT");

    d_sol_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::sol",1);
    d_rhs_var = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::rhs",1);
    d_F_var   = new SAMRAI::pdat::CellVariable<NDIM,double>(d_object_name+"::f",1);
    d_w_var   = new SAMRAI::pdat::FaceVariable<NDIM,double>(d_object_name+"::w",1);

    const SAMRAI::hier::IntVector<NDIM> ghosts = 1;
    const SAMRAI::hier::IntVector<NDIM> no_ghosts = 0;

    d_sol_idx = var_db->registerVariableAndContext(d_sol_var, d_context, ghosts);
    d_rhs_idx = var_db->registerVariableAndContext(d_rhs_var, d_context, ghosts);
    d_F_idx   = var_db->registerVariableAndContext(  d_F_var, d_context, ghosts);
    d_w_idx   = var_db->registerVariableAndContext(  d_w_var, d_context, ghosts);

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager =
        SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var =
        new SAMRAI::pdat::CellVariable<NDIM,double>("cc_var");
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > fc_var =
        new SAMRAI::pdat::FaceVariable<NDIM,double>("fc_var");

    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy);
    d_hier_fc_data_ops = hier_ops_manager->getOperationsDouble(fc_var, hierarchy);

    // Initialize the Poisson specifications.
    d_poisson_spec.setCZero();
    d_poisson_spec.setDConstant(-1.0);

    // Setup default velocity boundary condition objects to specify homogeneous
    // Dirichlet boundary conditions for all components of the velocity.
    for (int d = 0; d < NDIM; ++d)
    {
        std::ostringstream stream;
        stream << d;
        SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>* u_bc_coef = new SAMRAI::solv::LocationIndexRobinBcCoefs<NDIM>(
            d_object_name + "::default_u_bc_coef_" + stream.str(),
            SAMRAI::tbox::Pointer<SAMRAI::tbox::Database>(NULL));

        for (int location_index = 0; location_index < 2*NDIM; ++location_index)
        {
            u_bc_coef->setBoundaryValue(location_index,0.0);
        }

        d_default_u_bc_coefs.push_back(u_bc_coef);
    }

    // Initialize the boundary conditions objects.
    d_Phi_bc_coef = new INSProjectionBcCoef(d_w_idx,d_default_u_bc_coefs,true);
    setVelocityPhysicalBcCoefs(d_default_u_bc_coefs);

    // Get initialization data for the FAC ops and FAC preconditioners and
    // initialize them.
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
    d_poisson_fac_op->setPhysicalBcCoef(d_Phi_bc_coef);

    d_poisson_fac_pc = new SAMRAI::solv::FACPreconditioner<NDIM>(
        d_object_name+"::FAC Preconditioner", *d_poisson_fac_op, fac_pc_db);
    d_poisson_fac_op->setPreconditioner(d_poisson_fac_pc);

    // Initialize the Poisson solver.
    static const bool homogeneous_bc = false;
    d_laplace_op = new STOOLS::CCLaplaceOperator(
        d_object_name+"::Laplace Operator",
        d_poisson_spec, d_Phi_bc_coef, homogeneous_bc);

    d_poisson_solver = new STOOLS::PETScKrylovLinearSolver(d_object_name+"::PETSc Krylov solver", "proj_");
    d_poisson_solver->setMaxIterations(d_max_iterations);
    d_poisson_solver->setAbsoluteTolerance(d_abs_residual_tol);
    d_poisson_solver->setRelativeTolerance(d_rel_residual_tol);
    d_poisson_solver->setInitialGuessNonzero(true);
    d_poisson_solver->setOperator(d_laplace_op);
    d_poisson_solver->setPreconditioner(new STOOLS::FACPreconditionerLSWrapper(d_poisson_fac_pc, fac_pc_db));

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_project_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::HierarchyProjector::projectHierarchy");
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
    delete d_Phi_bc_coef;
    for (int d = 0; d < NDIM; ++d)
    {
        delete d_default_u_bc_coefs[d];
    }
    d_default_u_bc_coefs.clear();

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
///  allow for the sharing of a single HierarchyMathOps object between mutiple
///  HierarchyIntegrator objects.
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
///      setVelocityPhysicalBcCoefs(),
///      getVelocityPhysicalBcCoefs(),
///      getPoissonSolver()
///
///  allow other objects to access the Poisson solver and related data used by
///  this integrator.
///

void
HierarchyProjector::setVelocityPhysicalBcCoefs(
    const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& u_bc_coefs)
{
    if (u_bc_coefs.size() != NDIM)
    {
        TBOX_ERROR(d_object_name << "::setVelocityPhysicalBcCoefs():\n"
                   << "  precisely NDIM boundary condiiton objects must be provided." << endl);
    }
#ifdef DEBUG_CHECK_ASSERTIONS
    for (unsigned l = 0; l < u_bc_coefs.size(); ++l)
    {
        assert(u_bc_coefs[l] != NULL);
    }
#endif
    d_u_bc_coefs = u_bc_coefs;
    d_Phi_bc_coef->setVelocityPhysicalBcCoefs(d_u_bc_coefs);
    return;
}// setVelocityPhysicalBcCoefs

const std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>&
HierarchyProjector::getVelocityPhysicalBcCoefs() const
{
    return d_u_bc_coefs;
}// getVelocityPhysicalBcCoefs

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
    const double rho,
    const double dt,
    const double time,
    const std::string& projection_type,
    const int u_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& u_var,
    const int Phi_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Phi_var,
    const int grad_Phi_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& grad_Phi_var,
    const int w_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> >& w_var,
    const int Q_idx,
    const SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> >& Q_var)
{
    t_project_hierarchy->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > > no_fill(finest_ln+1,NULL);

    // Allocate temporary data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_F_idx, time);
        level->allocatePatchData(d_w_idx, time);
    }

    // Setup the boundary coefficient specification object.
    //
    // NOTE: These boundary coefficients are also used by the linear operator
    // and by the FAC preconditioner.
    d_Phi_bc_coef->setProblemCoefs(rho, dt);
    d_Phi_bc_coef->setIntermediateVelocityPatchDataIndex(d_w_idx);
    d_Phi_bc_coef->setVelocityPhysicalBcCoefs(d_u_bc_coefs);

    // Setup the linear operator.
    d_laplace_op->setTime(time);
    d_laplace_op->setHierarchyMathOps(d_hier_math_ops);

    // Setup the preconditioner.
    d_poisson_fac_op->setTime(time);

    // Compute F = (rho/dt)*(Q - div w).
    const bool w_cf_bdry_synch = true;
    d_hier_math_ops->div(d_F_idx, d_F_var, // dst
                         -rho/dt,          // alpha
                         w_idx, w_var,     // src1
                         no_fill,          // src1_bdry_fill
                         time,             // src1_bdry_fill_time
                         w_cf_bdry_synch,  // src1_cf_bdry_synch
                         +rho/dt,          // beta
                         Q_idx, Q_var);    // src2

    // Copy the intermediate velocity w into the scratch data.
    d_hier_fc_data_ops->copyData(d_w_idx, w_idx, false);

    // Solve -div grad Phi = F = (rho/dt)*(Q - div w).
    SAMRAI::solv::SAMRAIVectorReal<NDIM,double> sol_vec(
        d_object_name+"::sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    sol_vec.addComponent(Phi_var, Phi_idx, d_wgt_idx, d_hier_cc_data_ops);

    SAMRAI::solv::SAMRAIVectorReal<NDIM,double> rhs_vec(
        d_object_name+"::rhs_vec", d_hierarchy, coarsest_ln, finest_ln);
    rhs_vec.addComponent(d_F_var, d_F_idx, d_wgt_idx, d_hier_cc_data_ops);

    if (d_do_log) SAMRAI::tbox::plog << "HierarchyProjector::projectHierarchy(): about to solve projection equation . . . \n";
    d_poisson_solver->solveSystem(sol_vec,rhs_vec);
    if (d_do_log) SAMRAI::tbox::plog << "HierarchyProjector::projectHierarchy(): number of iterations = " << d_poisson_solver->getNumIterations() << "\n";
    if (d_do_log) SAMRAI::tbox::plog << "HierarchyProjector::projectHierarchy(): residual norm        = " << d_poisson_solver->getResidualNorm()  << "\n";

    if (d_poisson_solver->getNumIterations() == d_poisson_solver->getMaxIterations())
    {
        SAMRAI::tbox::pout << d_object_name << "::projectHierarchy():"
                           <<"  WARNING: linear solver iterations == max iterations\n";
    }

    // Fill the physical boundary conditions for Phi.
    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineAlgorithm<NDIM> > bdry_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > Phi_bc_fill_op = d_grid_geom->
        lookupRefineOperator(Phi_var, "CONSTANT_REFINE");

    bdry_fill_alg->registerRefine(Phi_idx,  // destination
                                  Phi_idx,  // source
                                  Phi_idx,  // temporary work space
                                  Phi_bc_fill_op);

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > w_bc_fill_op = d_grid_geom->
        lookupRefineOperator(d_w_var, "CONSTANT_REFINE");

    bdry_fill_alg->registerRefine(d_w_idx,  // destination
                                  d_w_idx,  // source
                                  d_w_idx,  // temporary work space
                                  w_bc_fill_op);

    d_bc_op->setPatchDataIndex(Phi_idx);
    d_bc_op->setPhysicalBcCoef(d_Phi_bc_coef);
    d_bc_op->setHomogeneousBc(false);
    d_Phi_bc_coef->setHomogeneousBc(false);

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        bdry_fill_alg->resetSchedule(d_inhomogeneous_bc_fill_scheds[ln]);
        d_inhomogeneous_bc_fill_scheds[ln]->fillData(time);
        d_inhomogeneous_bc_fill_alg->resetSchedule(d_inhomogeneous_bc_fill_scheds[ln]);
    }

    d_Phi_bc_coef->setHomogeneousBc(true);

    // Set u = w - (dt/rho)*grad Phi.
    const bool grad_Phi_cf_bdry_synch = true;
    d_hier_math_ops->grad(grad_Phi_idx, grad_Phi_var, // dst
                          grad_Phi_cf_bdry_synch,     // dst_cf_bdry_synch
                          1.0,                        // alpha
                          Phi_idx, Phi_var,           // src
                          no_fill,                    // src_bdry_fill
                          time);                      // src_bdry_fill_time
    d_hier_fc_data_ops->axpy(u_idx, -dt/rho, grad_Phi_idx, w_idx);

    // Deallocate scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_F_idx);
        level->deallocatePatchData(d_w_idx);
    }

    t_project_hierarchy->stop();
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

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    d_hier_fc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_fc_data_ops->resetLevels(0, finest_hier_level);

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

    // Setup the communications algorithm.
    d_inhomogeneous_bc_fill_alg = new SAMRAI::xfer::RefineAlgorithm<NDIM>();

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > F_bc_fill_op = d_grid_geom->
        lookupRefineOperator(d_F_var, "CONSTANT_REFINE");
    d_inhomogeneous_bc_fill_alg->registerRefine(d_F_idx,  // destination
                                                d_F_idx,  // source
                                                d_F_idx,  // temporary work space
                                                F_bc_fill_op);

    SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > w_bc_fill_op = d_grid_geom->
        lookupRefineOperator(d_w_var, "CONSTANT_REFINE");
    d_inhomogeneous_bc_fill_alg->registerRefine(d_w_idx,  // destination
                                                d_w_idx,  // source
                                                d_w_idx,  // temporary work space
                                                w_bc_fill_op);

    d_bc_op = new STOOLS::CartRobinPhysBdryOp(d_F_idx, d_Phi_bc_coef, false);
    d_bc_refine_strategy = d_bc_op;

    // Setup the communications schedules.
    d_inhomogeneous_bc_fill_scheds.resize(finest_hier_level+1);
    for (int ln = 0; ln <= finest_hier_level; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        d_inhomogeneous_bc_fill_scheds[ln] = d_inhomogeneous_bc_fill_alg->
            createSchedule(level, ln-1, d_hierarchy, d_bc_refine_strategy);
    }

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
    os << "this = " << const_cast<HierarchyProjector*>(this) << endl;
    os << "d_object_name = " << d_object_name << "\n"
       << "d_registered_for_restart = " << d_registered_for_restart << endl;
    os << "d_do_log = " << d_do_log << endl;
    os << "d_hierarchy = " << d_hierarchy.getPointer() << "\n"
       << "d_grid_geom = " << d_grid_geom.getPointer() << endl;
    os << "d_hier_cc_data_ops = " << d_hier_cc_data_ops.getPointer() << "\n"
       << "d_hier_fc_data_ops = " << d_hier_fc_data_ops.getPointer() << "\n"
       << "d_hier_math_ops = " << d_hier_math_ops.getPointer() << "\n"
       << "d_is_managing_hier_math_ops = " << d_is_managing_hier_math_ops << endl;
    os << "d_wgt_var = " << d_wgt_var.getPointer() << "\n"
       << "d_wgt_idx = " << d_wgt_idx << "\n"
       << "d_volume = " << d_volume << endl;
    os << "Skipping variables, patch data descriptors, communications algorithms, etc." << endl;
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
