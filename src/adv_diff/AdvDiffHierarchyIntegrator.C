// Filename: AdvDiffHierarchyIntegrator.C
// Last modified: <07.Apr.2007 16:04:42 griffith@box221.cims.nyu.edu>
// Created on 17 Mar 2004 by Boyce Griffith (boyce@bigboy.speakeasy.net)

#include "AdvDiffHierarchyIntegrator.h"

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
#include <stools/FACPreconditionerLSWrapper.h>
#include <stools/PETScKrylovLinearSolver.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CellDataFactory.h>
#include <CoarsenOperator.h>
#include <HierarchyDataOpsManager.h>
#include <RefineOperator.h>
#include <VariableDatabase.h>
#include <tbox/NullDatabase.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <cassert>
#include <limits>
#include <set>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy_integrator;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_advance_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_regrid_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_integrate_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_hierarchy;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_synchronize_new_levels;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_time_dependent_hier_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_hier_data_to_preadvance_state;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_initialize_level_data;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_reset_hierarchy_configuration;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_apply_gradient_detector;
static SAMRAI::tbox::Pointer<SAMRAI::tbox::Timer> t_put_to_database;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;

// Version of AdvDiffHierarchyIntegrator restart file data.
static const int ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffHierarchyIntegrator::AdvDiffHierarchyIntegrator(
    const std::string& object_name,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    SAMRAI::tbox::Pointer<GodunovAdvector> explicit_predictor,
    bool register_for_restart)
    : d_Q_vars(),
      d_F_vars(),
      d_Psi_vars(),
      d_grad_vars(),
      d_Q_inits(),
      d_Q_bc_coefs(),
      d_F_sets(),
      d_Q_mus(),
      d_u_var(NULL),
      d_u_set(NULL),
      d_u_is_div_free(false),
      d_object_name(object_name),
      d_registered_for_restart(register_for_restart),
      d_hierarchy(hierarchy),
      d_gridding_alg(NULL),
      d_hyp_level_integrator(NULL),
      d_hyp_patch_ops(NULL),
      d_start_time(0.0),
      d_end_time(std::numeric_limits<double>::max()),
      d_grow_dt(2.0),
      d_max_integrator_steps(std::numeric_limits<int>::max()),
      d_regrid_interval(1),
      d_using_default_tag_buffer(true),
      d_tag_buffer(),
      d_old_dt(-1.0),
      d_integrator_time(std::numeric_limits<double>::quiet_NaN()),
      d_integrator_step(std::numeric_limits<int>::max()),
      d_is_initialized(false),
      d_do_log(false),
      d_hier_cc_data_ops(NULL),
      d_hier_math_ops(NULL),
      d_is_managing_hier_math_ops(true),
      d_wgt_var(NULL),
      d_wgt_idx(-1),
      d_temp_context(),
      d_ralgs(),
      d_rstrategies(),
      d_rscheds(),
      d_calgs(),
      d_cstrategies(),
      d_cscheds(),
      d_bc_op(new STOOLS::CartRobinPhysBdryOp()),
      d_sol_vecs(),
      d_rhs_vecs(),
      d_max_iterations(25),
      d_abs_residual_tol(1.0e-30),
      d_rel_residual_tol(1.0e-8),
      d_using_FAC(true),
      d_helmholtz_ops(),
      d_helmholtz_specs(),
      d_helmholtz_solvers(),
      d_helmholtz_fac_ops(),
      d_helmholtz_fac_pcs(),
      d_helmholtz_solvers_need_init(),
      d_coarsest_reset_ln(-1),
      d_finest_reset_ln(-1),
      d_fac_ops_db(NULL),
      d_fac_pcs_db(NULL)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!object_name.empty());
    assert(!input_db.isNull());
    assert(!hierarchy.isNull());
    assert(!explicit_predictor.isNull());
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

    // Initialize the SAMRAI::algs::HyperbolicPatchStrategy and
    // SAMRAI::algs::HyperbolicLevelIntegrator objects used to provide
    // the numerical routines for explicitly integrating the advective
    // terms.
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> hyp_patch_ops_input_db;
    if (input_db->keyExists("AdvDiffHypPatchOps"))
    {
        hyp_patch_ops_input_db = input_db->getDatabase("AdvDiffHypPatchOps");
    }
    else
    {
        hyp_patch_ops_input_db = new SAMRAI::tbox::NullDatabase();
    }

    d_hyp_patch_ops = new AdvDiffHypPatchOps(
        object_name+"::AdvDiffHypPatchOps",
        hyp_patch_ops_input_db,
        explicit_predictor,
        d_hierarchy->getGridGeometry(),
        d_registered_for_restart);

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> hyp_level_integrator_input_db;
    if (input_db->keyExists("HyperbolicLevelIntegrator"))
    {
        hyp_level_integrator_input_db = input_db->getDatabase("HyperbolicLevelIntegrator");
    }
    else
    {
        hyp_level_integrator_input_db = new SAMRAI::tbox::NullDatabase();
    }

    static const bool using_time_refinement = false;
    d_hyp_level_integrator = new SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>(
        object_name+"::HyperbolicLevelIntegrator",
        hyp_level_integrator_input_db,
        d_hyp_patch_ops,
        register_for_restart,
        using_time_refinement);

    // Get initialization data for the FAC ops and FAC preconditioners.
    if (d_using_FAC)
    {
        if (input_db->keyExists("FACOps"))
        {
            d_fac_ops_db = input_db->getDatabase("FACOps");
        }
        if (input_db->keyExists("FACPreconditioners"))
        {
            d_fac_pcs_db = input_db->getDatabase("FACPreconditioners");
        }
    }

    // Obtain the Hierarchy data operations objects.
    SAMRAI::math::HierarchyDataOpsManager<NDIM>* hier_ops_manager =
        SAMRAI::math::HierarchyDataOpsManager<NDIM>::getManager();

    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > cc_var =
        new SAMRAI::pdat::CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(
        cc_var, hierarchy);

    // Setup Timers.
    static bool timers_need_init = true;
    if (timers_need_init)
    {
        t_initialize_hierarchy_integrator = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::advanceHierarchy()");
        t_regrid_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::regridHierarchy()");
        t_integrate_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::integrateHierarchy()");
        t_synchronize_hierarchy = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_hier_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_hier_data_to_preadvance_state = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database = SAMRAI::tbox::TimerManager::getManager()->
            getTimer("IBAMR::AdvDiffHierarchyIntegrator::putToDatabase()");
        timers_need_init = false;
    }
    return;
}// AdvDiffHierarchyIntegrator

AdvDiffHierarchyIntegrator::~AdvDiffHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        SAMRAI::tbox::RestartManager::getManager()->unregisterRestartItem(d_object_name);
    }

    for (RefinePatchStrategyMap::iterator it = d_rstrategies.begin();
         it != d_rstrategies.end(); ++it)
    {
        if ((*it).second != d_bc_op.getPointer()) delete (*it).second;
    }

    for (CoarsenPatchStrategyMap::iterator it = d_cstrategies.begin();
         it != d_cstrategies.end(); ++it)
    {
        delete (*it).second;
    }

    // Deallocate all solver components.
    //
    // NOTE: The following code ensures that the solver components are
    // deallocated in the correct order.
    d_helmholtz_solvers.clear();
    d_helmholtz_fac_pcs.clear();
    d_helmholtz_fac_ops.clear();
    d_helmholtz_ops.clear();
    d_helmholtz_specs.clear();

    return;
}// ~AdvDiffHierarchyIntegrator

const std::string&
AdvDiffHierarchyIntegrator::getName() const
{
    return d_object_name;
}// getName

void
AdvDiffHierarchyIntegrator::registerVisItDataWriter(
    SAMRAI::tbox::Pointer<SAMRAI::appu::VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!visit_writer.isNull());
#endif
    d_hyp_patch_ops->registerVisItDataWriter(visit_writer);
    return;
}// registerVisItDataWriter

///
///  The following routines:
///
///      registerAdvectedAndDiffusedQuantity(),
///      registerAdvectedAndDiffusedQuantityWithSourceTerm(),
///      registerAdvectionVelocity()
///
///  allow the specification of quantities to be advected and
///  diffused.
///

void
AdvDiffHierarchyIntegrator::registerAdvectedAndDiffusedQuantity(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
    const double Q_mu,
    const bool conservation_form,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> Q_init,
    const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const Q_bc_coef,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var)
{
    registerAdvectedAndDiffusedQuantity(
        Q_var, Q_mu, conservation_form, Q_init,
        std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(1,Q_bc_coef),
        grad_var);
    return;
}// registerAdvectedAndDiffusedQuantity

void
AdvDiffHierarchyIntegrator::registerAdvectedAndDiffusedQuantity(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
    const double Q_mu,
    const bool conservation_form,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> Q_init,
    const std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& Q_bc_coefs,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!Q_var.isNull());
    assert(Q_mu >= 0.0);
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > Q_factory =
        Q_var->getPatchDataFactory();
    const int Q_depth = Q_factory->getDefaultDepth();

    std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> Q_bc_coefs_local = Q_bc_coefs;
    if (Q_bc_coefs_local.empty())
    {
        Q_bc_coefs_local = std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(
            Q_depth,static_cast<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NULL));
    }

    if (Q_depth != static_cast<int>(Q_bc_coefs_local.size()))
    {
        TBOX_ERROR(d_object_name << "::registerAdvectedAndDiffusedQuantity():\n"
                   << "  data depth for variable " << Q_var->getName() << " is " << Q_depth << "\n"
                   << "  but " << Q_bc_coefs_local.size() << " boundary condition coefficient objects were provided to the class constructor." << endl);
    }

    if (!SAMRAI::tbox::Utilities::deq(Q_mu,0.0))
    {
        d_Q_vars    .push_back(Q_var);
        d_Q_inits   .push_back(Q_init);
        d_Q_bc_coefs.push_back(Q_bc_coefs_local);
        d_Q_mus     .push_back(Q_mu);

        d_grad_vars.push_back(grad_var);

        d_F_vars.push_back(NULL);
        d_F_sets.push_back(NULL);

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > Q_factory =
            Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Psi_var =
            new SAMRAI::pdat::CellVariable<NDIM,double>(Q_var->getName()+"::Psi",Q_depth);

        d_Psi_vars.push_back(Psi_var);

        d_hyp_patch_ops->registerAdvectedQuantityWithSourceTerm(
            Q_var, Psi_var, conservation_form, Q_init, Q_bc_coefs_local,
            SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy>(NULL), grad_var);
    }
    else
    {
        d_hyp_patch_ops->registerAdvectedQuantity(
            Q_var, conservation_form, Q_init, Q_bc_coefs_local, grad_var);
    }
    return;
}// registerAdvectedAndDiffusedQuantity

void
AdvDiffHierarchyIntegrator::registerAdvectedAndDiffusedQuantityWithSourceTerm(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
    const double Q_mu,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var,
    const bool conservation_form,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> Q_init,
    const SAMRAI::solv::RobinBcCoefStrategy<NDIM>* const Q_bc_coef,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> F_set,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var)
{
    registerAdvectedAndDiffusedQuantityWithSourceTerm(
        Q_var, Q_mu, F_var, conservation_form, Q_init,
        std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(1,Q_bc_coef),
        F_set, grad_var);
    return;
}// registerAdvectedAndDiffusedQuantityWithSourceTerm

void
AdvDiffHierarchyIntegrator::registerAdvectedAndDiffusedQuantityWithSourceTerm(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var,
    const double Q_mu,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var,
    const bool conservation_form,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> Q_init,
    const std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& Q_bc_coefs,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> F_set,
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > grad_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!Q_var.isNull());
    assert(Q_mu >= 0.0);
    assert(!F_var.isNull());
#endif
    SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > Q_factory =
        Q_var->getPatchDataFactory();
    const int Q_depth = Q_factory->getDefaultDepth();

    std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> Q_bc_coefs_local = Q_bc_coefs;
    if (Q_bc_coefs_local.empty())
    {
        Q_bc_coefs_local = std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(
            Q_depth,static_cast<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>(NULL));
    }

    if (Q_depth != static_cast<int>(Q_bc_coefs_local.size()))
    {
        TBOX_ERROR(d_object_name << "::registerAdvectedAndDiffusedQuantityWithSourceTerm():\n"
                   << "  data depth for variable " << Q_var->getName() << " is " << Q_depth << "\n"
                   << "  but " << Q_bc_coefs_local.size() << " boundary condition coefficient objects were provided to the class constructor." << endl);
    }

    if (!SAMRAI::tbox::Utilities::deq(Q_mu,0.0))
    {
        d_Q_vars    .push_back(Q_var);
        d_Q_inits   .push_back(Q_init);
        d_Q_bc_coefs.push_back(Q_bc_coefs_local);
        d_Q_mus     .push_back(Q_mu);

        d_grad_vars.push_back(grad_var);

        d_F_vars.push_back(F_var);
        d_F_sets.push_back(F_set);

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > Q_factory =
            Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();

#ifdef DEBUG_CHECK_ASSERTIONS
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > F_factory =
            F_var->getPatchDataFactory();
        assert(Q_depth == F_factory->getDefaultDepth());
#endif

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Psi_var =
            new SAMRAI::pdat::CellVariable<NDIM,double>(Q_var->getName()+"::Psi",Q_depth);

        d_Psi_vars.push_back(Psi_var);

        d_hyp_patch_ops->registerAdvectedQuantityWithSourceTerm(
            Q_var, Psi_var, conservation_form, Q_init, Q_bc_coefs_local,
            SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy>(NULL), grad_var);
    }
    else
    {
        d_hyp_patch_ops->registerAdvectedQuantityWithSourceTerm(
            Q_var, F_var, conservation_form, Q_init, Q_bc_coefs_local, F_set, grad_var);
    }
    return;
}// registerAdvectedAndDiffusedQuantityWithSourceTerm

void
AdvDiffHierarchyIntegrator::registerAdvectionVelocity(
    SAMRAI::tbox::Pointer<SAMRAI::pdat::FaceVariable<NDIM,double> > u_var,
    const bool u_is_div_free,
    SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> u_set)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!u_var.isNull());
#endif
    d_hyp_patch_ops->registerAdvectionVelocity(
        u_var,u_is_div_free,u_set);
    d_u_var = u_var;
    d_u_set = u_set;
    d_u_is_div_free = u_is_div_free;
    return;
}// registerAdvectionVelocity

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
AdvDiffHierarchyIntegrator::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!d_hier_math_ops.isNull());
#endif
    return d_hier_math_ops;
}// getHierarchyMathOps

void
AdvDiffHierarchyIntegrator::setHierarchyMathOps(
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
AdvDiffHierarchyIntegrator::isManagingHierarchyMathOps() const
{
    return d_is_managing_hier_math_ops;
}// isManagingHierarchyMathOps

///
///  The following routines:
///
///      initializeHierarchyIntegrator(),
///      initializeHierarchy(),
///      advanceHierarchy(),
///      atRegridPoint(),
///      getIntegratorTime(),
///      getStartTime(),
///      getEndTime(),
///      getIntegratorStep(),
///      getMaxIntegratorSteps(),
///      stepsRemaining(),
///      getPatchHierarchy(),
///      getGriddingAlgorithm(),
///      getHyperbolicLevelIntegrator(),
///      getHyperbolicPatchStrategy()
///
///  allow the AdvDiffHierarchyIntegrator to be used as a
///  hierarchy integrator.
///

void
AdvDiffHierarchyIntegrator::initializeHierarchyIntegrator(
    SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!gridding_alg.isNull());
#endif
    d_gridding_alg = gridding_alg;

    // Setup the tag buffer.
    if (d_using_default_tag_buffer)
    {
        d_tag_buffer.resizeArray(d_gridding_alg->getMaxLevels());
        for (int i = 0; i < d_gridding_alg->getMaxLevels(); ++i)
        {
            d_tag_buffer[i] = d_regrid_interval;
        }
    }
    else
    {
        if (d_tag_buffer.getSize() < d_gridding_alg->getMaxLevels())
        {
            int tsize = d_tag_buffer.getSize();
            d_tag_buffer.resizeArray(d_gridding_alg->getMaxLevels());
            for (int i = tsize; i < d_gridding_alg->getMaxLevels(); ++i)
            {
                d_tag_buffer[i] = d_tag_buffer[tsize-1];
            }
        }
    }

    // Register forcing term data with the hyperbolic level
    // integrator.
    const SAMRAI::hier::IntVector<NDIM> cell_ghosts = CELLG;
    for (unsigned l = 0; l < d_F_vars.size(); ++l)
    {
        // Advected and diffused quantities do not necessarily have
        // source terms.
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var = d_F_vars[l];

        if (!F_var.isNull())
        {
            d_hyp_level_integrator->registerVariable(
                F_var, cell_ghosts,
                SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                d_hierarchy->getGridGeometry(),
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
        }
    }

    // Register temporary data with the variable database.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    d_temp_context = var_db->getContext(d_object_name+"::TEMP_CONTEXT");

    for (unsigned l = 0; l < d_Q_vars.size(); ++l)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
        const int Q_temp_idx = var_db->registerVariableAndContext(
            Q_var, d_temp_context, CELLG);
        (void) Q_temp_idx;

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Psi_var = d_Psi_vars[l];
        const int Psi_temp_idx = var_db->registerVariableAndContext(
            Psi_var, d_temp_context, CELLG);
        (void) Psi_temp_idx;
    }

    // Initialize the SAMRAI::algs::HyperbolicLevelIntegrator.
    //
    // NOTE: This must be done AFTER all variables have been registered.
    d_hyp_level_integrator->initializeLevelIntegrator(d_gridding_alg);

    // Create communications schedules used in setting up the rhs
    // terms and in synchronizing solution data.
    SAMRAI::tbox::Pointer<SAMRAI::geom::CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    for (unsigned l = 0; l < d_Q_vars.size(); ++l)
    {
        ostringstream stream;
        stream << l;
        const std::string& name = stream.str();

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
        const int Q_current_idx = var_db->mapVariableAndContextToIndex(
            Q_var, getCurrentContext());
        const int Q_temp_idx = var_db->mapVariableAndContextToIndex(
            Q_var, d_temp_context);

        SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineOperator<NDIM> > refine_operator =
            grid_geom->lookupRefineOperator(
                Q_var, "CONSTANT_REFINE");
        d_ralgs[name+"::C->T::CONSTANT_REFINE"] = new SAMRAI::xfer::RefineAlgorithm<NDIM>();
        d_ralgs[name+"::C->T::CONSTANT_REFINE"]->
            registerRefine(Q_temp_idx,    // destination
                           Q_current_idx, // source
                           Q_temp_idx,    // temporary work space
                           refine_operator);
        d_rstrategies[name+"::C->T::CONSTANT_REFINE"] = d_bc_op;
    }

    d_calgs["SYNCH_NEW_STATE_DATA"] = new SAMRAI::xfer::CoarsenAlgorithm<NDIM>();
    for (unsigned l = 0; l < d_Q_vars.size(); ++l)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(
            Q_var, getNewContext());

        SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenOperator<NDIM> > coarsen_operator =
            grid_geom->lookupCoarsenOperator(
                Q_var, "CONSERVATIVE_COARSEN");
        d_calgs["SYNCH_NEW_STATE_DATA"]->
            registerCoarsen(Q_new_idx, // destination
                            Q_new_idx, // source
                            coarsen_operator);
    }

    // Setup the Hierarchy math operations object.
    if (d_hier_math_ops.isNull())
    {
        d_hier_math_ops = new STOOLS::HierarchyMathOps(
            d_object_name+"::HierarchyMathOps", d_hierarchy);
        d_is_managing_hier_math_ops = true;
    }

    // One set of operators/specifications is maintained for each
    // variable registered with the integrator.
    for (unsigned l = 0; l < d_Q_vars.size(); ++l)
    {
        ostringstream stream;
        stream << l;
        const std::string& name = stream.str();

        d_helmholtz_specs[l] = new SAMRAI::solv::PoissonSpecifications(
            d_object_name+"::Helmholtz Specs::"+name);
        d_helmholtz_ops[l] = new STOOLS::CCLaplaceOperator(
            d_object_name+"::Helmholtz Operator::"+name,
            *d_helmholtz_specs[l], d_Q_bc_coefs[l]);
    }

    // One set of solvers/preconditioners is maintained for each
    // variable registered with the integrator
    for (unsigned l = 0; l < d_Q_vars.size(); ++l)
    {
        ostringstream stream;
        stream << l;
        const std::string& name = stream.str();

        // Helmholtz solver/operator data.
        if (d_using_FAC)
        {
            d_helmholtz_fac_ops[l] = new STOOLS::CCPoissonFACOperator(
                d_object_name+"::FAC Ops::"+name, d_fac_ops_db);
            d_helmholtz_fac_ops[l]->setPoissonSpecifications(
                *d_helmholtz_specs[l]);
            d_helmholtz_fac_ops[l]->setPhysicalBcCoefs(d_Q_bc_coefs[l]);

            d_helmholtz_fac_pcs[l] = new SAMRAI::solv::FACPreconditioner<NDIM>(
                d_object_name+"::FAC Preconditioner::"+name,
                *d_helmholtz_fac_ops[l], d_fac_pcs_db);
            d_helmholtz_fac_ops[l]->setPreconditioner(
                d_helmholtz_fac_pcs[l]);
        }

        d_helmholtz_solvers[l] = new STOOLS::PETScKrylovLinearSolver(
            d_object_name+"::PETSc Krylov solver::"+name, "adv_diff_");
        d_helmholtz_solvers[l]->setMaxIterations(d_max_iterations);
        d_helmholtz_solvers[l]->setAbsoluteTolerance(d_abs_residual_tol);
        d_helmholtz_solvers[l]->setRelativeTolerance(d_rel_residual_tol);
        d_helmholtz_solvers[l]->setInitialGuessNonzero(true);
        d_helmholtz_solvers[l]->setOperator(d_helmholtz_ops[l]);
        if (d_using_FAC)
        {
            d_helmholtz_solvers[l]->setPreconditioner(
                new STOOLS::FACPreconditionerLSWrapper(
                    d_helmholtz_fac_pcs[l], d_fac_pcs_db));
        }

        // Indicate the the solvers need to be initialized.
        d_helmholtz_solvers_need_init[l] = true;
    }

    // Set the current integration time.
    if (!SAMRAI::tbox::RestartManager::getManager()->isFromRestart())
    {
        d_integrator_time = d_start_time;
        d_integrator_step = 0;
    }

    // Indicate that the integrator has been initialized.
    d_is_initialized = true;

    t_initialize_hierarchy_integrator->stop();
    return;
}// initializeHierarchyIntegrator

double
AdvDiffHierarchyIntegrator::initializeHierarchy()
{
    t_initialize_hierarchy->start();

    if (!d_is_initialized)
    {
        TBOX_ERROR(d_object_name << "::initializeHierarchy()\n" <<
                   "  must initialize the integrator prior to call to initializeHierarchy()." << endl);
    }

    // Initialize the patch hierarchy.
    const bool initial_time = !SAMRAI::tbox::RestartManager::getManager()->isFromRestart();

    if (!initial_time)
    {
        d_hierarchy->getFromRestart(d_gridding_alg->getMaxLevels());
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        d_gridding_alg->getTagAndInitializeStrategy()->
            resetHierarchyConfiguration(d_hierarchy, coarsest_ln, finest_ln);
    }
    else
    {
        d_gridding_alg->makeCoarsestLevel(d_hierarchy,d_start_time);

        int level_number = 0;
        bool done = false;
        while (!done && (d_gridding_alg->levelCanBeRefined(level_number)))
        {
            d_gridding_alg->
                makeFinerLevel(d_hierarchy,
                               d_integrator_time, initial_time,
                               d_tag_buffer[level_number]);

            done = !d_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        // After data on each level is initialized at simulation start
        // time, coarser levels are synchronized with finer levels
        // that didn't exist when the coarser level initial data was
        // set.
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();

        if (finest_ln > 0)
        {
            synchronizeNewLevels(d_hierarchy, coarsest_ln, finest_ln,
                                 d_start_time, initial_time);
        }
    }

    // The timestep is given by the minimum allowable timestep over
    // all levels in the patch hierarchy.
    double dt_next = std::numeric_limits<double>::max();
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt_next = SAMRAI::tbox::Utilities::
            dmin(dt_next, d_hyp_level_integrator->
                 getLevelDt(level, d_integrator_time, initial_time));
    }
    if (d_integrator_time+dt_next > d_end_time)
    {
        dt_next = d_end_time - d_integrator_time;
    }

    t_initialize_hierarchy->stop();
    return dt_next;
}// initializeHierarchy

double
AdvDiffHierarchyIntegrator::advanceHierarchy(
    const double dt)
{
    t_advance_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(d_end_time >= d_integrator_time+dt);
#endif

    // First: regrid (when appropriate).
    const bool do_regrid = ((d_regrid_interval == 0)
                            ? false
                            : (d_integrator_step % d_regrid_interval == 0));
    if (do_regrid) regridHierarchy();

    // Second: integrate the data and synchronize the hierarchy.
    const double dt_next = integrateHierarchy(d_integrator_time, d_integrator_time+dt);
    synchronizeHierarchy();

    // Third: reset all time dependent data.
    resetTimeDependentHierData(d_integrator_time+dt);

    t_advance_hierarchy->stop();
    return dt_next;
}// advanceHierarchy

bool
AdvDiffHierarchyIntegrator::atRegridPoint() const
{
    const int level_number = 0;

    return ((d_integrator_step > 0)
            && d_gridding_alg->levelCanBeRefined(level_number)
            && (d_regrid_interval == 0
                ? false
                : (d_integrator_step % d_regrid_interval == 0)));
}// atRegridPoint

double
AdvDiffHierarchyIntegrator::getIntegratorTime() const
{
    return d_integrator_time;
}// getIntegratorTime

double
AdvDiffHierarchyIntegrator::getStartTime() const
{
    return d_start_time;
}// getStartTime

double
AdvDiffHierarchyIntegrator::getEndTime() const
{
    return d_end_time;
}// getEndTime

int
AdvDiffHierarchyIntegrator::getIntegratorStep() const
{
    return d_integrator_step;
}// getIntegratorStep

int
AdvDiffHierarchyIntegrator::getMaxIntegratorSteps() const
{
    return d_max_integrator_steps;
}// getMaxIntegratorSteps

bool
AdvDiffHierarchyIntegrator::stepsRemaining() const
{
    return (d_integrator_step < d_max_integrator_steps);
}// stepsRemaining

const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> >
AdvDiffHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

SAMRAI::tbox::Pointer<SAMRAI::mesh::GriddingAlgorithm<NDIM> >
AdvDiffHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

SAMRAI::tbox::Pointer<SAMRAI::algs::HyperbolicLevelIntegrator<NDIM> >
AdvDiffHierarchyIntegrator::getHyperbolicLevelIntegrator() const
{
    return d_hyp_level_integrator;
}// getHyperbolicLevelIntegrator

SAMRAI::tbox::Pointer<AdvDiffHypPatchOps>
AdvDiffHierarchyIntegrator::getHyperbolicPatchStrategy() const
{
    return d_hyp_patch_ops;
}// getHyperbolicPatchStrategy

///
///  The following routines:
///
///      regridHierarchy(),
///      integrateHierarchy(),
///      synchronizeHierarchy(),
///      synchronizeNewLevels(),
///      resetTimeDependentHierData(),
///      resetHierDataToPreadvanceState()
///
///  allow the AdvDiffHierarchyIntegrator to provide data
///  management for a time integrator which making use of this class.
///

void
AdvDiffHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    const int coarsest_ln = 0;
    d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln,
                                         d_integrator_time, d_tag_buffer);

    t_regrid_hierarchy->stop();
    return;
}// regridHierarchy

double
AdvDiffHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time)
{
    t_integrate_hierarchy->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(current_time <= new_time);
    assert(d_end_time > d_integrator_time);
    assert(SAMRAI::tbox::Utilities::deq(d_integrator_time,current_time));
#endif

    const double dt = new_time - current_time;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    //////////////////////////////////////////////////////////////////////
    // Predict the advective terms and synchronize them across all
    // levels of the patch hierarchy.
    //////////////////////////////////////////////////////////////////////

    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();

    for (unsigned l = 0; l < d_Q_vars.size(); ++l)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var   = d_Q_vars[l];
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var   = d_F_vars[l];
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Psi_var = d_Psi_vars[l];
        const double mu = d_Q_mus[l];
        const std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& Q_bc_coefs = d_Q_bc_coefs[l];
        SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> F_set = d_F_sets[l];

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > Q_factory =
            Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();

        const int Q_temp_idx = var_db->mapVariableAndContextToIndex(
            Q_var, d_temp_context);
        const int F_current_idx =
            (F_var.isNull()
             ? -1
             : var_db->mapVariableAndContextToIndex(
                 F_var, getCurrentContext()));
        const int Psi_current_idx = var_db->mapVariableAndContextToIndex(
            Psi_var, getCurrentContext());

        // Compute the time-dependent forcing term F(n).
        if (!F_set.isNull() && F_set->isTimeDependent())
        {
            F_set->setDataOnPatchHierarchy(
                F_var, getCurrentContext(),
                d_hierarchy, current_time);
        }

        // Allocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_temp_idx, current_time);
        }

        // Setup the right hand side for the advective flux
        // prediction.
        SAMRAI::solv::PoissonSpecifications mu_spec("mu_spec");
        mu_spec.setCConstant(0.0);
        mu_spec.setDConstant(mu);

        d_bc_op->setPatchDataIndex(Q_temp_idx);
        d_bc_op->setPhysicalBcCoefs(Q_bc_coefs);
        d_bc_op->setHomogeneousBc(false);
        d_bc_op->setExtrapolationType("LINEAR");

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            ostringstream stream;
            stream << l;
            const std::string& name = stream.str();
            d_rscheds[name+"::C->T::CONSTANT_REFINE"][ln]->fillData(current_time);
        }

        for (int depth = 0; depth < Q_depth; ++depth)
        {
            d_hier_math_ops->laplace(
                Psi_current_idx, Psi_var,  // Psi(n)
                mu_spec,                   // Poisson spec
                Q_temp_idx     , Q_var  ,  // Q(n)
                d_rscheds["NONE"],         // don't need to re-fill Q(n) data
                current_time,              // Q(n) bdry fill time
                1.0,                       // gamma
                F_current_idx  , F_var  ,  // F(n)
                depth, depth, depth);      // dst_depth, src1_depth, src2_depth
        }

        // Deallocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_temp_idx);
        }
    }

    double dt_next = std::numeric_limits<double>::max();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);

        const bool first_step = true;
        const bool last_step = false;

        const double level_dt_next = d_hyp_level_integrator->advanceLevel(
            level, d_hierarchy, current_time, new_time,
            first_step, last_step);

        dt_next = SAMRAI::tbox::Utilities::dmin(dt_next,level_dt_next);
    }

    dt_next = SAMRAI::tbox::Utilities::dmin(dt_next,d_grow_dt*dt);

    if (new_time+dt_next >= d_end_time)
    {
        dt_next = d_end_time - new_time;
    }

    if (finest_ln > 0)
    {
        d_hyp_level_integrator->standardLevelSynchronization(
            d_hierarchy, coarsest_ln, finest_ln,
            new_time, current_time);
    }

    //////////////////////////////////////////////////////////////////////
    // Solve one or more Helmholtz equations for Q(n+1).
    //////////////////////////////////////////////////////////////////////

    // Indicate that all solvers need to be reinitialized if the
    // current timestep size is different from the previous one.
    if (!SAMRAI::tbox::Utilities::deq(dt,d_old_dt))
    {
        for (std::map<unsigned,bool>::iterator it = d_helmholtz_solvers_need_init.begin();
             it != d_helmholtz_solvers_need_init.end(); ++it)
        {
            (*it).second = true;
        }
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
    }

    for (unsigned l = 0; l < d_Q_vars.size(); ++l)
    {
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var   = d_Q_vars[l];
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var   = d_F_vars[l];
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Psi_var = d_Psi_vars[l];
        const double mu = d_Q_mus[l];
        const std::vector<const SAMRAI::solv::RobinBcCoefStrategy<NDIM>*>& Q_bc_coefs = d_Q_bc_coefs[l];
        SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> F_set = d_F_sets[l];

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellDataFactory<NDIM,double> > Q_factory =
            Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();

        const int Q_temp_idx = var_db->mapVariableAndContextToIndex(
            Q_var, d_temp_context);
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(
            Q_var, getNewContext());
        const int F_current_idx =
            (F_var.isNull()
             ? -1
             : var_db->mapVariableAndContextToIndex(
                 F_var, getCurrentContext()));
        const int Psi_temp_idx = var_db->mapVariableAndContextToIndex(
            Psi_var, d_temp_context);

        // Compute the time-dependent forcing term F(n+1/2).
        if (!F_set.isNull() && F_set->isTimeDependent())
        {
            F_set->setDataOnPatchHierarchy(
                F_var, getCurrentContext(),
                d_hierarchy, current_time+0.5*dt);
        }

        // Set problem coefficients corresponding to Crank-Nicolson.
        SAMRAI::tbox::Pointer<SAMRAI::solv::PoissonSpecifications> helmholtz_spec = d_helmholtz_specs[l];
        SAMRAI::tbox::Pointer<STOOLS::CCLaplaceOperator> helmholtz_op = d_helmholtz_ops[l];

        helmholtz_spec->setCConstant(1.0);
        helmholtz_spec->setDConstant(-0.5*dt*mu);

        SAMRAI::solv::PoissonSpecifications rhs_spec("rhs_spec");
        rhs_spec.setCConstant(1.0);
        rhs_spec.setDConstant(+0.5*dt*mu);

        // Allocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_temp_idx, current_time);
            level->allocatePatchData(Psi_temp_idx, new_time);
        }

        // Setup the right hand side for the linear solve for Q(n+1).
        d_bc_op->setPatchDataIndex(Q_temp_idx);
        d_bc_op->setPhysicalBcCoefs(Q_bc_coefs);
        d_bc_op->setHomogeneousBc(false);
        d_bc_op->setExtrapolationType("LINEAR");

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            ostringstream stream;
            stream << l;
            const std::string& name = stream.str();
            d_rscheds[name+"::C->T::CONSTANT_REFINE"][ln]->fillData(current_time);
        }

        if (!F_var.isNull())
        {
            d_hier_cc_data_ops->add(Q_new_idx, F_current_idx, Q_new_idx);
        }

        for (int depth = 0; depth < Q_depth; ++depth)
        {
            d_hier_math_ops->laplace(
                Psi_temp_idx, Psi_var,  // Psi(n+1/2)
                rhs_spec,               // Poisson spec
                Q_temp_idx  , Q_var  ,  // Q(n)
                d_rscheds["NONE"],      // don't need to re-fill Q(n) data
                current_time,           // Q(n) bdry fill time
                dt,                     // gamma
                Q_new_idx   , Q_var  ,  // N(n+1/2) = (u*grad Q)(n+1/2)
                depth, depth, depth);   // dst_depth, src1_depth, src2_depth
        }

        // Initialize the linear solver.
        helmholtz_op->setPoissonSpecifications(*helmholtz_spec);
        helmholtz_op->setPhysicalBcCoefs(Q_bc_coefs);
        helmholtz_op->setHomogeneousBc(false);
        helmholtz_op->setTime(new_time);
        helmholtz_op->setHierarchyMathOps(d_hier_math_ops);

        SAMRAI::tbox::Pointer<STOOLS::CCPoissonFACOperator>           helmholtz_fac_op = d_helmholtz_fac_ops[l];
        SAMRAI::tbox::Pointer<SAMRAI::solv::FACPreconditioner<NDIM> > helmholtz_fac_pc = d_helmholtz_fac_pcs[l];
        SAMRAI::tbox::Pointer<STOOLS::KrylovLinearSolver>             helmholtz_solver = d_helmholtz_solvers[l];

        if (d_helmholtz_solvers_need_init[l])
        {
            if (d_do_log) SAMRAI::tbox::plog << d_object_name << ": "
                                             << "Initializing Helmholtz solvers for variable number " << l
                                             << ", dt = " << dt << "\n";
            if (d_using_FAC)
            {
                helmholtz_fac_op->setPoissonSpecifications(*helmholtz_spec);
                helmholtz_fac_op->setTime(new_time);
                helmholtz_fac_op->setResetLevels(d_coarsest_reset_ln, d_finest_reset_ln);
            }
            helmholtz_solver->initializeSolverState(*d_sol_vecs[l],*d_rhs_vecs[l]);

            // Indicate that the solvers do not presently require
            // re-initialization.
            d_helmholtz_solvers_need_init[l] = false;
        }

        // Solve for Q(n+1).
        if (d_using_FAC) helmholtz_fac_op->setTime(new_time);
        helmholtz_solver->solveSystem(*d_sol_vecs[l],*d_rhs_vecs[l]);
        d_hier_cc_data_ops->copyData(Q_new_idx, Q_temp_idx);

        if (d_do_log) SAMRAI::tbox::plog << "AdvDiffHierarchyIntegrator::integrateHierarchy(): linear solve number of iterations = " << helmholtz_solver->getNumIterations() << "\n";
        if (d_do_log) SAMRAI::tbox::plog << "AdvDiffHierarchyIntegrator::integrateHierarchy(): linear solve residual norm        = " << helmholtz_solver->getResidualNorm()  << "\n";

        // Deallocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_temp_idx);
            level->deallocatePatchData(Psi_temp_idx);
        }

    }

    t_integrate_hierarchy->stop();
    return dt_next;
}// integrateHierarchy

void
AdvDiffHierarchyIntegrator::synchronizeHierarchy()
{
    t_synchronize_hierarchy->start();

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = finest_ln; ln > coarsest_ln; --ln)
    {
        d_cscheds["SYNCH_NEW_STATE_DATA"][ln]->coarsenData();
    }

    t_synchronize_hierarchy->stop();
    return;
}// synchronizeHierarchy

void
AdvDiffHierarchyIntegrator::synchronizeNewLevels(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const bool initial_time)
{
    t_synchronize_new_levels->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((coarsest_level >= 0)
           && (coarsest_level < finest_level)
           && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        assert(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    // We use the HyperbolicLevelIntegrator to handle as much data
    // management as possible.
    d_hyp_level_integrator->
        synchronizeNewLevels(hierarchy, coarsest_level, finest_level,
                             sync_time, initial_time);

    t_synchronize_new_levels->stop();
    return;
}// synchronizeNewLevels

void
AdvDiffHierarchyIntegrator::resetTimeDependentHierData(
    const double new_time)
{
    t_reset_time_dependent_hier_data->start();

    // Advance the simulation time.
    d_old_dt = new_time - d_integrator_time;
    d_integrator_time = new_time;
    ++d_integrator_step;

    // Reset the time dependent data.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hyp_level_integrator->resetTimeDependentData(
            d_hierarchy->getPatchLevel(ln),
            d_integrator_time,
            d_gridding_alg->levelCanBeRefined(ln));
    }

    t_reset_time_dependent_hier_data->stop();
    return;
}// resetTimeDependentHierData

void
AdvDiffHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    t_reset_hier_data_to_preadvance_state->start();

    // We use the HyperbolicLevelIntegrator to handle as much data
    // management as possible.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hyp_level_integrator->resetDataToPreadvanceState(d_hierarchy->getPatchLevel(ln));
    }

    t_reset_hier_data_to_preadvance_state->stop();
    return;
}// resetHierDataToPreadvanceState

///
///  The following routines:
///
///      initializeLevelData(),
///      resetHierarchyConfiguration(),
///      applyGradientDetector()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::mesh::StandardTagAndInitStrategy abstract base class.
///

void
AdvDiffHierarchyIntegrator::initializeLevelData(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchLevel<NDIM> > base_old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > old_level = base_old_level;
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
    // We use the HyperbolicLevelIntegrator to handle as much data
    // management as possible.
    d_hyp_level_integrator->
        initializeLevelData(hierarchy, level_number, init_data_time,
                            can_be_refined, initial_time, old_level,
                            allocate_data);

    // Set the initial values of any forcing terms.
    if (initial_time)
    {
        SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
        SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
        for (unsigned l = 0; l < d_F_vars.size(); ++l)
        {
            SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > F_var = d_F_vars[l];

            if (!F_var.isNull())
            {
                const int F_idx = var_db->mapVariableAndContextToIndex(
                    F_var, getCurrentContext());
                SAMRAI::tbox::Pointer<STOOLS::SetDataStrategy> F_set = d_F_sets [l];

                if (!F_set.isNull())
                {
                    F_set->setDataOnPatchLevel(F_idx, F_var, level,
                                               init_data_time, initial_time);
                }
                else
                {
                    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());

                        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,double> > F_data =
                            patch->getPatchData(F_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                        assert(!F_data.isNull());
#endif
                        F_data->fillAll(0.0);
                    }
                }
            }
        }
    }

    t_initialize_level_data->stop();
    return;
}// initializeLevelData

void
AdvDiffHierarchyIntegrator::resetHierarchyConfiguration(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy = base_hierarchy;
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

    // We use the HyperbolicLevelIntegrator to handle as much data
    // management as possible.
    d_hyp_level_integrator->
        resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Reset the Hierarchy data operations for the new hierarchy
    // configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    // Reset the Hierarchy math operations for the new configuration.
    if (d_is_managing_hier_math_ops)
    {
        d_hier_math_ops->setPatchHierarchy(hierarchy);
        d_hier_math_ops->resetLevels(0, finest_hier_level);
    }

    // Get the cell weights data.
    d_wgt_var = d_hier_math_ops->getCellWeightVariable();
    d_wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();

    // Reset the solution and rhs vectors.
    SAMRAI::hier::VariableDatabase<NDIM>* var_db = SAMRAI::hier::VariableDatabase<NDIM>::getDatabase();
    for (unsigned l = 0; l < d_Q_vars.size(); ++l)
    {
        ostringstream stream;
        stream << l;
        const std::string& name = stream.str();

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Q_var = d_Q_vars[l];
        const int Q_temp_idx = var_db->mapVariableAndContextToIndex(
            Q_var, d_temp_context);
        d_sol_vecs[l] = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
            d_object_name+"::sol_vec::"+name, d_hierarchy, 0, finest_hier_level);
        d_sol_vecs[l]->addComponent(Q_var,Q_temp_idx,d_wgt_idx,d_hier_cc_data_ops);

        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellVariable<NDIM,double> > Psi_var = d_Psi_vars[l];
        const int Psi_temp_idx = var_db->mapVariableAndContextToIndex(
            Psi_var, d_temp_context);
        d_rhs_vecs[l] = new SAMRAI::solv::SAMRAIVectorReal<NDIM,double>(
            d_object_name+"::rhs_vec::"+name, d_hierarchy, 0, finest_hier_level);
        d_rhs_vecs[l]->addComponent(Psi_var,Psi_temp_idx,d_wgt_idx,d_hier_cc_data_ops);
    }

    // Indicate that all linear solvers must be re-initialized.
    for (std::map<unsigned,bool>::iterator it = d_helmholtz_solvers_need_init.begin();
         it != d_helmholtz_solvers_need_init.end(); ++it)
    {
        (*it).second = true;
    }

    d_coarsest_reset_ln = coarsest_level;
    d_finest_reset_ln = finest_level;

    // If we have added or removed a level, resize the schedule
    // vectors.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it != d_ralgs.end(); ++it)
    {
        d_rscheds[(*it).first].resize(finest_hier_level+1);
    }

    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it != d_calgs.end(); ++it)
    {
        d_cscheds[(*it).first].resize(finest_hier_level+1);
    }

    // (Re)build refine communication schedules.  These are created
    // for all levels in the hierarchy.
    for (RefineAlgMap::const_iterator it = d_ralgs.begin();
         it != d_ralgs.end(); ++it)
    {
        for (int ln = coarsest_level; ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);

            d_rscheds[(*it).first][ln] = (*it).second->
                createSchedule(
                    level, ln-1, hierarchy, d_rstrategies[(*it).first]);
        }
    }

    // (Re)build coarsen communication schedules.  These are set only
    // for levels >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin();
         it != d_calgs.end(); ++it)
    {
        for (int ln = SAMRAI::tbox::Utilities::imax(coarsest_level,1);
             ln <= finest_hier_level; ++ln)
        {
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > coarser_level =
                hierarchy->getPatchLevel(ln-1);

            d_cscheds[(*it).first][ln] = (*it).second->
                createSchedule(
                    coarser_level, level, d_cstrategies[(*it).first]);
        }
    }

    // Reset the "empty" schedule vectors.
    d_rscheds["NONE"] = std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::RefineSchedule<NDIM> > >(finest_hier_level+1);
    d_cscheds["NONE"] = std::vector<SAMRAI::tbox::Pointer<SAMRAI::xfer::CoarsenSchedule<NDIM> > >(finest_hier_level+1);

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
AdvDiffHierarchyIntegrator::applyGradientDetector(
    const SAMRAI::tbox::Pointer<SAMRAI::hier::BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    t_apply_gradient_detector->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!hierarchy.isNull());
    assert((level_number >= 0)
           && (level_number <= hierarchy->getFinestLevelNumber()));
    assert(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Due to tag buffers, it is necessary to untag all cells prior to
    // tagging.
    for (SAMRAI::hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM> > patch = level->getPatch(p());
        SAMRAI::tbox::Pointer<SAMRAI::pdat::CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
        tags_data->fillAll(0);
    }

    // Tag cells for refinement according to the criteria specified by
    // the criteria specified by the SAMRAI::algs::HyperbolicLevelIntegrator<NDIM>.
    d_hyp_level_integrator->
        applyGradientDetector(hierarchy, level_number, error_data_time,
                              tag_index, initial_time,
                              uses_richardson_extrapolation_too);

    t_apply_gradient_detector->stop();
    return;
}// applyGradientDetector

///
///  The following routines:
///
///      getCurrentContext(),
///      getNewContext(),
///      getOldContext(),
///      getScratchContext(),
///      getPlotContext()
///
///  allow access to the various variable contexts maintained by the
///  integrator.
///

///
/// We simply reuse the SAMRAI::hier::VariableContext objects defined
/// in the SAMRAI::algs::HyperbolicLevelIntegrator<NDIM> object.
///

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
AdvDiffHierarchyIntegrator::getCurrentContext() const
{
    return d_hyp_level_integrator->getCurrentContext();
}// getCurrentContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
AdvDiffHierarchyIntegrator::getNewContext() const
{
    return d_hyp_level_integrator->getNewContext();
}// getNewContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
AdvDiffHierarchyIntegrator::getOldContext() const
{
    return d_hyp_level_integrator->getOldContext();
}// getOldContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
AdvDiffHierarchyIntegrator::getScratchContext() const
{
    return d_hyp_level_integrator->getScratchContext();
}// getScratchContext

SAMRAI::tbox::Pointer<SAMRAI::hier::VariableContext>
AdvDiffHierarchyIntegrator::getPlotContext() const
{
    return d_hyp_level_integrator->getPlotContext();
}// getPlotContext

///
///  The following routines:
///
///      putToDatabase()
///
///  are concrete implementations of functions declared in the
///  SAMRAI::tbox::Serializable abstract base class.
///

void
AdvDiffHierarchyIntegrator::putToDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif

    db->putInteger("ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION",
                   ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION);
    db->putDouble("d_start_time", d_start_time);
    db->putDouble("d_end_time", d_end_time);
    db->putDouble("d_grow_dt", d_grow_dt);
    db->putInteger("d_max_integrator_steps", d_max_integrator_steps);
    db->putInteger("d_regrid_interval", d_regrid_interval);
    db->putBool("d_using_default_tag_buffer", d_using_default_tag_buffer);
    db->putIntegerArray("d_tag_buffer", d_tag_buffer);
    db->putDouble("d_integrator_time", d_integrator_time);
    db->putInteger("d_integrator_step", d_integrator_step);

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
AdvDiffHierarchyIntegrator::printClassData(
    std::ostream& os) const
{
    os << "\nAdvDiffHierarchyIntegrator::printClassData..." << endl;
    os << "\nAdvDiffHierarchyIntegrator: this = " << const_cast<AdvDiffHierarchyIntegrator*>(this) << endl;
    os << "d_object_name = " << d_object_name << endl;
    os << "d_integrator_time = " << d_integrator_time << "\n"
       << "d_start_time = " << d_start_time << "\n"
       << "d_end_time = " << d_end_time << "\n"
       << "d_integrator_step = " << d_integrator_step << "\n"
       << "d_grow_dt = " << d_grow_dt << endl;
    os << "I AM INCOMPLETE!!!!!!!!!" << endl;
    return;
}// printClassData

/////////////////////////////// PRIVATE //////////////////////////////////////

void
AdvDiffHierarchyIntegrator::getFromInput(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    assert(!db.isNull());
#endif
    // Read in data members from input database.
    d_end_time = db->getDoubleWithDefault("end_time", d_end_time);
    d_grow_dt = db->getDoubleWithDefault("grow_dt", d_grow_dt);

    d_max_integrator_steps = db->getIntegerWithDefault(
        "max_integrator_steps", d_max_integrator_steps);

    d_regrid_interval = db->getIntegerWithDefault(
        "regrid_interval", d_regrid_interval);

    if (db->keyExists("tag_buffer"))
    {
        d_tag_buffer = db->getIntegerArray("tag_buffer");
        d_using_default_tag_buffer = false;
    }
    else
    {
        d_using_default_tag_buffer = true;
        TBOX_WARNING(d_object_name << ":  "
                     << "Key data `tag_buffer' not found in input.  "
                     << "Default values used.  See class header for details.");
    }

    d_do_log = db->getBoolWithDefault("enable_logging", d_do_log);

    d_max_iterations = db->getIntegerWithDefault("max_iterations", d_max_iterations);
    d_abs_residual_tol = db->getDoubleWithDefault("abs_residual_tol", d_abs_residual_tol);
    d_rel_residual_tol = db->getDoubleWithDefault("rel_residual_tol", d_rel_residual_tol);
    d_using_FAC = db->getBoolWithDefault("using_FAC", d_using_FAC);

    if (!is_from_restart)
    {
        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);
    }

    return;
}// getFromInput

void
AdvDiffHierarchyIntegrator::getFromRestart()
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
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart database corresponding to "
                   << d_object_name << " not found in restart file.");
    }

    int ver = db->getInteger("ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  "
                   << "Restart file version different than class version.");
    }

    d_start_time = db->getDouble("d_start_time");
    d_end_time = db->getDouble("d_end_time");
    d_grow_dt = db->getDouble("d_grow_dt");
    d_max_integrator_steps = db->getInteger("d_max_integrator_steps");
    d_regrid_interval = db->getInteger("d_regrid_interval");
    d_using_default_tag_buffer = db->getBool("d_using_default_tag_buffer");
    d_tag_buffer = db->getIntegerArray("d_tag_buffer");
    d_integrator_time = db->getDouble("d_integrator_time");
    d_integrator_step = db->getInteger("d_integrator_step");

    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

/////////////////////// TEMPLATE INSTANTIATION ///////////////////////////////

#include <tbox/Pointer.C>
template class SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
