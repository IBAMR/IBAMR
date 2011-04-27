// Filename: AdvDiffHierarchyIntegrator.C
// Created on 17 Mar 2004 by Boyce Griffith
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

// IBAMR INCLUDES
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/PETScKrylovLinearSolver.h>

// SAMRAI INCLUDES
#include <CartesianGridGeometry.h>
#include <CellDataFactory.h>
#include <CoarsenOperator.h>
#include <HierarchyDataOpsManager.h>
#include <VariableDatabase.h>
#include <tbox/MathUtilities.h>
#include <tbox/NullDatabase.h>
#include <tbox/RestartManager.h>
#include <tbox/Timer.h>
#include <tbox/TimerManager.h>
#include <tbox/Utilities.h>

// C++ STDLIB INCLUDES
#include <iterator>
#include <limits>
#include <set>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Timers.
static Pointer<Timer> t_initialize_hierarchy_integrator;
static Pointer<Timer> t_initialize_hierarchy;
static Pointer<Timer> t_advance_hierarchy;
static Pointer<Timer> t_regrid_hierarchy;
static Pointer<Timer> t_integrate_hierarchy;
static Pointer<Timer> t_synchronize_hierarchy;
static Pointer<Timer> t_synchronize_new_levels;
static Pointer<Timer> t_reset_time_dependent_hier_data;
static Pointer<Timer> t_reset_hier_data_to_preadvance_state;
static Pointer<Timer> t_initialize_level_data;
static Pointer<Timer> t_reset_hierarchy_configuration;
static Pointer<Timer> t_apply_gradient_detector;
static Pointer<Timer> t_put_to_database;

// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Version of AdvDiffHierarchyIntegrator restart file data.
static const int ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffHierarchyIntegrator::AdvDiffHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GodunovAdvector> explicit_predictor,
    bool register_for_restart)
    : d_viscous_timestepping_type(CRANK_NICOLSON),
      d_u_var(),
      d_manage_u_data(),
      d_u_is_div_free(),
      d_u_fcn(),
      d_F_var(),
      d_manage_F_data(),
      d_F_fcn(),
      d_grad_Phi_var(),
      d_manage_grad_Phi_data(),
      d_grad_Phi_fcn(),
      d_Q_var(),
      d_Psi_var(),
      d_manage_Q_data(),
      d_Q_u_map(),
      d_Q_grad_Phi_map(),
      d_Q_F_map(),
      d_Q_Psi_map(),
      d_Q_difference_form(),
      d_Q_diffusion_coef(),
      d_Q_damping_coef(),
      d_Q_init(),
      d_Q_bc_coef(),
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
      d_regrid_mode(STANDARD),
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
      d_calgs(),
      d_cstrategies(),
      d_cscheds(),
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
      d_fac_op_db(NULL),
      d_fac_pc_db(NULL)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT(!explicit_predictor.isNull());
#endif

    if (d_registered_for_restart)
    {
        RestartManager::getManager()->registerRestartItem(d_object_name, this);
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Initialize the HyperbolicPatchStrategy and
    // HyperbolicLevelIntegrator objects used to provide the
    // numerical routines for explicitly integrating the advective terms.
    Pointer<Database> hyp_patch_ops_input_db;
    if (input_db->keyExists("AdvDiffHypPatchOps"))
    {
        hyp_patch_ops_input_db = input_db->getDatabase("AdvDiffHypPatchOps");
    }
    else
    {
        hyp_patch_ops_input_db = new NullDatabase();
    }

    d_hyp_patch_ops = new AdvDiffHypPatchOps(
        object_name+"::AdvDiffHypPatchOps",
        hyp_patch_ops_input_db,
        explicit_predictor,
        d_hierarchy->getGridGeometry(),
        d_registered_for_restart);

    Pointer<Database> hyp_level_integrator_input_db;
    if (input_db->keyExists("HyperbolicLevelIntegrator"))
    {
        hyp_level_integrator_input_db = input_db->getDatabase("HyperbolicLevelIntegrator");
    }
    else
    {
        hyp_level_integrator_input_db = new NullDatabase();
    }

    static const bool using_time_refinement = false;
    d_hyp_level_integrator = new HyperbolicLevelIntegrator<NDIM>(
        object_name+"::HyperbolicLevelIntegrator",
        hyp_level_integrator_input_db,
        d_hyp_patch_ops,
        register_for_restart,
        using_time_refinement);

    // Get initialization data for the FAC ops and FAC preconditioners.
    if (d_using_FAC)
    {
        if (input_db->keyExists("FACOp"))
        {
            d_fac_op_db = input_db->getDatabase("FACOp");
        }
        else if (input_db->keyExists("FACOps"))
        {
            d_fac_op_db = input_db->getDatabase("FACOps");
        }

        if (input_db->keyExists("FACPreconditioner"))
        {
            d_fac_pc_db = input_db->getDatabase("FACPreconditioner");
        }
        else if (input_db->keyExists("FACPreconditioners"))
        {
            d_fac_pc_db = input_db->getDatabase("FACPreconditioners");
        }
    }

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<CellVariable<NDIM,double> > cc_var = new CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, hierarchy, true);

    // Setup Timers.
    IBAMR_DO_ONCE(
        t_initialize_hierarchy_integrator     = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::initializeHierarchyIntegrator()");
        t_initialize_hierarchy                = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::initializeHierarchy()");
        t_advance_hierarchy                   = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::advanceHierarchy()");
        t_regrid_hierarchy                    = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::regridHierarchy()");
        t_integrate_hierarchy                 = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::integrateHierarchy()");
        t_synchronize_hierarchy               = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::synchronizeHierarchy()");
        t_synchronize_new_levels              = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::synchronizeNewLevels()");
        t_reset_time_dependent_hier_data      = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::resetTimeDependentHierData()");
        t_reset_hier_data_to_preadvance_state = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::resetHierDataToPreadvanceState()");
        t_initialize_level_data               = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::initializeLevelData()");
        t_reset_hierarchy_configuration       = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::resetHierarchyConfiguration()");
        t_apply_gradient_detector             = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::applyGradientDetector()");
        t_put_to_database                     = TimerManager::getManager()->getTimer("IBAMR::AdvDiffHierarchyIntegrator::putToDatabase()");
                  );
    return;
}// AdvDiffHierarchyIntegrator

AdvDiffHierarchyIntegrator::~AdvDiffHierarchyIntegrator()
{
    if (d_registered_for_restart)
    {
        RestartManager::getManager()->unregisterRestartItem(d_object_name);
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

const ViscousTimesteppingType&
AdvDiffHierarchyIntegrator::getViscousTimesteppingType() const
{
    return d_viscous_timestepping_type;
}// getViscousTimesteppingType

void
AdvDiffHierarchyIntegrator::registerVisItDataWriter(
    Pointer<VisItDataWriter<NDIM> > visit_writer)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!visit_writer.isNull());
#endif
    d_hyp_patch_ops->registerVisItDataWriter(visit_writer);
    return;
}// registerVisItDataWriter

void
AdvDiffHierarchyIntegrator::registerAdvectionVelocity(
    Pointer<FaceVariable<NDIM,double> > u_var,
    const bool manage_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_var.isNull());
#endif
    d_u_var.insert(u_var);
    d_manage_u_data[u_var] = manage_data;
    d_u_is_div_free[u_var] = true;
    return;
}// registerAdvectionVelocity

void
AdvDiffHierarchyIntegrator::setAdvectionVelocityIsDivergenceFree(
    Pointer<FaceVariable<NDIM,double> > u_var,
    const bool is_div_free)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
#endif
    d_u_is_div_free[u_var] = is_div_free;
    return;
}// setAdvectionVelocityIsDivergenceFree

void
AdvDiffHierarchyIntegrator::setAdvectionVelocityFunction(
    Pointer<FaceVariable<NDIM,double> > u_var,
    Pointer<IBTK::CartGridFunction> u_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
#endif
    d_u_fcn[u_var] = u_fcn;
    return;
}// setAdvectionVelocityFunction

void
AdvDiffHierarchyIntegrator::registerSourceTerm(
    Pointer<CellVariable<NDIM,double> > F_var,
    const bool manage_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!F_var.isNull());
#endif
    d_F_var.insert(F_var);
    d_manage_F_data[F_var] = manage_data;
    return;
}// registerSourceTerm

void
AdvDiffHierarchyIntegrator::setSourceTermFunction(
    Pointer<CellVariable<NDIM,double> > F_var,
    Pointer<IBTK::CartGridFunction> F_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_F_var.find(F_var) != d_F_var.end());
#endif
    d_F_fcn[F_var] = F_fcn;
    return;
}// setSourceTermFunction

void
AdvDiffHierarchyIntegrator::registerIncompressibilityFixTerm(
    Pointer<FaceVariable<NDIM,double> > grad_Phi_var,
    const bool manage_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!grad_Phi_var.isNull());
#endif
    d_grad_Phi_var.insert(grad_Phi_var);
    d_manage_grad_Phi_data[grad_Phi_var] = manage_data;
    return;
}// registerIncompressibilityFixTerm

void
AdvDiffHierarchyIntegrator::setIncompressibilityFixTermFunction(
    Pointer<FaceVariable<NDIM,double> > grad_Phi_var,
    Pointer<IBTK::CartGridFunction> grad_Phi_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_grad_Phi_var.find(grad_Phi_var) != d_grad_Phi_var.end());
#endif
    d_grad_Phi_fcn[grad_Phi_var] = grad_Phi_fcn;
    return;
}// setIncompressibilityFixTermFunction

void
AdvDiffHierarchyIntegrator::registerTransportedQuantity(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const bool manage_data)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Q_var.isNull());
#endif
    d_Q_var.insert(Q_var);
    d_manage_Q_data[Q_var] = manage_data;
    d_Q_difference_form[Q_var] = CONSERVATIVE;
    d_Q_diffusion_coef[Q_var] = 0.0;
    d_Q_damping_coef[Q_var] = 0.0;
    return;
}// registerTransportedQuantity

void
AdvDiffHierarchyIntegrator::setAdvectionVelocity(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Pointer<FaceVariable<NDIM,double> > u_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    TBOX_ASSERT(d_u_var.find(u_var) != d_u_var.end());
#endif
    d_Q_u_map[Q_var] = u_var;
    return;
}// setAdvectionVelocity

void
AdvDiffHierarchyIntegrator::setSourceTerm(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Pointer<CellVariable<NDIM,double> > F_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    TBOX_ASSERT(d_F_var.find(F_var) != d_F_var.end());
#endif
    d_Q_F_map[Q_var] = F_var;
    return;
}// setSourceTerm

void
AdvDiffHierarchyIntegrator::setIncompressibilityFixTerm(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Pointer<FaceVariable<NDIM,double> > grad_Phi_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    TBOX_ASSERT(d_grad_Phi_var.find(grad_Phi_var) != d_grad_Phi_var.end());
#endif
    d_Q_grad_Phi_map[Q_var] = grad_Phi_var;
    return;
}// setIncompressibilityFixTerm

void
AdvDiffHierarchyIntegrator::setConvectiveDifferencingType(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const ConvectiveDifferencingType& difference_form)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_difference_form[Q_var] = difference_form;
    return;
}// setConvectiveDifferencingType

void
AdvDiffHierarchyIntegrator::setDiffusionCoefficient(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const double& mu)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_diffusion_coef[Q_var] = mu;
    return;
}// setDiffusionCoefficient

void
AdvDiffHierarchyIntegrator::setDampingCoefficient(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const double& lambda)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_damping_coef[Q_var] = lambda;
    return;
}// setDampingCoefficient

void
AdvDiffHierarchyIntegrator::setInitialConditions(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Pointer<IBTK::CartGridFunction> Q_init)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_init[Q_var] = Q_init;
    return;
}// setInitialConditions

void
AdvDiffHierarchyIntegrator::setBoundaryConditions(
    Pointer<CellVariable<NDIM,double> > Q_var,
    RobinBcCoefStrategy<NDIM>* Q_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
    const unsigned Q_depth = Q_factory->getDefaultDepth();
    TBOX_ASSERT(Q_depth == 1);
#endif
    d_Q_bc_coef[Q_var] = std::vector<RobinBcCoefStrategy<NDIM>*>(1,Q_bc_coef);
    return;
}// setBoundaryConditions

void
AdvDiffHierarchyIntegrator::setBoundaryConditions(
    Pointer<CellVariable<NDIM,double> > Q_var,
    std::vector<RobinBcCoefStrategy<NDIM>*> Q_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
    const unsigned Q_depth = Q_factory->getDefaultDepth();
    TBOX_ASSERT(Q_depth == Q_bc_coef.size());
#endif
    d_Q_bc_coef[Q_var] = Q_bc_coef;
    return;
}// setBoundaryConditions

///
///  The following routines:
///
///      getHierarchyMathOps(),
///      setHierarchyMathOps(),
///      isManagingHierarchyMathOps()
///
///  allow for the sharing of a single HierarchyMathOps object between multiple
///  HierarchyIntegrator objects.
///

Pointer<HierarchyMathOps>
AdvDiffHierarchyIntegrator::getHierarchyMathOps() const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_hier_math_ops.isNull());
#endif
    return d_hier_math_ops;
}// getHierarchyMathOps

void
AdvDiffHierarchyIntegrator::setHierarchyMathOps(
    Pointer<HierarchyMathOps> hier_math_ops,
    const bool manage_ops)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hier_math_ops.isNull());
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
///  allow the AdvDiffHierarchyIntegrator to be used as a hierarchy integrator.
///

void
AdvDiffHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    t_initialize_hierarchy_integrator->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!gridding_alg.isNull());
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

    // Setup an extra variable context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_temp_context = var_db->getContext(d_object_name+"::TEMP_CONTEXT");

    // Register variables with the hyperbolic level integrator.
    TBOX_ASSERT(false);
#if 0
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin();
         cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        d_hyp_patch_ops->registerAdvectionVelocity(u_var,d_manage_u_data[u_var]);
        d_hyp_patch_ops->setAdvectionVelocityIsDivergenceFree(d_u_is_div_free[u_var]);
        if (!d_u_fcn [u_var].isNull()) d_hyp_patch_ops->setAdvectionVelocityFunction(d_u_fcn[u_var]);
    }

    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_grad_Phi_var.begin();
         cit != d_grad_Phi_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > grad_Phi_var = *cit;
        d_hyp_patch_ops->registerIncompressibilityFix(grad_Phi_var,d_manage_grad_Phi_data[grad_Phi_var]);
        if (!d_grad_Phi_fcn [grad_Phi_var].isNull()) d_hyp_patch_ops->setIncompressibilityFixFunction(d_grad_Phi_fcn[grad_Phi_var]);
    }

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_F_var.begin();
         cit != d_F_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > F_var = *cit;
        if (d_manage_F_data[F_var])
        {
            d_hyp_level_integrator->registerVariable(
                F_var, cell_ghosts,
                HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
                d_hierarchy->getGridGeometry(),
                "CONSERVATIVE_COARSEN",
                "CONSERVATIVE_LINEAR_REFINE");
        }
    }

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Psi_var.begin();
         cit != d_Psi_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Psi_var = *cit;
        d_hyp_patch_ops->registerSourceTerm(Psi_var, true);
    }

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        d_hyp_patch_ops->registerTransportedQuantity(Q_var,d_manage_Q_data[Q_var]);
        d_hyp_patch_ops->setAdvectionVelocity(Q_var,d_Q_u_map[Q_var]);
        d_hyp_patch_ops->setSourceTerm(Q_var,d_Q_Psi_map[Q_var]);
        if (!d_Q_grad_Phi_map[Q_var].isNull()) d_hyp_patch_ops->setIncompressibilityFixTerm(d_Q_grad_Phi_map[Q_var]);
        d_hyp_patch_ops->setConvectiveDifferencingType(d_Q_difference_form[Q_var]);
        if (!d_Q_init[Q_var].isNull()) d_hyp_patch_ops->setAdvectionVelocityInitialConditions(d_Q_init[Q_var]);
        if (!d_Q_bc_coef[Q_var].empty()) d_hyp_patch_ops->setAdvectionVelocityBoundaryConditions(d_Q_bc_coef[Q_var]);
    }
#endif

    // Initialize the HyperbolicLevelIntegrator.
    //
    // NOTE: This must be done AFTER all variables have been registered.
    d_hyp_level_integrator->initializeLevelIntegrator(d_gridding_alg);

    // Create coarsening communications algorithms, used in synchronizing
    // refined regions of coarse data with the underlying fine data.
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();
    d_calgs["SYNCH_NEW_STATE_DATA"] = new CoarsenAlgorithm<NDIM>();
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        Pointer<CoarsenOperator<NDIM> > coarsen_operator = grid_geom->lookupCoarsenOperator(Q_var, "CONSERVATIVE_COARSEN");
        d_calgs["SYNCH_NEW_STATE_DATA"]->registerCoarsen(Q_new_idx, Q_new_idx, coarsen_operator);
    }

    // Setup the Hierarchy math operations object.
    if (d_hier_math_ops.isNull())
    {
        d_hier_math_ops = new HierarchyMathOps(d_object_name+"::HierarchyMathOps", d_hierarchy);
        d_is_managing_hier_math_ops = true;
    }

    // Operators and solvers are maintained for each variable registered with the
    // integrator.
    d_helmholtz_specs.reserve(d_Q_var.size());
    d_helmholtz_ops.resize(d_Q_var.size());
    d_helmholtz_fac_ops.resize(d_Q_var.size());
    d_helmholtz_fac_pcs.resize(d_Q_var.size());
    d_helmholtz_solvers.resize(d_Q_var.size());
    d_helmholtz_solvers_need_init.resize(d_Q_var.size());
    unsigned l = 0;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        const std::string& name = Q_var->getName();

        d_helmholtz_specs.push_back(PoissonSpecifications(d_object_name+"::Helmholtz Specs::"+name));
        d_helmholtz_ops[l] = new CCLaplaceOperator(d_object_name+"::Helmholtz Operator::"+name,
                                                   d_helmholtz_specs[l], d_Q_bc_coef[Q_var]);

        if (d_using_FAC)
        {
            d_helmholtz_fac_ops[l] = new CCPoissonFACOperator(d_object_name+"::FAC Ops::"+name, d_fac_op_db);
            d_helmholtz_fac_ops[l]->setPoissonSpecifications(d_helmholtz_specs[l]);
            d_helmholtz_fac_ops[l]->setPhysicalBcCoefs(d_Q_bc_coef[Q_var]);
            d_helmholtz_fac_pcs[l] = new FACPreconditioner(d_object_name+"::FAC Preconditioner::"+name,
                                                           *d_helmholtz_fac_ops[l], d_fac_pc_db);
        }

        d_helmholtz_solvers[l] = new PETScKrylovLinearSolver(d_object_name+"::PETSc Krylov solver::"+name, "adv_diff_");
        d_helmholtz_solvers[l]->setMaxIterations(d_max_iterations);
        d_helmholtz_solvers[l]->setAbsoluteTolerance(d_abs_residual_tol);
        d_helmholtz_solvers[l]->setRelativeTolerance(d_rel_residual_tol);
        d_helmholtz_solvers[l]->setInitialGuessNonzero(true);
        d_helmholtz_solvers[l]->setOperator(d_helmholtz_ops[l]);
        if (d_using_FAC)
        {
            d_helmholtz_solvers[l]->setPreconditioner(d_helmholtz_fac_pcs[l]);
        }

        // Indicate that the solvers need to be initialized.
        d_helmholtz_solvers_need_init[l] = true;
    }

    // Set the current integration time.
    if (!RestartManager::getManager()->isFromRestart())
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
                   "  must initialize the integrator prior to call to initializeHierarchy()." << std::endl);
    }

    // Initialize the patch hierarchy.
    const bool initial_time = !RestartManager::getManager()->isFromRestart();

    if (!initial_time)
    {
        d_hierarchy->getFromRestart(d_gridding_alg->getMaxLevels());
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        d_gridding_alg->getTagAndInitializeStrategy()->resetHierarchyConfiguration(d_hierarchy, coarsest_ln, finest_ln);
    }
    else
    {
        d_gridding_alg->makeCoarsestLevel(d_hierarchy,d_start_time);
        int level_number = 0;
        bool done = false;
        while (!done && (d_gridding_alg->levelCanBeRefined(level_number)))
        {
            d_gridding_alg->makeFinerLevel(d_hierarchy, d_integrator_time, initial_time, d_tag_buffer[level_number]);
            done = !d_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        // After data on each level is initialized at simulation start time,
        // coarser levels are synchronized with finer levels that didn't exist
        // when the coarser level initial data was set.
        const int coarsest_ln = 0;
        const int finest_ln = d_hierarchy->getFinestLevelNumber();
        if (finest_ln > 0)
        {
            synchronizeNewLevels(d_hierarchy, coarsest_ln, finest_ln, d_start_time, initial_time);
        }
    }

    // The timestep is given by the minimum allowable timestep over all levels
    // in the patch hierarchy.
    double dt_next = std::numeric_limits<double>::max();
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt_next = std::min(dt_next, d_hyp_level_integrator->getLevelDt(level, d_integrator_time, initial_time));
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
    TBOX_ASSERT(d_end_time >= d_integrator_time+dt);
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

const Pointer<PatchHierarchy<NDIM> >
AdvDiffHierarchyIntegrator::getPatchHierarchy() const
{
    return d_hierarchy;
}// getPatchHierarchy

Pointer<GriddingAlgorithm<NDIM> >
AdvDiffHierarchyIntegrator::getGriddingAlgorithm() const
{
    return d_gridding_alg;
}// getGriddingAlgorithm

Pointer<HyperbolicLevelIntegrator<NDIM> >
AdvDiffHierarchyIntegrator::getHyperbolicLevelIntegrator() const
{
    return d_hyp_level_integrator;
}// getHyperbolicLevelIntegrator

Pointer<AdvDiffHypPatchOps>
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
///  allow the AdvDiffHierarchyIntegrator to provide data management for a time
///  integrator which making use of this class.
///

void
AdvDiffHierarchyIntegrator::regridHierarchy()
{
    t_regrid_hierarchy->start();

    const int coarsest_ln = 0;
    switch (d_regrid_mode)
    {
        case STANDARD:
            d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln, d_integrator_time, d_tag_buffer);
            break;
        case AGGRESSIVE:
            for (int k = 0; k < std::max(1,d_hierarchy->getFinestLevelNumber()); ++k)
            {
                d_gridding_alg->regridAllFinerLevels(d_hierarchy, coarsest_ln, d_integrator_time, d_tag_buffer);
            }
            break;
        default:
            TBOX_ERROR(d_object_name << "::regridHierarchy():\n"
                       << "  unrecognized regrid mode: " << enum_to_string<RegridMode>(d_regrid_mode) << "." << std::endl);
    }

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
    TBOX_ASSERT(current_time <= new_time);
    TBOX_ASSERT(d_end_time > d_integrator_time);
    TBOX_ASSERT(MathUtilities<double>::equalEps(d_integrator_time,current_time));
#endif

    const double dt = new_time - current_time;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    ////////////////////////////////////////////////////////////////////////////
    // Compute any time-dependent source terms at time-level n.
    ////////////////////////////////////////////////////////////////////////////

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_F_var.begin();
         cit != d_F_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > F_var = *cit;
        Pointer<CartGridFunction> F_fcn = d_F_fcn[F_var];
        if (!F_fcn.isNull() && F_fcn->isTimeDependent())
        {
            const int F_current_idx = var_db->mapVariableAndContextToIndex(F_var, getCurrentContext());
            F_fcn->setDataOnPatchHierarchy(F_current_idx, F_var, d_hierarchy, current_time);
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Predict the advective terms and synchronize them across all levels of the
    // patch hierarchy.
    ////////////////////////////////////////////////////////////////////////////

    unsigned l = 0;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > Q_var   = *cit;
        Pointer<CellVariable<NDIM,double> > F_var   = d_Q_F_map[Q_var];
        Pointer<CellVariable<NDIM,double> > Psi_var = d_Q_Psi_map[Q_var];
        const double mu = d_Q_diffusion_coef[Q_var];
        const double lambda = d_Q_damping_coef[Q_var];

        Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();

        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_temp_idx = var_db->mapVariableAndContextToIndex(Q_var, d_temp_context);
        const int F_current_idx = (F_var.isNull() ? -1 : var_db->mapVariableAndContextToIndex(F_var, getCurrentContext()));
        const int Psi_current_idx = var_db->mapVariableAndContextToIndex(Psi_var, getCurrentContext());

        // Allocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_temp_idx, current_time);
        }

        // Setup the right hand side for the advective flux prediction.
        d_hier_cc_data_ops->copyData(Q_temp_idx, Q_current_idx, false);
        d_hier_bdry_fill_ops[l]->setHomogeneousBc(false);
        d_hier_bdry_fill_ops[l]->fillData(current_time);

        PoissonSpecifications mu_spec("mu_spec");
        mu_spec.setCConstant(-lambda);
        mu_spec.setDConstant(+mu    );

        for (int depth = 0; depth < Q_depth; ++depth)
        {
            d_hier_math_ops->laplace(
                Psi_current_idx, Psi_var,  // Psi(n)
                mu_spec,                   // Poisson spec
                Q_temp_idx     , Q_var  ,  // Q(n)
                d_no_fill_op,              // don't need to re-fill Q(n) data
                current_time,              // Q(n) bdry fill time
                1.0,                       // gamma
                F_current_idx  , F_var  ,  // F(n)
                depth, depth, depth);      // dst_depth, src1_depth, src2_depth
        }

        // Deallocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_temp_idx);
        }
    }

    double dt_next = std::numeric_limits<double>::max();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        const bool first_step = true;
        const bool last_step = false;
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        const double level_dt_next = d_hyp_level_integrator->advanceLevel(
            level, d_hierarchy, current_time, new_time, first_step, last_step);
        dt_next = std::min(dt_next,level_dt_next);
    }

    dt_next = std::min(dt_next,d_grow_dt*dt);

    if (new_time+dt_next >= d_end_time)
    {
        dt_next = d_end_time - new_time;
    }

    if (finest_ln > 0)
    {
        d_hyp_level_integrator->standardLevelSynchronization(
            d_hierarchy, coarsest_ln, finest_ln, new_time, current_time);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Compute any time-dependent source terms at time-level n+1/2.
    ////////////////////////////////////////////////////////////////////////////

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_F_var.begin();
         cit != d_F_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > F_var = *cit;
        Pointer<CartGridFunction> F_fcn = d_F_fcn[F_var];
        if (!F_fcn.isNull() && F_fcn->isTimeDependent())
        {
            const int F_current_idx = var_db->mapVariableAndContextToIndex(F_var, getCurrentContext());
            F_fcn->setDataOnPatchHierarchy(F_current_idx, F_var, d_hierarchy, current_time+0.5*dt);
        }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Solve for Q(n+1).
    ////////////////////////////////////////////////////////////////////////////

    // Indicate that all solvers need to be reinitialized if the current
    // timestep size is different from the previous one.
    if (!MathUtilities<double>::equalEps(dt,d_old_dt))
    {
        std::fill(d_helmholtz_solvers_need_init.begin(),d_helmholtz_solvers_need_init.end(), true);
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
    }

    l = 0;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > Q_var   = *cit;
        Pointer<CellVariable<NDIM,double> > F_var   = d_Q_F_map[Q_var];
        Pointer<CellVariable<NDIM,double> > Psi_var = d_Q_Psi_map[Q_var];
        const double mu = d_Q_diffusion_coef[Q_var];
        const double lambda = d_Q_damping_coef[Q_var];
        const std::vector<RobinBcCoefStrategy<NDIM>*>& Q_bc_coef = d_Q_bc_coef[Q_var];

        Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();

        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_temp_idx = var_db->mapVariableAndContextToIndex(Q_var, d_temp_context);
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        const int F_current_idx = (F_var.isNull()
                                   ? -1
                                   : var_db->mapVariableAndContextToIndex(F_var, getCurrentContext()));
        const int Psi_temp_idx = var_db->mapVariableAndContextToIndex(Psi_var, d_temp_context);

        // Allocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_temp_idx, current_time);
            level->allocatePatchData(Psi_temp_idx, new_time);
        }

        if (!F_var.isNull())
        {
            d_hier_cc_data_ops->add(Q_new_idx, F_current_idx, Q_new_idx);
        }

        // Setup the problem coefficients and right hand side for the linear
        // solve for Q(n+1).
        PoissonSpecifications& helmholtz_spec = d_helmholtz_specs[l];
        Pointer<CCLaplaceOperator> helmholtz_op = d_helmholtz_ops[l];
        switch (d_viscous_timestepping_type)
        {
            case BACKWARD_EULER:
            {
                // The backward Euler discretization is:
                //
                //     (I-dt*mu*L(t_new)) Q(n+1) = Q(n) + F(t_avg) dt
                //
                // where
                //
                //    t_new = (n+1) dt
                //    t_avg = (t_new+t_old)/2
                //
                // Note that for simplicity of implementation, we always use a
                // timestep-centered forcing term.
                helmholtz_spec.setCConstant(1.0+dt*lambda);
                helmholtz_spec.setDConstant(   -dt*mu    );

                PoissonSpecifications rhs_spec("rhs_spec");
                rhs_spec.setCConstant(1.0);
                rhs_spec.setDConstant(0.0);

                d_hier_cc_data_ops->copyData(Q_temp_idx, Q_current_idx, false);
                d_hier_bdry_fill_ops[l]->setHomogeneousBc(false);
                d_hier_bdry_fill_ops[l]->fillData(current_time);

                for (int depth = 0; depth < Q_depth; ++depth)
                {
                    d_hier_math_ops->laplace(
                        Psi_temp_idx, Psi_var,  // Psi(n+1/2)
                        rhs_spec,               // Poisson spec
                        Q_temp_idx  , Q_var  ,  // Q(n)
                        d_no_fill_op,           // don't need to re-fill Q(n) data
                        current_time,           // Q(n) bdry fill time
                        dt,                     // gamma
                        Q_new_idx   , Q_var  ,  // N(n+1/2) = (u*grad Q)(n+1/2)
                        depth, depth, depth);   // dst_depth, src1_depth, src2_depth
                }
                break;
            }
            case CRANK_NICOLSON:
            {
                // The Crank-Nicolson discretization is:
                //
                //     (I-0.5*dt*mu*L(t_new)) Q(n+1) = (I+0.5*dt*mu*L(t_old)) Q(n) + F(t_avg) dt
                //
                // where
                //
                //    t_old = n dt
                //    t_new = (n+1) dt
                //    t_avg = (t_new+t_old)/2
                helmholtz_spec.setCConstant(1.0+0.5*dt*lambda);
                helmholtz_spec.setDConstant(   -0.5*dt*mu    );

                PoissonSpecifications rhs_spec("rhs_spec");
                rhs_spec.setCConstant(1.0-0.5*dt*lambda);
                rhs_spec.setDConstant(   +0.5*dt*mu    );

                d_hier_cc_data_ops->copyData(Q_temp_idx, Q_current_idx, false);
                d_hier_bdry_fill_ops[l]->setHomogeneousBc(false);
                d_hier_bdry_fill_ops[l]->fillData(current_time);

                for (int depth = 0; depth < Q_depth; ++depth)
                {
                    d_hier_math_ops->laplace(
                        Psi_temp_idx, Psi_var,  // Psi(n+1/2)
                        rhs_spec,               // Poisson spec
                        Q_temp_idx  , Q_var  ,  // Q(n)
                        d_no_fill_op,           // don't need to re-fill Q(n) data
                        current_time,           // Q(n) bdry fill time
                        dt,                     // gamma
                        Q_new_idx   , Q_var  ,  // N(n+1/2) = (u*grad Q)(n+1/2)
                        depth, depth, depth);   // dst_depth, src1_depth, src2_depth
                }
                break;
            }
            default:
                TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                           << "  unrecognized viscous timestepping type: " << enum_to_string<ViscousTimesteppingType>(d_viscous_timestepping_type) << "." << std::endl);
        }

        // Initialize the linear solver.
        helmholtz_op->setPoissonSpecifications(helmholtz_spec);
        helmholtz_op->setPhysicalBcCoefs(Q_bc_coef);
        helmholtz_op->setHomogeneousBc(false);
        helmholtz_op->setTime(new_time);
        helmholtz_op->setHierarchyMathOps(d_hier_math_ops);

        Pointer<CCPoissonFACOperator> helmholtz_fac_op = d_helmholtz_fac_ops[l];
        Pointer<FACPreconditioner>    helmholtz_fac_pc = d_helmholtz_fac_pcs[l];
        Pointer<KrylovLinearSolver>   helmholtz_solver = d_helmholtz_solvers[l];

        if (d_helmholtz_solvers_need_init[l])
        {
            if (d_do_log) plog << d_object_name << ": "
                               << "Initializing Helmholtz solvers for variable number " << l
                               << ", dt = " << dt << "\n";
            if (d_using_FAC)
            {
                helmholtz_fac_op->setPoissonSpecifications(helmholtz_spec);
                helmholtz_fac_op->setTime(new_time);
                helmholtz_fac_op->setResetLevels(d_coarsest_reset_ln, d_finest_reset_ln);
            }
            helmholtz_solver->initializeSolverState(*d_sol_vecs[l],*d_rhs_vecs[l]);

            // Indicate that the solvers do not presently require
            // re-initialization.
            d_helmholtz_solvers_need_init[l] = false;
        }

        // Solve for Q(n+1).
        switch (d_viscous_timestepping_type)
        {
            case BACKWARD_EULER:
            case CRANK_NICOLSON:
                helmholtz_op->setTime(new_time);
                if (d_using_FAC) helmholtz_fac_op->setTime(new_time);
                helmholtz_solver->solveSystem(*d_sol_vecs[l],*d_rhs_vecs[l]);
                d_hier_cc_data_ops->copyData(Q_new_idx, Q_temp_idx);

                if (d_do_log) plog << "AdvDiffHierarchyIntegrator::integrateHierarchy(): linear solve number of iterations = " << helmholtz_solver->getNumIterations() << "\n";
                if (d_do_log) plog << "AdvDiffHierarchyIntegrator::integrateHierarchy(): linear solve residual norm        = " << helmholtz_solver->getResidualNorm()  << "\n";

                if (helmholtz_solver->getNumIterations() == helmholtz_solver->getMaxIterations())
                {
                    pout << d_object_name << "::integrateHierarchy():"
                         <<"  WARNING: linear solver iterations == max iterations\n";
                }
                break;
            default:
                TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                           << "  unrecognized viscous timestepping type: " << enum_to_string<ViscousTimesteppingType>(d_viscous_timestepping_type) << "." << std::endl);
        }

        // Deallocate temporary data.
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
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
    const Pointer<PatchHierarchy<NDIM> > hierarchy,
    const int coarsest_level,
    const int finest_level,
    const double sync_time,
    const bool initial_time)
{
    t_synchronize_new_levels->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level < finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = coarsest_level; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    // We use the HyperbolicLevelIntegrator to handle as much data management as
    // possible.
    d_hyp_level_integrator->synchronizeNewLevels(hierarchy, coarsest_level, finest_level, sync_time, initial_time);

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
        d_hyp_level_integrator->resetTimeDependentData(d_hierarchy->getPatchLevel(ln), d_integrator_time, d_gridding_alg->levelCanBeRefined(ln));
    }

    t_reset_time_dependent_hier_data->stop();
    return;
}// resetTimeDependentHierData

void
AdvDiffHierarchyIntegrator::resetHierDataToPreadvanceState()
{
    t_reset_hier_data_to_preadvance_state->start();

    // We use the HyperbolicLevelIntegrator to handle as much data management as
    // possible.
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
///  StandardTagAndInitStrategy abstract base class.
///

void
AdvDiffHierarchyIntegrator::initializeLevelData(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > base_old_level,
    const bool allocate_data)
{
    t_initialize_level_data->start();

    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const Pointer<PatchLevel<NDIM> > old_level = base_old_level;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    if (!old_level.isNull())
    {
        TBOX_ASSERT(level_number == old_level->getLevelNumber());
    }
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    // We use the HyperbolicLevelIntegrator to handle as much data management as
    // possible.
    d_hyp_level_integrator->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level, allocate_data);

    // Set the initial values of any forcing terms.  All other variables are
    // initialized by the hyperbolic level integrator.
    if (initial_time)
    {
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);
        for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_F_var.begin();
             cit != d_F_var.end(); ++cit)
        {
            Pointer<CellVariable<NDIM,double> > F_var = *cit;
            if (d_manage_F_data[F_var])
            {
                const int F_idx = var_db->mapVariableAndContextToIndex(F_var, getCurrentContext());
                Pointer<CartGridFunction> F_fcn = d_F_fcn[F_var];
                if (!F_fcn.isNull())
                {
                    F_fcn->setDataOnPatchLevel(F_idx, F_var, level, init_data_time, initial_time);
                }
                else
                {
                    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                    {
                        Pointer<Patch<NDIM> > patch = level->getPatch(p());
                        Pointer<CellData<NDIM,double> > F_data = patch->getPatchData(F_idx);
#ifdef DEBUG_CHECK_ASSERTIONS
                        TBOX_ASSERT(!F_data.isNull());
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
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    t_reset_hierarchy_configuration->start();

    const Pointer<BasePatchHierarchy<NDIM> > hierarchy = base_hierarchy;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((coarsest_level >= 0)
                && (coarsest_level <= finest_level)
                && (finest_level <= hierarchy->getFinestLevelNumber()));
    for (int ln = 0; ln <= finest_level; ++ln)
    {
        TBOX_ASSERT(!(hierarchy->getPatchLevel(ln)).isNull());
    }
#endif
    const int finest_hier_level = hierarchy->getFinestLevelNumber();

    // We use the HyperbolicLevelIntegrator to handle as much data management as
    // possible.
    d_hyp_level_integrator->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);

    // Reset the Hierarchy data operations for the new hierarchy configuration.
    d_hier_cc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_cc_data_ops->resetLevels(0, finest_hier_level);

    // Reset the interpolation operators.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_hier_bdry_fill_ops.resize(d_Q_var.size());
    unsigned l = 0;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        const int Q_temp_idx = var_db->mapVariableAndContextToIndex(Q_var, d_temp_context);

        // Setup the interpolation transaction information.
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        InterpolationTransactionComponent transaction_comp(Q_temp_idx, DATA_COARSEN_TYPE, BDRY_EXTRAP_TYPE, CONSISTENT_TYPE_2_BDRY, d_Q_bc_coef[Q_var]);

        // Initialize the interpolation operators.
        d_hier_bdry_fill_ops[l] = new HierarchyGhostCellInterpolation();
        d_hier_bdry_fill_ops[l]->initializeOperatorState(transaction_comp, d_hierarchy);
    }

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
    d_sol_vecs.resize(d_Q_var.size());
    d_rhs_vecs.resize(d_Q_var.size());
    l = 0;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        const std::string& name = Q_var->getName();

        const int Q_temp_idx = var_db->mapVariableAndContextToIndex(Q_var, d_temp_context);
        d_sol_vecs[l] = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::sol_vec::"+name, d_hierarchy, 0, finest_hier_level);
        d_sol_vecs[l]->addComponent(Q_var,Q_temp_idx,d_wgt_idx,d_hier_cc_data_ops);

        Pointer<CellVariable<NDIM,double> > Psi_var = d_Q_Psi_map[Q_var];
        const int Psi_temp_idx = var_db->mapVariableAndContextToIndex(Psi_var, d_temp_context);
        d_rhs_vecs[l] = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::rhs_vec::"+name, d_hierarchy, 0, finest_hier_level);
        d_rhs_vecs[l]->addComponent(Psi_var,Psi_temp_idx,d_wgt_idx,d_hier_cc_data_ops);
    }

    // Indicate that all linear solvers must be re-initialized.
    std::fill(d_helmholtz_solvers_need_init.begin(), d_helmholtz_solvers_need_init.end(), true);
    d_coarsest_reset_ln = coarsest_level;
    d_finest_reset_ln = finest_level;

    // If we have added or removed a level, resize the schedule vectors.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin(); it != d_calgs.end(); ++it)
    {
        d_cscheds[(*it).first].resize(finest_hier_level+1);
    }

    // (Re)build coarsen communication schedules.  These are set only for levels
    // >= 1.
    for (CoarsenAlgMap::const_iterator it = d_calgs.begin(); it != d_calgs.end(); ++it)
    {
        for (int ln = std::max(coarsest_level,1); ln <= finest_hier_level; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            Pointer<PatchLevel<NDIM> > coarser_level = hierarchy->getPatchLevel(ln-1);
            d_cscheds[(*it).first][ln] = (*it).second->createSchedule(coarser_level, level, d_cstrategies[(*it).first]);
        }
    }

    t_reset_hierarchy_configuration->stop();
    return;
}// resetHierarchyConfiguration

void
AdvDiffHierarchyIntegrator::applyGradientDetector(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    t_apply_gradient_detector->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0)
                && (level_number <= hierarchy->getFinestLevelNumber()));
    TBOX_ASSERT(!(hierarchy->getPatchLevel(level_number)).isNull());
#endif
    Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(level_number);

    // Due to tag buffers, it is necessary to untag all cells prior to tagging.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        Pointer<CellData<NDIM,int> > tags_data = patch->getPatchData(tag_index);
        tags_data->fillAll(0);
    }

    // Tag cells for refinement according to the criteria specified by the
    // criteria specified by the HyperbolicLevelIntegrator<NDIM>.
    d_hyp_level_integrator->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);

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
///  allow access to the various variable contexts maintained by the integrator.
///

///
/// We simply reuse the VariableContext objects defined in the
/// HyperbolicLevelIntegrator<NDIM> object.
///

Pointer<VariableContext>
AdvDiffHierarchyIntegrator::getCurrentContext() const
{
    return d_hyp_level_integrator->getCurrentContext();
}// getCurrentContext

Pointer<VariableContext>
AdvDiffHierarchyIntegrator::getNewContext() const
{
    return d_hyp_level_integrator->getNewContext();
}// getNewContext

Pointer<VariableContext>
AdvDiffHierarchyIntegrator::getOldContext() const
{
    return d_hyp_level_integrator->getOldContext();
}// getOldContext

Pointer<VariableContext>
AdvDiffHierarchyIntegrator::getScratchContext() const
{
    return d_hyp_level_integrator->getScratchContext();
}// getScratchContext

Pointer<VariableContext>
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
///  Serializable abstract base class.
///

void
AdvDiffHierarchyIntegrator::putToDatabase(
    Pointer<Database> db)
{
    t_put_to_database->start();

#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif

    db->putInteger("ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION", ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION);
    db->putString("d_viscous_timestepping_type", enum_to_string<ViscousTimesteppingType>(d_viscous_timestepping_type));
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

/////////////////////////////// PRIVATE //////////////////////////////////////

void
AdvDiffHierarchyIntegrator::getFromInput(
    Pointer<Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    // Read in data members from input database.
    d_end_time = db->getDoubleWithDefault("end_time", d_end_time);
    d_grow_dt = db->getDoubleWithDefault("grow_dt", d_grow_dt);

    d_max_integrator_steps = db->getIntegerWithDefault("max_integrator_steps", d_max_integrator_steps);

    d_regrid_interval = db->getIntegerWithDefault("regrid_interval", d_regrid_interval);
    d_regrid_mode = string_to_enum<RegridMode>(db->getStringWithDefault("regrid_mode", enum_to_string<RegridMode>(d_regrid_mode)));

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
        d_viscous_timestepping_type = string_to_enum<ViscousTimesteppingType>(
            db->getStringWithDefault("viscous_timestepping_type", enum_to_string<ViscousTimesteppingType>(d_viscous_timestepping_type)));
        d_start_time = db->getDoubleWithDefault("start_time", d_start_time);
    }
    return;
}// getFromInput

void
AdvDiffHierarchyIntegrator::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();

    Pointer<Database> db;
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

    d_viscous_timestepping_type = string_to_enum<ViscousTimesteppingType>(db->getString("d_viscous_timestepping_type"));
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
template class Pointer<IBAMR::AdvDiffHierarchyIntegrator>;

//////////////////////////////////////////////////////////////////////////////
