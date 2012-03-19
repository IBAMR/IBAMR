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
#include <HierarchyDataOpsManager.h>
#include <tbox/NullDatabase.h>

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

// Type of coarsening to perform prior to setting coarse-fine boundary and
// physical boundary ghost cell values.
static const std::string DATA_COARSEN_TYPE = "CUBIC_COARSEN";

// Type of extrapolation to use at physical boundaries.
static const std::string BDRY_EXTRAP_TYPE = "LINEAR";

// Whether to enforce consistent interpolated values at Type 2 coarse-fine
// interface ghost cells.
static const bool CONSISTENT_TYPE_2_BDRY = false;

// Version of AdvDiffHierarchyIntegrator restart file data.
static const int ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION = 2;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffHierarchyIntegrator::AdvDiffHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<GodunovAdvector> explicit_predictor,
    bool register_for_restart)
    : HierarchyIntegrator(object_name, input_db, register_for_restart),
      d_viscous_timestepping_type(CRANK_NICOLSON),
      d_u_var(),
      d_u_is_div_free(),
      d_u_fcn(),
      d_F_var(),
      d_F_fcn(),
      d_Q_var(),
      d_Psi_var(),
      d_Q_u_map(),
      d_Q_F_map(),
      d_Q_Psi_map(),
      d_Q_difference_form(),
      d_Q_diffusion_coef(),
      d_Q_damping_coef(),
      d_Q_init(),
      d_Q_bc_coef(),
      d_hyp_level_integrator(NULL),
      d_hyp_level_integrator_db(NULL),
      d_hyp_patch_ops(NULL),
      d_hyp_patch_ops_db(NULL),
      d_explicit_predictor(explicit_predictor),
      d_integrator_is_initialized(false),
      d_hier_cc_data_ops(NULL),
      d_temp_context(),
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
    TBOX_ASSERT(!explicit_predictor.isNull());
#endif
    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Get initialization data for the hyperbolic patch strategy objects.
    if (input_db->keyExists("HyperbolicLevelIntegrator"))
    {
        d_hyp_level_integrator_db = input_db->getDatabase("HyperbolicLevelIntegrator");
    }
    else
    {
        d_hyp_level_integrator_db = new NullDatabase();
    }
    if (input_db->keyExists("AdvDiffHypPatchOps"))
    {
        d_hyp_patch_ops_db = input_db->getDatabase("AdvDiffHypPatchOps");
    }
    else
    {
        d_hyp_patch_ops_db = new NullDatabase();
    }

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
        else
        {
            d_fac_op_db = new NullDatabase();
        }

        if (input_db->keyExists("FACPreconditioner"))
        {
            d_fac_pc_db = input_db->getDatabase("FACPreconditioner");
        }
        else if (input_db->keyExists("FACPreconditioners"))
        {
            d_fac_pc_db = input_db->getDatabase("FACPreconditioners");
        }
        else
        {
            d_fac_pc_db = new NullDatabase();
        }
    }
    return;
}// AdvDiffHierarchyIntegrator

AdvDiffHierarchyIntegrator::~AdvDiffHierarchyIntegrator()
{
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

ViscousTimesteppingType
AdvDiffHierarchyIntegrator::getViscousTimesteppingType() const
{
    return d_viscous_timestepping_type;
}// getViscousTimesteppingType

void
AdvDiffHierarchyIntegrator::registerAdvectionVelocity(
    Pointer<FaceVariable<NDIM,double> > u_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!u_var.isNull());
#endif
    d_u_var.insert(u_var);
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
    Pointer<CellVariable<NDIM,double> > F_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!F_var.isNull());
#endif
    d_F_var.insert(F_var);
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
AdvDiffHierarchyIntegrator::registerTransportedQuantity(
    Pointer<CellVariable<NDIM,double> > Q_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!Q_var.isNull());
#endif
    d_Q_var.insert(Q_var);
    d_Q_difference_form[Q_var] = CONSERVATIVE;
    d_Q_diffusion_coef[Q_var] = 0.0;
    d_Q_damping_coef[Q_var] = 0.0;

    Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
    const int Q_depth = Q_factory->getDefaultDepth();
    d_Q_bc_coef[Q_var] = std::vector<RobinBcCoefStrategy<NDIM>*>(Q_depth,static_cast<RobinBcCoefStrategy<NDIM>*>(NULL));

    Pointer<CellVariable<NDIM,double> > Psi_var = new CellVariable<NDIM,double>(Q_var->getName()+"::Psi",Q_depth);
    d_Psi_var.insert(Psi_var);
    d_Q_Psi_map[Q_var] = Psi_var;
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
AdvDiffHierarchyIntegrator::setConvectiveDifferencingType(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const ConvectiveDifferencingType difference_form)
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
    const double kappa)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_diffusion_coef[Q_var] = kappa;
    return;
}// setDiffusionCoefficient

void
AdvDiffHierarchyIntegrator::setDampingCoefficient(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const double lambda)
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
AdvDiffHierarchyIntegrator::setPhysicalBcCoefs(
    Pointer<CellVariable<NDIM,double> > Q_var,
    RobinBcCoefStrategy<NDIM>* Q_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
    const unsigned int Q_depth = Q_factory->getDefaultDepth();
    TBOX_ASSERT(Q_depth == 1);
#endif
    d_Q_bc_coef[Q_var] = std::vector<RobinBcCoefStrategy<NDIM>*>(1,Q_bc_coef);
    return;
}// setPhysicalBcCoefs

void
AdvDiffHierarchyIntegrator::setPhysicalBcCoefs(
    Pointer<CellVariable<NDIM,double> > Q_var,
    std::vector<RobinBcCoefStrategy<NDIM>*> Q_bc_coef)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
    const unsigned int Q_depth = Q_factory->getDefaultDepth();
    TBOX_ASSERT(Q_depth == Q_bc_coef.size());
#endif
    d_Q_bc_coef[Q_var] = Q_bc_coef;
    return;
}// setPhysicalBcCoefs

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

void
AdvDiffHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Initialize the HyperbolicPatchStrategy and HyperbolicLevelIntegrator
    // objects that provide numerical routines for explicitly integrating the
    // advective terms.
    d_hyp_patch_ops = new AdvDiffHypPatchOps(
        d_object_name+"::AdvDiffHypPatchOps",
        d_hyp_patch_ops_db,
        d_explicit_predictor,
        grid_geom,
        d_registered_for_restart);
    d_hyp_level_integrator = new HyperbolicLevelIntegrator<NDIM>(
        d_object_name+"::HyperbolicLevelIntegrator",
        d_hyp_level_integrator_db,
        d_hyp_patch_ops,
        d_registered_for_restart,
        /*using_time_refinement*/ false);

    // Setup variable contexts.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_temp_context = var_db->getContext(d_object_name+"::TEMP_CONTEXT");
    d_current_context = d_hyp_level_integrator->getCurrentContext();
    d_scratch_context = d_hyp_level_integrator->getScratchContext();
    d_new_context = d_hyp_level_integrator->getNewContext();

    // Register the VisIt data writer with the patch strategy object, which is
    // the object that actually registers variables for plotting.
    if (!d_visit_writer.isNull())
    {
        d_hyp_patch_ops->registerVisItDataWriter(d_visit_writer);
    }

    // Setup hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<CellVariable<NDIM,double> > cc_var = new CellVariable<NDIM,double>("cc_var");
    d_hier_cc_data_ops = hier_ops_manager->getOperationsDouble(cc_var, d_hierarchy, true);

    // Register variables with the hyperbolic level integrator.
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin();
         cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        d_hyp_patch_ops->registerAdvectionVelocity(u_var);
        d_hyp_patch_ops->setAdvectionVelocityIsDivergenceFree(u_var,d_u_is_div_free[u_var]);
        if (!d_u_fcn[u_var].isNull()) d_hyp_patch_ops->setAdvectionVelocityFunction(u_var,d_u_fcn[u_var]);
    }

    const IntVector<NDIM> cell_ghosts = CELLG;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_F_var.begin();
         cit != d_F_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > F_var = *cit;
        d_hyp_level_integrator->registerVariable(
            F_var, cell_ghosts,
            HyperbolicLevelIntegrator<NDIM>::TIME_DEP,
            d_hierarchy->getGridGeometry(),
            "CONSERVATIVE_COARSEN",
            "CONSERVATIVE_LINEAR_REFINE");
    }

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Psi_var.begin();
         cit != d_Psi_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Psi_var = *cit;
        d_hyp_patch_ops->registerSourceTerm(Psi_var);
        var_db->registerVariableAndContext(Psi_var, d_temp_context, CELLG);
    }

    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        d_hyp_patch_ops->registerTransportedQuantity(Q_var);
        if (!d_Q_u_map[Q_var].isNull()) d_hyp_patch_ops->setAdvectionVelocity(Q_var,d_Q_u_map[Q_var]);
        d_hyp_patch_ops->setSourceTerm(Q_var,d_Q_Psi_map[Q_var]);
        d_hyp_patch_ops->setConvectiveDifferencingType(Q_var,d_Q_difference_form[Q_var]);
        if (!d_Q_init[Q_var].isNull()) d_hyp_patch_ops->setInitialConditions(Q_var,d_Q_init[Q_var]);
        if (!d_Q_bc_coef[Q_var].empty()) d_hyp_patch_ops->setPhysicalBcCoefs(Q_var,d_Q_bc_coef[Q_var]);
        var_db->registerVariableAndContext(Q_var, d_temp_context, CELLG);
    }

    // Initialize the HyperbolicLevelIntegrator.
    //
    // NOTE: This must be done AFTER all variables have been registered with the
    // level integrator.
    d_hyp_level_integrator->initializeLevelIntegrator(d_gridding_alg);

    // Setup coarsening communications algorithms, used in synchronizing refined
    // regions of coarse data with the underlying fine data.
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        Pointer<CoarsenOperator<NDIM> > coarsen_operator = grid_geom->lookupCoarsenOperator(Q_var, "CONSERVATIVE_COARSEN");
        getCoarsenAlgorithm(SYNCH_CURRENT_DATA_ALG)->registerCoarsen(Q_current_idx, Q_current_idx, coarsen_operator);
        getCoarsenAlgorithm(SYNCH_NEW_DATA_ALG    )->registerCoarsen(Q_new_idx    , Q_new_idx    , coarsen_operator);
    }

    // Operators and solvers are maintained for each variable registered with the
    // integrator.
    d_helmholtz_specs.reserve(d_Q_var.size());
    d_helmholtz_ops.resize(d_Q_var.size());
    d_helmholtz_fac_ops.resize(d_Q_var.size());
    d_helmholtz_fac_pcs.resize(d_Q_var.size());
    d_helmholtz_solvers.resize(d_Q_var.size());
    d_helmholtz_solvers_need_init.resize(d_Q_var.size());
    unsigned int l = 0;
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
                                                           d_helmholtz_fac_ops[l], d_fac_pc_db);
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

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
}// initializeHierarchyIntegrator

void
AdvDiffHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time,
    const int cycle_num)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(current_time <= new_time);
    TBOX_ASSERT(d_end_time > d_integrator_time);
    TBOX_ASSERT(MathUtilities<double>::equalEps(d_integrator_time,current_time));
#endif
    const double dt = new_time - current_time;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    ////////////////////////////////////////////////////////////////////////////
    // Reset time-dependent data when necessary.
    ////////////////////////////////////////////////////////////////////////////

    if (cycle_num > 0)
    {
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            d_hyp_level_integrator->resetDataToPreadvanceState(d_hierarchy->getPatchLevel(ln));
        }
    }

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

    unsigned int l = 0;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > Q_var   = *cit;
        Pointer<CellVariable<NDIM,double> > F_var   = d_Q_F_map[Q_var];
        Pointer<CellVariable<NDIM,double> > Psi_var = d_Q_Psi_map[Q_var];
        const double kappa = d_Q_diffusion_coef[Q_var];
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

        PoissonSpecifications kappa_spec("kappa_spec");
        kappa_spec.setCConstant(-lambda);
        kappa_spec.setDConstant(+kappa );

        for (int depth = 0; depth < Q_depth; ++depth)
        {
            d_hier_math_ops->laplace(
                Psi_current_idx, Psi_var,  // Psi(n)
                kappa_spec,                // Poisson spec
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

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        static const bool first_step = true;
        static const bool last_step = false;
        d_hyp_level_integrator->advanceLevel(
            d_hierarchy->getPatchLevel(ln), d_hierarchy, current_time, new_time, first_step, last_step);
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
    if (initial_time || !MathUtilities<double>::equalEps(dt,d_dt_previous[0]))
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
        const double kappa = d_Q_diffusion_coef[Q_var];
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
                //     (I-dt*kappa*L(t_new)) Q(n+1) = Q(n) + F(t_avg) dt
                //
                // where
                //
                //    t_new = (n+1) dt
                //    t_avg = (t_new+t_old)/2
                //
                // Note that for simplicity of implementation, we always use a
                // timestep-centered forcing term.
                helmholtz_spec.setCConstant(1.0+dt*lambda);
                helmholtz_spec.setDConstant(   -dt*kappa );

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
                //     (I-0.5*dt*kappa*L(t_new)) Q(n+1) = (I+0.5*dt*kappa*L(t_old)) Q(n) + F(t_avg) dt
                //
                // where
                //
                //    t_old = n dt
                //    t_new = (n+1) dt
                //    t_avg = (t_new+t_old)/2
                helmholtz_spec.setCConstant(1.0+0.5*dt*lambda);
                helmholtz_spec.setDConstant(   -0.5*dt*kappa );

                PoissonSpecifications rhs_spec("rhs_spec");
                rhs_spec.setCConstant(1.0-0.5*dt*lambda);
                rhs_spec.setDConstant(   +0.5*dt*kappa );

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
    return;
}// integrateHierarchy

/////////////////////////////// PROTECTED ////////////////////////////////////

double
AdvDiffHierarchyIntegrator::getTimeStepSizeSpecialized()
{
    double dt = d_dt_max;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt = std::min(dt, d_hyp_level_integrator->getLevelDt(level, d_integrator_time, initial_time));
    }
    if (!initial_time && d_dt_growth_factor >= 1.0)
    {
        dt = std::min(dt,d_dt_growth_factor*d_dt_previous[0]);
    }
    return dt;
}// getTimeStepSizeSpecialized

void
AdvDiffHierarchyIntegrator::resetTimeDependentHierarchyDataSpecialized(
    const double new_time)
{
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Advance the simulation time.
    d_dt_previous.push_front(new_time - d_integrator_time);
    static const unsigned int MAX_DT_PREVIOUS_SIZE = 31;
    if (d_dt_previous.size() > MAX_DT_PREVIOUS_SIZE) d_dt_previous.pop_back();
    d_integrator_time = new_time;
    ++d_integrator_step;

    // Reset the time dependent data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hyp_level_integrator->resetTimeDependentData(
            d_hierarchy->getPatchLevel(ln), d_integrator_time, d_gridding_alg->levelCanBeRefined(ln));
    }
    return;
}// resetTimeDependentHierarchyDataSpecialized

void
AdvDiffHierarchyIntegrator::resetIntegratorToPreadvanceStateSpecialized()
{
    // We use the HyperbolicLevelIntegrator to handle most data management.
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        d_hyp_level_integrator->resetDataToPreadvanceState(d_hierarchy->getPatchLevel(ln));
    }
    return;
}// resetIntegratorToPreadvanceStateSpecialized

void
AdvDiffHierarchyIntegrator::initializeLevelDataSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int level_number,
    const double init_data_time,
    const bool can_be_refined,
    const bool initial_time,
    const Pointer<BasePatchLevel<NDIM> > base_old_level,
    const bool allocate_data)
{
    const Pointer<PatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const Pointer<PatchLevel<NDIM> > old_level = base_old_level;
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!hierarchy.isNull());
    TBOX_ASSERT((level_number >= 0) && (level_number <= hierarchy->getFinestLevelNumber()));
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
    return;
}// initializeLevelDataSpecialized

void
AdvDiffHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
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
    unsigned int l = 0;
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

    // Reset the solution and rhs vectors.
    d_sol_vecs.resize(d_Q_var.size());
    d_rhs_vecs.resize(d_Q_var.size());
    l = 0;
    const int wgt_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin();
         cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        const std::string& name = Q_var->getName();

        const int Q_temp_idx = var_db->mapVariableAndContextToIndex(Q_var, d_temp_context);
        d_sol_vecs[l] = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::sol_vec::"+name, d_hierarchy, 0, finest_hier_level);
        d_sol_vecs[l]->addComponent(Q_var,Q_temp_idx,wgt_idx,d_hier_cc_data_ops);

        Pointer<CellVariable<NDIM,double> > Psi_var = d_Q_Psi_map[Q_var];
        const int Psi_temp_idx = var_db->mapVariableAndContextToIndex(Psi_var, d_temp_context);
        d_rhs_vecs[l] = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::rhs_vec::"+name, d_hierarchy, 0, finest_hier_level);
        d_rhs_vecs[l]->addComponent(Psi_var,Psi_temp_idx,wgt_idx,d_hier_cc_data_ops);
    }

    // Indicate that all linear solvers must be re-initialized.
    std::fill(d_helmholtz_solvers_need_init.begin(), d_helmholtz_solvers_need_init.end(), true);
    d_coarsest_reset_ln = coarsest_level;
    d_finest_reset_ln = finest_level;
    return;
}// resetHierarchyConfigurationSpecialized

void
AdvDiffHierarchyIntegrator::applyGradientDetectorSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy,
    const int level_number,
    const double error_data_time,
    const int tag_index,
    const bool initial_time,
    const bool uses_richardson_extrapolation_too)
{
    // Tag cells for refinement according to the criteria specified by the
    // criteria specified by the level integrator.
    d_hyp_level_integrator->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time, uses_richardson_extrapolation_too);
    return;
}// applyGradientDetectorSpecialized

void
AdvDiffHierarchyIntegrator::putToDatabaseSpecialized(
    Pointer<Database> db)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    db->putInteger("ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION", ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION);
    db->putString("d_viscous_timestepping_type", enum_to_string<ViscousTimesteppingType>(d_viscous_timestepping_type));
    return;
}// putToDatabaseSpecialized

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
    if (!is_from_restart)
    {
        if (db->keyExists("viscous_timestepping_type")) d_viscous_timestepping_type = string_to_enum<ViscousTimesteppingType>(db->getString("viscous_timestepping_type"));
    }
    if (db->keyExists("max_iterations")) d_max_iterations = db->getInteger("max_iterations");
    if (db->keyExists("abs_residual_tol")) d_abs_residual_tol = db->getDouble("abs_residual_tol");
    if (db->keyExists("rel_residual_tol")) d_rel_residual_tol = db->getDouble("rel_residual_tol");
    if (db->keyExists("using_FAC")) d_using_FAC = db->getBool("using_FAC");
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
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to "
                   << d_object_name << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != ADV_DIFF_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    d_viscous_timestepping_type = string_to_enum<ViscousTimesteppingType>(db->getString("d_viscous_timestepping_type"));
    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
