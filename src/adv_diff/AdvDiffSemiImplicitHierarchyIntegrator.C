// Filename: AdvDiffSemiImplicitHierarchyIntegrator.C
// Created on 22 May 2012 by Boyce Griffith
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

#include "AdvDiffSemiImplicitHierarchyIntegrator.h"

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
#include <ibamr/AdvDiffConvectiveOperatorManager.h>
#include <ibamr/ibamr_utilities.h>
#include <ibamr/namespaces.h>

// SAMRAI INCLUDES
#include <HierarchyDataOpsManager.h>
#include <tbox/NullDatabase.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

// Number of ghosts cells used for each variable quantity.
static const int CELLG = (USING_LARGE_GHOST_CELL_WIDTH ? 2 : 1);

/////////////////////////////// PUBLIC ///////////////////////////////////////

AdvDiffSemiImplicitHierarchyIntegrator::AdvDiffSemiImplicitHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    bool register_for_restart)
    : AdvDiffHierarchyIntegrator(object_name, input_db, register_for_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!object_name.empty());
    TBOX_ASSERT(!input_db.isNull());
#endif
    // Default values.
    d_default_convective_time_stepping_type = MIDPOINT_RULE;
    d_default_init_convective_time_stepping_type = MIDPOINT_RULE;
    d_default_convective_op_type = AdvDiffConvectiveOperatorManager::DEFAULT;
    d_default_convective_bdry_extrap_type = "LINEAR";

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Check to make sure the time stepping types are supported.
    switch (d_default_diffusion_time_stepping_type)
    {
        case BACKWARD_EULER:
        case FORWARD_EULER:
        case TRAPEZOIDAL_RULE:
            break;
        default:
            TBOX_ERROR(d_object_name << "::AdvDiffSemiImplicitHierarchyIntegrator():\n"
                       << "  unsupported default diffusion time stepping type: " << enum_to_string<TimeSteppingType>(d_default_diffusion_time_stepping_type) << " \n"
                       << "  valid choices are: BACKWARD_EULER, FORWARD_EULER, TRAPEZOIDAL_RULE\n");
    }

    switch (d_default_convective_time_stepping_type)
    {
        case ADAMS_BASHFORTH:
        case FORWARD_EULER:
        case MIDPOINT_RULE:
        case TRAPEZOIDAL_RULE:
            break;
        default:
            TBOX_ERROR(d_object_name << "::AdvDiffSemiImplicitHierarchyIntegrator():\n"
                       << "  unsupported default convective time stepping type: " << enum_to_string<TimeSteppingType>(d_default_convective_time_stepping_type) << " \n"
                       << "  valid choices are: ADAMS_BASHFORTH, FORWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    if (is_multistep_time_stepping_type(d_default_convective_time_stepping_type))
    {
        switch (d_default_init_convective_time_stepping_type)
        {
            case FORWARD_EULER:
            case MIDPOINT_RULE:
            case TRAPEZOIDAL_RULE:
                break;
            default:
                TBOX_ERROR(d_object_name << "::AdvDiffSemiImplicitHierarchyIntegrator():\n"
                           << "  unsupported default initial convective time stepping type: " << enum_to_string<TimeSteppingType>(d_default_init_convective_time_stepping_type) << " \n"
                           << "  valid choices are: FORWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
        }
    }
    return;
}// AdvDiffSemiImplicitHierarchyIntegrator

AdvDiffSemiImplicitHierarchyIntegrator::~AdvDiffSemiImplicitHierarchyIntegrator()
{
    // intentionally blank
    return;
}// ~AdvDiffSemiImplicitHierarchyIntegrator

void
AdvDiffSemiImplicitHierarchyIntegrator::setDefaultConvectiveTimeSteppingType(
    TimeSteppingType default_convective_time_stepping_type)
{
    d_default_convective_time_stepping_type = default_convective_time_stepping_type;
    return;
}// setDefaultConvectiveTimeSteppingType

TimeSteppingType
AdvDiffSemiImplicitHierarchyIntegrator::getDefaultConvectiveTimeSteppingType() const
{
    return d_default_convective_time_stepping_type;
}// getDefaultConvectiveTimeSteppingType

void
AdvDiffSemiImplicitHierarchyIntegrator::setDefaultInitialConvectiveTimeSteppingType(
    TimeSteppingType default_init_convective_time_stepping_type)
{
    d_default_init_convective_time_stepping_type = default_init_convective_time_stepping_type;
    return;
}// setDefaultInitialConvectiveTimeSteppingType

TimeSteppingType
AdvDiffSemiImplicitHierarchyIntegrator::getDefaultInitialConvectiveTimeSteppingType() const
{
    return d_default_init_convective_time_stepping_type;
}// getDefaultInitialConvectiveTimeSteppingType

void
AdvDiffSemiImplicitHierarchyIntegrator::setDefaultConvectiveOperatorType(
    const std::string& op_type)
{
    d_default_convective_op_type = op_type;
    return;
}// setDefaultConvectiveOperatorType

const std::string&
AdvDiffSemiImplicitHierarchyIntegrator::getDefaultConvectiveOperatorType() const
{
    return d_default_convective_op_type;
}// getDefaultConvectiveOperatorType

void
AdvDiffSemiImplicitHierarchyIntegrator::setDefaultConvectiveOperatorBoundaryExtrapolation(
    const std::string& bdry_extrap_type)
{
    d_default_convective_bdry_extrap_type = bdry_extrap_type;
    return;
}// setDefaultConvectiveOperatorBoundaryExtrapolation

const std::string&
AdvDiffSemiImplicitHierarchyIntegrator::getDefaultConvectiveOperatorBoundaryExtrapolation() const
{
    return d_default_convective_bdry_extrap_type;
}// getDefaultConvectiveOperatorBoundaryExtrapolation

void
AdvDiffSemiImplicitHierarchyIntegrator::registerTransportedQuantity(
    Pointer<CellVariable<NDIM,double> > Q_var)
{
    AdvDiffHierarchyIntegrator::registerTransportedQuantity(Q_var);

    // Set default values.
    d_Q_convective_time_stepping_type[Q_var] = d_default_convective_time_stepping_type;
    d_Q_init_convective_time_stepping_type[Q_var] = d_default_init_convective_time_stepping_type;
    d_Q_convective_op_type[Q_var] = d_default_convective_op_type;
    d_Q_convective_bdry_extrap_type[Q_var] = d_default_convective_bdry_extrap_type;
    d_Q_convective_op[Q_var] = NULL;
    d_Q_convective_op_needs_init[Q_var] = false;
    return;
}// registerTransportedQuantity

void
AdvDiffSemiImplicitHierarchyIntegrator::setConvectiveTimeSteppingType(
    Pointer<CellVariable<NDIM,double> > Q_var,
    TimeSteppingType convective_time_stepping_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_convective_time_stepping_type[Q_var] = convective_time_stepping_type;
    return;
}// setConvectiveTimeSteppingType

TimeSteppingType
AdvDiffSemiImplicitHierarchyIntegrator::getConvectiveTimeSteppingType(
    Pointer<CellVariable<NDIM,double> > Q_var) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    return d_Q_convective_time_stepping_type.find(Q_var)->second;
}// getConvectiveTimeSteppingType

void
AdvDiffSemiImplicitHierarchyIntegrator::setInitialConvectiveTimeSteppingType(
    Pointer<CellVariable<NDIM,double> > Q_var,
    TimeSteppingType init_convective_time_stepping_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_init_convective_time_stepping_type[Q_var] = init_convective_time_stepping_type;
    return;
}// setInitialConvectiveTimeSteppingType

TimeSteppingType
AdvDiffSemiImplicitHierarchyIntegrator::getInitialConvectiveTimeSteppingType(
    Pointer<CellVariable<NDIM,double> > Q_var) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    return d_Q_init_convective_time_stepping_type.find(Q_var)->second;
}// getInitialConvectiveTimeSteppingType

void
AdvDiffSemiImplicitHierarchyIntegrator::setConvectiveOperatorType(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const std::string& op_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_convective_op_type[Q_var] = op_type;
    return;
}// setConvectiveOperatorType

const std::string&
AdvDiffSemiImplicitHierarchyIntegrator::getConvectiveOperatorType(
    Pointer<CellVariable<NDIM,double> > Q_var) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    return d_Q_convective_op_type.find(Q_var)->second;
}// getConvectiveOperatorType

void
AdvDiffSemiImplicitHierarchyIntegrator::setConvectiveBoundaryExtrapolation(
    Pointer<CellVariable<NDIM,double> > Q_var,
    const std::string& bdry_extrap_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    d_Q_convective_bdry_extrap_type[Q_var] = bdry_extrap_type;
    return;
}// setConvectiveBoundaryExtrapolation

const std::string&
AdvDiffSemiImplicitHierarchyIntegrator::getConvectiveBoundaryExtrapolation(
    Pointer<CellVariable<NDIM,double> > Q_var) const
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    return d_Q_convective_bdry_extrap_type.find(Q_var)->second;
}// getConvectiveBoundaryExtrapolation

void
AdvDiffSemiImplicitHierarchyIntegrator::setConvectiveOperator(
    Pointer<CellVariable<NDIM,double> > Q_var,
    Pointer<ConvectiveOperator> convective_op)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_Q_convective_op[Q_var] = convective_op;
    d_Q_convective_op_needs_init[Q_var] = true;
    return;
}// setConvectiveOperator

Pointer<ConvectiveOperator>
AdvDiffSemiImplicitHierarchyIntegrator::getConvectiveOperator(
    Pointer<CellVariable<NDIM,double> > Q_var)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(d_Q_var.find(Q_var) != d_Q_var.end());
#endif
    if (d_Q_convective_op[Q_var].isNull())
    {
        AdvDiffConvectiveOperatorManager* convective_op_manager = AdvDiffConvectiveOperatorManager::getManager();
        d_Q_convective_op[Q_var] = convective_op_manager->allocateOperator(
            d_Q_convective_op_type[Q_var], d_object_name+"::"+Q_var->getName()+"::ConvectiveOperator", Q_var, d_Q_difference_form[Q_var], d_Q_convective_bdry_extrap_type[Q_var]);
        d_Q_convective_op_needs_init[Q_var] = true;
    }
    return d_Q_convective_op[Q_var];
}// getConvectiveOperator

void
AdvDiffSemiImplicitHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    d_hierarchy = hierarchy;
    d_gridding_alg = gridding_alg;
    Pointer<CartesianGridGeometry<NDIM> > grid_geom = d_hierarchy->getGridGeometry();

    // Obtain the Hierarchy data operations objects.
    HierarchyDataOpsManager<NDIM>* hier_ops_manager = HierarchyDataOpsManager<NDIM>::getManager();
    d_hier_fc_data_ops = hier_ops_manager->getOperationsDouble(new FaceVariable<NDIM,double>("fc_var"), hierarchy, true);

    // Register variables using the default variable registration routine.
    AdvDiffHierarchyIntegrator::registerVariables();

    // Setup the convective operators.
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        getConvectiveOperator(Q_var);
    }

    // Register additional variables required for present time stepping algorithm.
    const IntVector<NDIM> cell_ghosts = CELLG;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();

        Pointer<CellVariable<NDIM,double> > N_var = new CellVariable<NDIM,double>(Q_var->getName()+"::N",Q_depth);
        d_N_var.insert(N_var);
        d_Q_N_map[Q_var] = N_var;
        int N_scratch_idx;
        registerVariable(N_scratch_idx, N_var, cell_ghosts, getScratchContext());

        Pointer<CellVariable<NDIM,double> > N_old_var = new CellVariable<NDIM,double>(Q_var->getName()+"::N_old",Q_depth);
        d_N_old_var.insert(N_old_var);
        d_Q_N_old_map[Q_var] = N_old_var;
        int N_old_current_idx, N_old_new_idx, N_old_scratch_idx;
        registerVariable(N_old_current_idx, N_old_new_idx, N_old_scratch_idx, N_old_var, cell_ghosts, "CONSERVATIVE_COARSEN", "CONSERVATIVE_LINEAR_REFINE");
    }

    // Perform hierarchy initialization operations common to all implementations
    // of AdvDiffHierarchyIntegrator.
    AdvDiffHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);

    // Indicate that the integrator has been initialized.
    d_integrator_is_initialized = true;
    return;
}// initializeHierarchyIntegrator

int
AdvDiffSemiImplicitHierarchyIntegrator::getNumberOfCycles() const
{
    int num_cycles = d_num_cycles;
    if (MathUtilities<double>::equalEps(d_integrator_time, d_start_time))
    {
        for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit)
        {
            Pointer<CellVariable<NDIM,double> > Q_var = *cit;
            if (d_Q_u_map.find(Q_var)->second.isNull()) continue;
            if (is_multistep_time_stepping_type(d_Q_convective_time_stepping_type.find(Q_var)->second) && d_Q_init_convective_time_stepping_type.find(Q_var)->second != FORWARD_EULER)
            {
                num_cycles = std::max(2,num_cycles);
            }
        }
    }
    return num_cycles;
}// getNumberOfCycles

void
AdvDiffSemiImplicitHierarchyIntegrator::preprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const int num_cycles)
{
    AdvDiffHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Indicate that all solvers need to be reinitialized if the current
    // timestep size is different from the previous one.
    const bool dt_change = initial_time || !MathUtilities<double>::equalEps(dt,d_dt_previous[0]);
    if (dt_change)
    {
        std::fill(d_helmholtz_solvers_need_init.begin(),d_helmholtz_solvers_need_init.end(), true);
        d_coarsest_reset_ln = 0;
        d_finest_reset_ln = finest_ln;
    }

    // Allocate the scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data    ,     new_time);
    }

    // Update the advection velocity.
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin(); cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        if (!d_u_fcn[u_var].isNull())
        {
            const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, getCurrentContext());
            d_u_fcn[u_var]->setDataOnPatchHierarchy(u_current_idx, u_var, d_hierarchy, current_time);
        }
    }

    // Setup the operators and solvers and compute the right-hand-side terms.
    unsigned int l = 0;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > Q_var     = *cit;
        Pointer<CellVariable<NDIM,double> > Q_rhs_var = d_Q_Q_rhs_map[Q_var];
        TimeSteppingType diffusion_time_stepping_type = d_Q_diffusion_time_stepping_type[Q_var];
        const double kappa  = d_Q_diffusion_coef[Q_var];
        const double lambda = d_Q_damping_coef  [Q_var];
        const std::vector<RobinBcCoefStrategy<NDIM>*>& Q_bc_coef = d_Q_bc_coef[Q_var];

        Pointer<CellDataFactory<NDIM,double> > Q_factory = Q_var->getPatchDataFactory();
        const int Q_depth = Q_factory->getDefaultDepth();

        const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
        const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
        const int Q_new_idx = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        const int Q_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(Q_rhs_var, getScratchContext());

        // Setup the problem coefficients and right-hand-side for the linear
        // solve for Q(n+1).
        double K = 0.0;
        switch (diffusion_time_stepping_type)
        {
            case BACKWARD_EULER:
                K = 1.0;
                break;
            case FORWARD_EULER:
                K = 0.0;
                break;
            case TRAPEZOIDAL_RULE:
                K = 0.5;
                break;
            default:
                TBOX_ERROR("this statment should not be reached");
        }
        PoissonSpecifications solver_spec(d_object_name+"::solver_spec::"+Q_var->getName());
        solver_spec.setCConstant(1.0+K*dt*lambda);
        solver_spec.setDConstant(   -K*dt*kappa );
        PoissonSpecifications rhs_spec(d_object_name+"::rhs_spec::"+Q_var->getName());
        rhs_spec.setCConstant(1.0-(1.0-K)*dt*lambda);
        rhs_spec.setDConstant(   +(1.0-K)*dt*kappa );
        d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_current_idx, false);
        d_hier_bdry_fill_ops[l]->setHomogeneousBc(false);
        d_hier_bdry_fill_ops[l]->fillData(current_time);
        for (int depth = 0; depth < Q_depth; ++depth)
        {
            d_hier_math_ops->laplace(Q_rhs_scratch_idx, Q_rhs_var, rhs_spec, Q_scratch_idx, Q_var, d_no_fill_op, current_time, 0.0, -1, NULL, depth, depth, depth);
        }

        // Initialize the linear solver.
        Pointer<PoissonSolver> helmholtz_solver = d_helmholtz_solvers[l];
        helmholtz_solver->setPoissonSpecifications(solver_spec);
        helmholtz_solver->setPhysicalBcCoefs(Q_bc_coef);
        helmholtz_solver->setHomogeneousBc(false);
        helmholtz_solver->setSolutionTime(new_time);
        helmholtz_solver->setTimeInterval(current_time, new_time);
        if (d_helmholtz_solvers_need_init[l])
        {
            if (d_enable_logging)
            {
                plog << d_object_name << ": "
                     << "Initializing Helmholtz solvers for variable number " << l
                     << ", dt = " << dt << "\n";
            }
            helmholtz_solver->initializeSolverState(*d_sol_vecs[l],*d_rhs_vecs[l]);
            d_helmholtz_solvers_need_init[l] = false;
        }

        // Account for the convective difference term.
        TimeSteppingType convective_time_stepping_type = d_Q_convective_time_stepping_type[Q_var];
        if (getIntegratorStep() == 0 && is_multistep_time_stepping_type(convective_time_stepping_type))
        {
            convective_time_stepping_type = d_Q_init_convective_time_stepping_type[Q_var];
        }
        Pointer<FaceVariable<NDIM,double> > u_var = d_Q_u_map[Q_var];
        Pointer<CellVariable<NDIM,double> > N_var = d_Q_N_map[Q_var];
        Pointer<CellVariable<NDIM,double> > N_old_var = d_Q_N_old_map[Q_var];
        if (!u_var.isNull())
        {
            if ((num_cycles == 1) && (convective_time_stepping_type == MIDPOINT_RULE || convective_time_stepping_type == TRAPEZOIDAL_RULE))
            {
                TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                           << "  time stepping type: " << enum_to_string<TimeSteppingType>(convective_time_stepping_type) << " requires num_cycles > 1.\n"
                           << "  at current time step, num_cycles = " << num_cycles << "\n");
            }
            if (d_Q_convective_op_needs_init[Q_var])
            {
                d_Q_convective_op[Q_var]->initializeOperatorState(*d_sol_vecs[l],*d_rhs_vecs[l]);
                d_Q_convective_op_needs_init[Q_var] = false;
            }
            const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, getCurrentContext());
            d_Q_convective_op[Q_var]->setAdvectionVelocity(u_current_idx);
            const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
            const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
            const int N_scratch_idx = var_db->mapVariableAndContextToIndex(N_var, getScratchContext());
            d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_current_idx);
            d_Q_convective_op[Q_var]->applyConvectiveOperator(Q_scratch_idx, N_scratch_idx);
            const int N_old_new_idx = var_db->mapVariableAndContextToIndex(N_old_var, getNewContext());
            d_hier_cc_data_ops->copyData(N_old_new_idx, N_scratch_idx);
            if (convective_time_stepping_type == FORWARD_EULER)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -1.0*dt, N_scratch_idx, Q_rhs_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -0.5*dt, N_scratch_idx, Q_rhs_scratch_idx);
            }
        }

        // Set the initial guess.
        d_hier_cc_data_ops->copyData(Q_new_idx, Q_current_idx);
    }
    return;
}// preprocessIntegrateHierarchy

void
AdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time,
    const int cycle_num)
{
    AdvDiffHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    const double dt  = new_time-current_time;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

    // Check to make sure that the number of cycles is what we expect it to be.
    const int expected_num_cycles = getNumberOfCycles();
    if (d_current_num_cycles != expected_num_cycles)
    {
        IBAMR_DO_ONCE(
            {
                pout << "AdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchy():\n"
                     << "  WARNING: num_cycles = " << d_current_num_cycles << " but expected num_cycles = " << expected_num_cycles << ".\n";
            }
                      );
    }

    // Perform a single step of fixed point iteration.
    unsigned int l = 0;
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit, ++l)
    {
        Pointer<CellVariable<NDIM,double> > Q_var     = *cit;
        Pointer<CellVariable<NDIM,double> > F_var     = d_Q_F_map    [Q_var];
        Pointer<CellVariable<NDIM,double> > Q_rhs_var = d_Q_Q_rhs_map[Q_var];

        const int Q_scratch_idx     = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
        const int Q_new_idx         = var_db->mapVariableAndContextToIndex(Q_var, getNewContext());
        const int F_scratch_idx     = d_F_fcn[F_var].isNull() ? -1 : var_db->mapVariableAndContextToIndex(F_var, getScratchContext());
        const int F_new_idx         = d_F_fcn[F_var].isNull() ? -1 : var_db->mapVariableAndContextToIndex(F_var, getNewContext());
        const int Q_rhs_scratch_idx = var_db->mapVariableAndContextToIndex(Q_rhs_var, getScratchContext());

        // Update the advection velocity.
        if (cycle_num > 0)
        {
            for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin(); cit != d_u_var.end(); ++cit)
            {
                Pointer<FaceVariable<NDIM,double> > u_var = *cit;
                const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, getCurrentContext());
                const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, getScratchContext());
                const int u_new_idx     = var_db->mapVariableAndContextToIndex(u_var, getNewContext()    );
                if (!d_u_fcn[u_var].isNull())
                {
                    d_u_fcn[u_var]->setDataOnPatchHierarchy(u_new_idx, u_var, d_hierarchy, new_time);
                }
                d_hier_fc_data_ops->linearSum(u_scratch_idx, 0.5, u_current_idx, 0.5, u_new_idx);
            }
        }

        // Account for the convective difference term.
        TimeSteppingType convective_time_stepping_type = d_Q_convective_time_stepping_type[Q_var];
        if (is_multistep_time_stepping_type(convective_time_stepping_type))
        {
#ifdef DEBUG_CHECK_ASSERTIONS
            TBOX_ASSERT(convective_time_stepping_type == ADAMS_BASHFORTH);
#endif
            if (getIntegratorStep() == 0)
            {
                convective_time_stepping_type = d_Q_init_convective_time_stepping_type[Q_var];
            }
            else if (cycle_num > 0)
            {
                convective_time_stepping_type = MIDPOINT_RULE;
                IBAMR_DO_ONCE(
                    {
                        pout << "AdvDiffSemiImplicitHierarchyIntegrator::integrateHierarchy():\n"
                             << "  WARNING: convective_time_stepping_type = " << enum_to_string<TimeSteppingType>(d_Q_convective_time_stepping_type[Q_var]) << " but num_cycles = " << d_current_num_cycles << " > 1.\n"
                             << "           using " << enum_to_string<TimeSteppingType>(d_Q_convective_time_stepping_type[Q_var]) << " only for the first cycle in each time step;\n"
                             << "           using " << enum_to_string<TimeSteppingType>(    convective_time_stepping_type       ) << " for subsequent cycles.\n";
                    }
                              );
            }
        }
        Pointer<FaceVariable<NDIM,double> > u_var = d_Q_u_map[Q_var];
        Pointer<CellVariable<NDIM,double> > N_var = d_Q_N_map[Q_var];
        Pointer<CellVariable<NDIM,double> > N_old_var = d_Q_N_old_map[Q_var];
        if (!u_var.isNull() && convective_time_stepping_type != FORWARD_EULER)
        {
            const int N_scratch_idx = var_db->mapVariableAndContextToIndex(N_var, getScratchContext());
            if (cycle_num > 0)
            {
                if (convective_time_stepping_type == MIDPOINT_RULE)
                {
                    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, getScratchContext());
                    d_Q_convective_op[Q_var]->setAdvectionVelocity(u_scratch_idx);
                    const int Q_current_idx = var_db->mapVariableAndContextToIndex(Q_var, getCurrentContext());
                    const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
                    const int Q_new_idx     = var_db->mapVariableAndContextToIndex(Q_var, getNewContext()    );
                    d_hier_cc_data_ops->linearSum(Q_scratch_idx, 0.5, Q_current_idx, 0.5, Q_new_idx);
                    d_Q_convective_op[Q_var]->applyConvectiveOperator(Q_scratch_idx, N_scratch_idx);
                }
                else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
                {
                    const int u_new_idx = var_db->mapVariableAndContextToIndex(u_var, getNewContext());
                    d_Q_convective_op[Q_var]->setAdvectionVelocity(u_new_idx);
                    const int Q_scratch_idx = var_db->mapVariableAndContextToIndex(Q_var, getScratchContext());
                    const int Q_new_idx     = var_db->mapVariableAndContextToIndex(Q_var, getNewContext()    );
                    d_hier_cc_data_ops->copyData(Q_scratch_idx, Q_new_idx);
                    d_Q_convective_op[Q_var]->applyConvectiveOperator(Q_scratch_idx, N_scratch_idx);
                }
            }
            if (convective_time_stepping_type == ADAMS_BASHFORTH)
            {
#ifdef DEBUG_CHECK_ASSERTIONS
                TBOX_ASSERT(cycle_num == 0);
#endif
                const int N_old_current_idx = var_db->mapVariableAndContextToIndex(N_old_var, getCurrentContext());
                const double omega = dt / d_dt_previous[0];
                d_hier_cc_data_ops->linearSum(N_scratch_idx, 1.0 + 0.5*omega, N_scratch_idx, -0.5*omega, N_old_current_idx);
            }
            if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -1.0*dt, N_scratch_idx, Q_rhs_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -0.5*dt, N_scratch_idx, Q_rhs_scratch_idx);
            }
        }

        // Account for forcing terms.
        if (!d_F_fcn[F_var].isNull())
        {
            d_F_fcn[F_var]->setDataOnPatchHierarchy(F_scratch_idx, F_var, d_hierarchy, current_time+0.5*dt);
            d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, dt, F_scratch_idx, Q_rhs_scratch_idx);
        }

        // Solve for Q(n+1).
        Pointer<PoissonSolver> helmholtz_solver = d_helmholtz_solvers[l];
        helmholtz_solver->solveSystem(*d_sol_vecs[l],*d_rhs_vecs[l]);
        d_hier_cc_data_ops->copyData(Q_new_idx, Q_scratch_idx);
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): diffusion solve number of iterations = " << helmholtz_solver->getNumIterations() << "\n";
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): diffusion solve residual norm        = " << helmholtz_solver->getResidualNorm()  << "\n";
        if (helmholtz_solver->getNumIterations() == helmholtz_solver->getMaxIterations())
        {
            pout << d_object_name << "::integrateHierarchy():"
                 <<"  WARNING: linear solver iterations == max iterations\n";
        }

        // Reset the right-hand side vector.
        if (!u_var.isNull())
        {
            const int N_scratch_idx = var_db->mapVariableAndContextToIndex(N_var, getScratchContext());
            if (convective_time_stepping_type == ADAMS_BASHFORTH || convective_time_stepping_type == MIDPOINT_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, +1.0*dt, N_scratch_idx, Q_rhs_scratch_idx);
            }
            else if (convective_time_stepping_type == TRAPEZOIDAL_RULE)
            {
                d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, +0.5*dt, N_scratch_idx, Q_rhs_scratch_idx);
            }
        }
        if (!d_F_fcn[F_var].isNull())
        {
            d_hier_cc_data_ops->axpy(Q_rhs_scratch_idx, -dt, F_scratch_idx, Q_rhs_scratch_idx);
            d_hier_cc_data_ops->copyData(F_new_idx, F_scratch_idx);
        }
    }
    return;
}// integrateHierarchy

void
AdvDiffSemiImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const bool skip_synchronize_new_state_data,
    const int num_cycles)
{
    AdvDiffHierarchyIntegrator::postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    // Update the advection velocity.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    for (std::set<Pointer<FaceVariable<NDIM,double> > >::const_iterator cit = d_u_var.begin(); cit != d_u_var.end(); ++cit)
    {
        Pointer<FaceVariable<NDIM,double> > u_var = *cit;
        if (!d_u_fcn[u_var].isNull() && d_u_fcn[u_var]->isTimeDependent())
        {
            const int u_idx = var_db->mapVariableAndContextToIndex(u_var, getNewContext());
            d_u_fcn[u_var]->setDataOnPatchHierarchy(u_idx, u_var, d_hierarchy, new_time);
        }
    }
    return;
}// postprocessIntegrateHierarchy

/////////////////////////////// PROTECTED ////////////////////////////////////

void
AdvDiffSemiImplicitHierarchyIntegrator::resetHierarchyConfigurationSpecialized(
    const Pointer<BasePatchHierarchy<NDIM> > base_hierarchy,
    const int coarsest_level,
    const int finest_level)
{
    const Pointer<BasePatchHierarchy<NDIM> > hierarchy = base_hierarchy;
    const int finest_hier_level = hierarchy->getFinestLevelNumber();
    d_hier_fc_data_ops->setPatchHierarchy(hierarchy);
    d_hier_fc_data_ops->resetLevels(0, finest_hier_level);
    for (std::set<Pointer<CellVariable<NDIM,double> > >::const_iterator cit = d_Q_var.begin(); cit != d_Q_var.end(); ++cit)
    {
        Pointer<CellVariable<NDIM,double> > Q_var = *cit;
        d_Q_convective_op_needs_init[Q_var] = true;
    }
    AdvDiffHierarchyIntegrator::resetHierarchyConfigurationSpecialized(base_hierarchy, coarsest_level, finest_level);
    return;
}// resetHierarchyConfigurationSpecialized

void
AdvDiffSemiImplicitHierarchyIntegrator::putToDatabaseSpecialized(
    Pointer<Database> db)
{
    db->putString("d_default_convective_time_stepping_type", enum_to_string<TimeSteppingType>(d_default_convective_time_stepping_type));
    db->putString("d_default_init_convective_time_stepping_type", enum_to_string<TimeSteppingType>(d_default_init_convective_time_stepping_type));
    db->putString("d_default_convective_op_type", d_default_convective_op_type);
    db->putString("d_default_convective_bdry_extrap_type", d_default_convective_bdry_extrap_type);
    AdvDiffHierarchyIntegrator::putToDatabaseSpecialized(db);
    return;
}// putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
AdvDiffSemiImplicitHierarchyIntegrator::getFromInput(
    Pointer<Database> db,
    bool is_from_restart)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!db.isNull());
#endif
    if (!is_from_restart)
    {
        if      (db->keyExists("convective_time_stepping_type")) d_default_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("convective_time_stepping_type"));
        else if (db->keyExists("convective_timestepping_type" )) d_default_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("convective_timestepping_type") );
        else if (db->keyExists("default_convective_time_stepping_type")) d_default_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("default_convective_time_stepping_type"));
        else if (db->keyExists("default_convective_timestepping_type" )) d_default_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("default_convective_timestepping_type") );

        if      (db->keyExists("init_convective_time_stepping_type")) d_default_init_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("init_convective_time_stepping_type"));
        else if (db->keyExists("init_convective_timestepping_type" )) d_default_init_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("init_convective_timestepping_type") );
        else if (db->keyExists("default_init_convective_time_stepping_type")) d_default_init_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("default_init_convective_time_stepping_type"));
        else if (db->keyExists("default_init_convective_timestepping_type" )) d_default_init_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("default_init_convective_timestepping_type") );

        if      (db->keyExists("convective_op_type"))               d_default_convective_op_type = db->getString("convective_op_type");
        else if (db->keyExists("convective_operator_type"))         d_default_convective_op_type = db->getString("convective_operator_type");
        else if (db->keyExists("default_convective_op_type"))       d_default_convective_op_type = db->getString("default_convective_op_type");
        else if (db->keyExists("default_convective_operator_type")) d_default_convective_op_type = db->getString("default_convective_operator_type");

        if      (db->keyExists("convective_bdry_extrap_type"))         d_default_convective_bdry_extrap_type = db->getString("convective_bdry_extrap_type");
        else if (db->keyExists("default_convective_bdry_extrap_type")) d_default_convective_bdry_extrap_type = db->getString("default_convective_bdry_extrap_type");
    }
    return;
}// getFromInput

void
AdvDiffSemiImplicitHierarchyIntegrator::getFromRestart()
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
    d_default_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("d_default_convective_time_stepping_type"));
    d_default_init_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("d_default_init_convective_time_stepping_type"));
    d_default_convective_op_type = db->getString("d_default_convective_op_type");
    d_default_convective_bdry_extrap_type = db->getString("d_default_convective_bdry_extrap_type");
    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
