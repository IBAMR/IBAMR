// Filename: INSHierarchyIntegrator.C
// Created on 10 Aug 2011 by Boyce Griffith
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

#include "INSHierarchyIntegrator.h"

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
#include <ibamr/INSIntermediateVelocityBcCoef.h>
#include <ibamr/INSProjectionBcCoef.h>
#include <ibamr/namespaces.h>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of INSHierarchyIntegrator restart file data.
static const int INS_HIERARCHY_INTEGRATOR_VERSION = 2;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

INSHierarchyIntegrator::~INSHierarchyIntegrator()
{
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        delete d_U_star_bc_coefs[d];
        d_U_star_bc_coefs[d] = NULL;
    }
    delete d_Phi_bc_coef;
    d_Phi_bc_coef = NULL;
    return;
}// ~INSHierarchyIntegrator

TimeSteppingType
INSHierarchyIntegrator::getViscousTimeSteppingType() const
{
    return d_viscous_time_stepping_type;
}// getViscousTimeSteppingType

TimeSteppingType
INSHierarchyIntegrator::getConvectiveTimeSteppingType() const
{
    return d_convective_time_stepping_type;
}// getConvectiveTimeSteppingType

TimeSteppingType
INSHierarchyIntegrator::getInitialConvectiveTimeSteppingType() const
{
    return d_init_convective_time_stepping_type;
}// getInitialConvectiveTimeSteppingType

void
INSHierarchyIntegrator::registerAdvDiffHierarchyIntegrator(
    Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!adv_diff_hier_integrator.isNull());
#endif
    d_adv_diff_hier_integrator = adv_diff_hier_integrator;
    registerChildHierarchyIntegrator(d_adv_diff_hier_integrator);
    d_adv_diff_hier_integrator->registerAdvectionVelocity(d_U_adv_diff_var);
    d_adv_diff_hier_integrator->setAdvectionVelocityIsDivergenceFree(d_U_adv_diff_var, d_Q_fcn.isNull());
    return;
}// registerAdvDiffHierarchyIntegrator

void
INSHierarchyIntegrator::setStokesSpecifications(
    StokesSpecifications problem_coefs)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_problem_coefs = problem_coefs;
    return;
}// setStokesSpecifications

const StokesSpecifications*
INSHierarchyIntegrator::getStokesSpecifications() const
{
    return &d_problem_coefs;
}// getStokesSpecifications

void
INSHierarchyIntegrator::registerPhysicalBoundaryConditions(
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    d_bc_coefs = bc_coefs;
    return;
}// registerPhysicalBoundaryConditions

void
INSHierarchyIntegrator::registerVelocityInitialConditions(
    Pointer<CartGridFunction> U_init)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_U_init = U_init;
    return;
}// registerVelocityInitialConditions

void
INSHierarchyIntegrator::registerPressureInitialConditions(
    Pointer<CartGridFunction> P_init)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_P_init = P_init;
    return;
}// registerPressureInitialConditions

void
INSHierarchyIntegrator::registerBodyForceFunction(
    Pointer<CartGridFunction> F_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_F_fcn = F_fcn;
    return;
}// registerBodyForceFunction

void
INSHierarchyIntegrator::registerFluidSourceFunction(
    Pointer<CartGridFunction> Q_fcn)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_Q_fcn = Q_fcn;
    if (!d_adv_diff_hier_integrator.isNull())
    {
        Pointer<FaceVariable<NDIM,double> > u_var = getAdvectionVelocityVariable();
        d_adv_diff_hier_integrator->setAdvectionVelocityIsDivergenceFree(u_var, d_Q_fcn.isNull());
    }
    return;
}// registerFluidSourceFunction

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getVelocityVariable() const
{
    return d_U_var;
}// getVelocityVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getPressureVariable() const
{
    return d_P_var;
}// getPressureVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getBodyForceVariable() const
{
    return d_F_var;
}// getBodyForceVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getFluidSourceVariable() const
{
    return d_Q_var;
}// getFluidSourceVariable

Pointer<FaceVariable<NDIM,double> >
INSHierarchyIntegrator::getAdvectionVelocityVariable() const
{
    return d_U_adv_diff_var;
}// getAdvectionVelocityVariable

std::vector<RobinBcCoefStrategy<NDIM>*>
INSHierarchyIntegrator::getIntermediateVelocityBoundaryConditions() const
{
    return d_U_star_bc_coefs;
}// getIntermediateVelocityBoundaryConditions

RobinBcCoefStrategy<NDIM>*
INSHierarchyIntegrator::getProjectionBoundaryConditions() const
{
    return d_Phi_bc_coef;
}// getProjectionBoundaryConditions

void
INSHierarchyIntegrator::setCreepingFlow(
    bool creeping_flow)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_creeping_flow = creeping_flow;
    d_convective_op.setNull();
    d_default_convective_difference_form = UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE;
    return;
}// setCreepingFlow

bool
INSHierarchyIntegrator::getCreepingFlow() const
{
    return d_creeping_flow;
}// getCreepingFlow

void
INSHierarchyIntegrator::setDefaultConvectiveOperatorType(
    ConvectiveOperatorType op_type)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(d_convective_op.isNull());
    TBOX_ASSERT(!d_creeping_flow);
#endif
    d_default_convective_op_type = op_type;
    return;
}// setDefaultConvectiveOperatorType

ConvectiveOperatorType
INSHierarchyIntegrator::getDefaultConvectiveOperatorType() const
{
    return d_default_convective_op_type;
}// getDefaultConvectiveOperatorType

void
INSHierarchyIntegrator::setDefaultConvectiveDifferencingType(
    ConvectiveDifferencingType difference_form)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(d_convective_op.isNull());
    TBOX_ASSERT(!d_creeping_flow);
#endif
    d_default_convective_difference_form = difference_form;
    return;
}// setDefaultConvectiveDifferencingType

ConvectiveDifferencingType
INSHierarchyIntegrator::getDefaultConvectiveDifferencingType() const
{
    return d_default_convective_difference_form;
}// getDefaultConvectiveDifferencingType

void
INSHierarchyIntegrator::setConvectiveOperator(
    Pointer<ConvectiveOperator> convective_op)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_convective_op = convective_op;
    d_creeping_flow = d_convective_op.isNull();
    return;
}// setConvectiveOperator

void
INSHierarchyIntegrator::setVelocitySubdomainSolver(
    Pointer<PoissonSolver> velocity_solver)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(d_velocity_solver.isNull());
#endif
    d_velocity_solver = velocity_solver;
    return;
}// setVelocitySubdomainSolver

void
INSHierarchyIntegrator::setPressureSubdomainSolver(
    Pointer<PoissonSolver> pressure_solver)
{
#ifdef DEBUG_CHECK_ASSERTIONS
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(d_pressure_solver.isNull());
#endif
    d_pressure_solver = pressure_solver;
    return;
}// setPressureSubdomainSolver

int
INSHierarchyIntegrator::getNumberOfCycles() const
{
    int num_cycles = d_num_cycles;
    if (!d_creeping_flow && MathUtilities<double>::equalEps(d_integrator_time, d_start_time) && is_multistep_time_stepping_type(d_convective_time_stepping_type) && d_init_convective_time_stepping_type != FORWARD_EULER)
    {
        num_cycles = std::max(2,num_cycles);
    }
    return num_cycles;
}// getNumberOfCycles

/////////////////////////////// PROTECTED ////////////////////////////////////

INSHierarchyIntegrator::INSHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<Variable<NDIM> > U_var,
    Pointer<Variable<NDIM> > P_var,
    Pointer<Variable<NDIM> > F_var,
    Pointer<Variable<NDIM> > Q_var,
    bool register_for_restart)
    : HierarchyIntegrator(object_name, input_db, register_for_restart),
      d_U_var(U_var),
      d_P_var(P_var),
      d_F_var(F_var),
      d_Q_var(Q_var),
      d_U_init(NULL),
      d_P_init(NULL),
      d_default_bc_coefs(d_object_name+"::default_bc_coefs", Pointer<Database>(NULL)),
      d_bc_coefs(NDIM,static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_F_fcn(NULL),
      d_Q_fcn(NULL)
{
    // Set some default values.
    d_integrator_is_initialized = false;
    d_viscous_time_stepping_type = TRAPEZOIDAL_RULE;
    d_convective_time_stepping_type = ADAMS_BASHFORTH;
    d_init_convective_time_stepping_type = MIDPOINT_RULE;
    d_num_cycles = 1;
    d_cfl_max = 1.0;
    d_using_vorticity_tagging = false;
    d_Omega_max = 0.0;
    d_normalize_pressure = false;
    d_default_convective_op_type = UNKNOWN_CONVECTIVE_OPERATOR_TYPE;
    d_default_convective_difference_form = ADVECTIVE;
    d_default_convective_bdry_extrap_type = "LINEAR";
    d_creeping_flow = false;
    d_regrid_max_div_growth_factor = 1.1;
    d_U_scale = 1.0;
    d_P_scale = 1.0;
    d_F_scale = 1.0;
    d_Q_scale = 1.0;
    d_Omega_scale = 1.0;
    d_Div_U_scale = 1.0;
    d_output_U = true;
    d_output_P = true;
    d_output_F = false;
    d_output_Q = false;
    d_output_Omega = true;
    d_output_Div_U = true;
    d_velocity_solver = NULL;
    d_pressure_solver = NULL;

    // Setup default boundary condition objects that specify homogeneous
    // Dirichlet (solid-wall) boundary conditions for the velocity.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_default_bc_coefs.setBoundaryValue(2*d  ,0.0);
        d_default_bc_coefs.setBoundaryValue(2*d+1,0.0);
    }
    registerPhysicalBoundaryConditions(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM,&d_default_bc_coefs));

    // Setup physical boundary conditions objects.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_U_star_bc_coefs[d] = new INSIntermediateVelocityBcCoef(d,&d_problem_coefs,d_bc_coefs);
    }
    d_Phi_bc_coef = new INSProjectionBcCoef(&d_problem_coefs,d_bc_coefs);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (!input_db.isNull()) getFromInput(input_db, from_restart);

    // Initialize an advection velocity variable.  NOTE: Patch data are
    // allocated for this variable only when an advection-diffusion solver is
    // registered with the INSHierarchyIntegrator.
    d_U_adv_diff_var = new FaceVariable<NDIM,double>(d_object_name+"::U_adv_diff");
    return;
}// INSHierarchyIntegrator

double
INSHierarchyIntegrator::getMaximumTimeStepSizeSpecialized()
{
    double dt = d_dt_max;
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt = std::min(dt, d_cfl_max*getStableTimestep(level));
    }
    const bool initial_time = MathUtilities<double>::equalEps(d_integrator_time, d_start_time);
    if (!initial_time && d_dt_growth_factor >= 1.0)
    {
        dt = std::min(dt,d_dt_growth_factor*d_dt_previous[0]);
    }
    return dt;
}// getMaximumTimeStepSizeSpecialized

double
INSHierarchyIntegrator::getStableTimestep(
    Pointer<PatchLevel<NDIM> > level) const
{
    double stable_dt = std::numeric_limits<double>::max();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        stable_dt = std::min(stable_dt,getStableTimestep(patch));
    }
    stable_dt = SAMRAI_MPI::minReduction(stable_dt);
    return stable_dt;
}// getStableTimestep

void
INSHierarchyIntegrator::putToDatabaseSpecialized(
    Pointer<Database> db)
{
    db->putInteger("INS_HIERARCHY_INTEGRATOR_VERSION",INS_HIERARCHY_INTEGRATOR_VERSION);
    db->putString("d_viscous_time_stepping_type", enum_to_string<TimeSteppingType>(d_viscous_time_stepping_type));
    db->putString("d_convective_time_stepping_type", enum_to_string<TimeSteppingType>(d_convective_time_stepping_type));
    db->putString("d_init_convective_time_stepping_type", enum_to_string<TimeSteppingType>(d_init_convective_time_stepping_type));
    db->putDouble("d_rho",d_problem_coefs.getRho());
    db->putDouble("d_mu",d_problem_coefs.getMu());
    db->putDouble("d_lambda",d_problem_coefs.getLambda());
    db->putDouble("d_cfl_max",d_cfl_max);
    db->putBool("d_using_vorticity_tagging",d_using_vorticity_tagging);
    if (d_Omega_rel_thresh.size() > 0) db->putDoubleArray("d_Omega_rel_thresh",d_Omega_rel_thresh);
    if (d_Omega_abs_thresh.size() > 0) db->putDoubleArray("d_Omega_abs_thresh",d_Omega_abs_thresh);
    db->putDouble("d_Omega_max",d_Omega_max);
    db->putBool("d_normalize_pressure",d_normalize_pressure);
    db->putString("d_default_convective_op_type",enum_to_string<ConvectiveOperatorType>(d_default_convective_op_type));
    db->putString("d_default_convective_difference_form",enum_to_string<ConvectiveDifferencingType>(d_default_convective_difference_form));
    db->putBool("d_creeping_flow",d_creeping_flow);
    db->putDouble("d_regrid_max_div_growth_factor",d_regrid_max_div_growth_factor);
    db->putDouble("d_U_scale",d_U_scale);
    db->putDouble("d_P_scale",d_P_scale);
    db->putDouble("d_F_scale",d_F_scale);
    db->putDouble("d_Q_scale",d_Q_scale);
    db->putDouble("d_Omega_scale",d_Omega_scale);
    db->putDouble("d_Div_U_scale",d_Div_U_scale);
    db->putBool("d_output_U",d_output_U);
    db->putBool("d_output_P",d_output_P);
    db->putBool("d_output_F",d_output_F);
    db->putBool("d_output_Q",d_output_Q);
    db->putBool("d_output_Omega",d_output_Omega);
    db->putBool("d_output_Div_U",d_output_Div_U);
    return;
}// putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSHierarchyIntegrator::getFromInput(
    Pointer<Database> db,
    const bool is_from_restart)
{
    if (!is_from_restart)
    {
        if      (db->keyExists("viscous_time_stepping_type")) d_viscous_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("viscous_time_stepping_type"));
        else if (db->keyExists("viscous_timestepping_type") ) d_viscous_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("viscous_timestepping_type") );
        if      (db->keyExists("convective_time_stepping_type")) d_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("convective_time_stepping_type"));
        else if (db->keyExists("convective_timestepping_type") ) d_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("convective_timestepping_type") );
        if      (db->keyExists("init_convective_time_stepping_type")) d_init_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("init_convective_time_stepping_type"));
        else if (db->keyExists("init_convective_timestepping_type") ) d_init_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("init_convective_timestepping_type") );
        if (db->keyExists("rho"))
        {
            d_problem_coefs.setRho(db->getDouble("rho"));
        }
        else
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `rho' not found in input.");
        }

        if (db->keyExists("mu"))
        {
            d_problem_coefs.setMu(db->getDouble("mu"));
        }
        else
        {
            TBOX_ERROR(d_object_name << ":  "
                       << "Key data `mu' not found in input.");
        }

        if (db->keyExists("lambda"))
        {
            d_problem_coefs.setLambda(db->getDouble("lambda"));
        }
        else
        {
            d_problem_coefs.setLambda(0.0);
        }
    }
    if (db->keyExists("num_cycles")) d_num_cycles = db->getInteger("num_cycles");
    if (db->keyExists("cfl")) d_cfl_max = db->getDouble("cfl");
    else if (db->keyExists("cfl_max")) d_cfl_max = db->getDouble("cfl_max");
    else if (db->keyExists("CFL")) d_cfl_max = db->getDouble("CFL");
    else if (db->keyExists("CFL_max")) d_cfl_max = db->getDouble("CFL_max");
    if (db->keyExists("using_vorticity_tagging")) d_using_vorticity_tagging = db->getBool("using_vorticity_tagging");
    if (db->keyExists("Omega_rel_thresh")) d_Omega_rel_thresh = db->getDoubleArray("Omega_rel_thresh");
    else if (db->keyExists("omega_rel_thresh")) d_Omega_rel_thresh = db->getDoubleArray("omega_rel_thresh");
    else if (db->keyExists("vorticity_rel_thresh")) d_Omega_rel_thresh = db->getDoubleArray("vorticity_rel_thresh");
    if (db->keyExists("Omega_abs_thresh")) d_Omega_abs_thresh = db->getDoubleArray("Omega_abs_thresh");
    else if (db->keyExists("omega_abs_thresh")) d_Omega_abs_thresh = db->getDoubleArray("omega_abs_thresh");
    else if (db->keyExists("vorticity_abs_thresh")) d_Omega_abs_thresh = db->getDoubleArray("vorticity_abs_thresh");
    if (db->keyExists("normalize_pressure")) d_normalize_pressure = db->getBool("normalize_pressure");
    if      (db->keyExists("convective_op_type"))               d_default_convective_op_type = string_to_enum<ConvectiveOperatorType>(db->getString("convective_op_type"));
    else if (db->keyExists("convective_operator_type"))         d_default_convective_op_type = string_to_enum<ConvectiveOperatorType>(db->getString("convective_operator_type"));
    else if (db->keyExists("default_convective_op_type"))       d_default_convective_op_type = string_to_enum<ConvectiveOperatorType>(db->getString("default_convective_op_type"));
    else if (db->keyExists("default_convective_operator_type")) d_default_convective_op_type = string_to_enum<ConvectiveOperatorType>(db->getString("default_convective_operator_type"));
    if      (db->keyExists("convective_difference_form"))         d_default_convective_difference_form = string_to_enum<ConvectiveDifferencingType>(db->getString("convective_difference_form"));
    else if (db->keyExists("convective_difference_type"))         d_default_convective_difference_form = string_to_enum<ConvectiveDifferencingType>(db->getString("convective_difference_type"));
    else if (db->keyExists("default_convective_difference_form")) d_default_convective_difference_form = string_to_enum<ConvectiveDifferencingType>(db->getString("default_convective_difference_form"));
    else if (db->keyExists("default_convective_difference_type")) d_default_convective_difference_form = string_to_enum<ConvectiveDifferencingType>(db->getString("default_convective_difference_type"));
    if      (db->keyExists("convective_bdry_extrap_type"))         d_default_convective_bdry_extrap_type = db->getString("convective_bdry_extrap_type");
    else if (db->keyExists("default_convective_bdry_extrap_type")) d_default_convective_bdry_extrap_type = db->getString("default_convective_bdry_extrap_type");
    if (db->keyExists("creeping_flow")) d_creeping_flow = db->getBool("creeping_flow");
    if (db->keyExists("regrid_max_div_growth_factor")) d_regrid_max_div_growth_factor = db->getDouble("regrid_max_div_growth_factor");
    if (db->keyExists("U_scale")) d_U_scale = db->getDouble("U_scale");
    if (db->keyExists("P_scale")) d_P_scale = db->getDouble("P_scale");
    if (db->keyExists("F_scale")) d_F_scale = db->getDouble("F_scale");
    if (db->keyExists("Q_scale")) d_Q_scale = db->getDouble("Q_scale");
    if (db->keyExists("Omega_scale")) d_Omega_scale = db->getDouble("Omega_scale");
    if (db->keyExists("Div_U_scale")) d_Div_U_scale = db->getDouble("Div_U_scale");
    if (db->keyExists("output_U")) d_output_U = db->getBool("output_U");
    if (db->keyExists("output_P")) d_output_P = db->getBool("output_P");
    if (db->keyExists("output_F")) d_output_F = db->getBool("output_F");
    if (db->keyExists("output_Q")) d_output_Q = db->getBool("output_Q");
    if (db->keyExists("output_Omega")) d_output_Omega = db->getBool("output_Omega");
    if (db->keyExists("output_Div_U")) d_output_Div_U = db->getBool("output_Div_U");

    if (db->keyExists("velocity_solver_type")) d_velocity_solver_type = db->getString("velocity_solver_type");
    if (db->keyExists("velocity_precond_type")) d_velocity_precond_type = db->getString("velocity_precond_type");
    if (db->keyExists("velocity_solver_db")) d_velocity_solver_db = db->getDatabase("velocity_solver_db");
    if (db->keyExists("velocity_precond_db")) d_velocity_precond_db = db->getDatabase("velocity_precond_db");
    if (!d_velocity_solver_db.isNull() && !d_velocity_solver_db->keyExists("options_prefix"))
    {
        d_velocity_solver_db->putString("options_prexix", "velocity_");
    }
    if (!d_velocity_precond_db.isNull() && !d_velocity_precond_db->keyExists("options_prefix"))
    {
        d_velocity_precond_db->putString("options_prexix", "velocity_pc_");
    }

    if (db->keyExists("pressure_solver_type")) d_pressure_solver_type = db->getString("pressure_solver_type");
    if (db->keyExists("pressure_precond_type")) d_pressure_precond_type = db->getString("pressure_precond_type");
    if (db->keyExists("pressure_solver_db")) d_pressure_solver_db = db->getDatabase("pressure_solver_db");
    if (db->keyExists("pressure_precond_db")) d_pressure_precond_db = db->getDatabase("pressure_precond_db");
    if (!d_pressure_solver_db.isNull() && !d_pressure_solver_db->keyExists("options_prefix"))
    {
        d_pressure_solver_db->putString("options_prexix", "pressure_");
    }
    if (!d_pressure_precond_db.isNull() && !d_pressure_precond_db->keyExists("options_prefix"))
    {
        d_pressure_precond_db->putString("options_prexix", "pressure_pc_");
    }

    if (db->keyExists("regrid_projection_solver_type")) d_regrid_projection_solver_type = db->getString("regrid_projection_solver_type");
    if (db->keyExists("regrid_projection_precond_type")) d_regrid_projection_precond_type = db->getString("regrid_projection_precond_type");
    if (db->keyExists("regrid_projection_solver_db")) d_regrid_projection_solver_db = db->getDatabase("regrid_projection_solver_db");
    if (db->keyExists("regrid_projection_precond_db")) d_regrid_projection_precond_db = db->getDatabase("regrid_projection_precond_db");
    if (!d_regrid_projection_solver_db.isNull() && !d_regrid_projection_solver_db->keyExists("options_prefix"))
    {
        d_regrid_projection_solver_db->putString("options_prexix", "regrid_projection_");
    }
    if (!d_regrid_projection_precond_db.isNull() && !d_regrid_projection_precond_db->keyExists("options_prefix"))
    {
        d_regrid_projection_precond_db->putString("options_prexix", "regrid_projection_pc_");
    }
    return;
}// getFromInput

void
INSHierarchyIntegrator::getFromRestart()
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
    int ver = db->getInteger("INS_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != INS_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    d_viscous_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("d_viscous_time_stepping_type"));
    d_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("d_convective_time_stepping_type"));
    d_init_convective_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("d_init_convective_time_stepping_type"));
    d_problem_coefs.setRho(db->getDouble("d_rho"));
    d_problem_coefs.setMu(db->getDouble("d_mu"));
    d_problem_coefs.setLambda(db->getDouble("d_lambda"));
    d_num_cycles = db->getInteger("d_num_cycles");
    d_cfl_max = db->getDouble("d_cfl_max");
    d_using_vorticity_tagging = db->getBool("d_using_vorticity_tagging");
    if (db->keyExists("d_Omega_rel_thresh")) d_Omega_rel_thresh = db->getDoubleArray("d_Omega_rel_thresh");
    else d_Omega_rel_thresh.resizeArray(0);
    if (db->keyExists("d_Omega_abs_thresh")) d_Omega_abs_thresh = db->getDoubleArray("d_Omega_abs_thresh");
    else d_Omega_abs_thresh.resizeArray(0);
    d_Omega_max = db->getDouble("d_Omega_max");
    d_normalize_pressure = db->getBool("d_normalize_pressure");
    d_default_convective_op_type = string_to_enum<ConvectiveOperatorType>(db->getString("d_default_convective_op_type"));
    d_default_convective_difference_form = string_to_enum<ConvectiveDifferencingType>(db->getString("d_default_convective_difference_form"));
    d_creeping_flow = db->getBool("d_creeping_flow");
    d_regrid_max_div_growth_factor = db->getDouble("d_regrid_max_div_growth_factor");
    d_U_scale = db->getDouble("d_U_scale");
    d_P_scale = db->getDouble("d_P_scale");
    d_F_scale = db->getDouble("d_F_scale");
    d_Q_scale = db->getDouble("d_Q_scale");
    d_Omega_scale = db->getDouble("d_Omega_scale");
    d_Div_U_scale = db->getDouble("d_Div_U_scale");
    d_output_U = db->getBool("d_output_U");
    d_output_P = db->getBool("d_output_P");
    d_output_F = db->getBool("d_output_F");
    d_output_Q = db->getBool("d_output_Q");
    d_output_Omega = db->getBool("d_output_Omega");
    d_output_Div_U = db->getBool("d_output_Div_U");
    return;
}// getFromRestart

//////////////////////////////////////////////////////////////////////////////

}// namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
