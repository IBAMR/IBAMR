// Filename: INSHierarchyIntegrator.cpp
// Created on 10 Aug 2011 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
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
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <algorithm>
#include <deque>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include "FaceVariable.h"
#include "IntVector.h"
#include "LocationIndexRobinBcCoefs.h"
#include "MultiblockDataTranslator.h"
#include "Patch.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "RobinBcCoefStrategy.h"
#include "Variable.h"
#include "ibamr/AdvDiffHierarchyIntegrator.h"
#include "ibamr/ConvectiveOperator.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSIntermediateVelocityBcCoef.h"
#include "ibamr/INSProjectionBcCoef.h"
#include "ibamr/StokesSpecifications.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/CartGridFunction.h"
#include "ibtk/CartGridFunctionSet.h"
#include "ibtk/HierarchyGhostCellInterpolation.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/PoissonSolver.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/MathUtilities.h"
#include "tbox/MemoryDatabase.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

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
} // ~INSHierarchyIntegrator

TimeSteppingType
INSHierarchyIntegrator::getViscousTimeSteppingType() const
{
    return d_viscous_time_stepping_type;
} // getViscousTimeSteppingType

TimeSteppingType
INSHierarchyIntegrator::getConvectiveTimeSteppingType() const
{
    return d_convective_time_stepping_type;
} // getConvectiveTimeSteppingType

TimeSteppingType
INSHierarchyIntegrator::getInitialConvectiveTimeSteppingType() const
{
    return d_init_convective_time_stepping_type;
} // getInitialConvectiveTimeSteppingType

void
INSHierarchyIntegrator::registerAdvDiffHierarchyIntegrator(Pointer<AdvDiffHierarchyIntegrator> adv_diff_hier_integrator)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(adv_diff_hier_integrator);
#endif
    d_adv_diff_hier_integrator = adv_diff_hier_integrator;
    registerChildHierarchyIntegrator(d_adv_diff_hier_integrator);
    d_adv_diff_hier_integrator->registerAdvectionVelocity(d_U_adv_diff_var);
    d_adv_diff_hier_integrator->setAdvectionVelocityIsDivergenceFree(d_U_adv_diff_var, !d_Q_fcn);
    return;
} // registerAdvDiffHierarchyIntegrator

void
INSHierarchyIntegrator::setStokesSpecifications(StokesSpecifications problem_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_problem_coefs = problem_coefs;
    return;
} // setStokesSpecifications

const StokesSpecifications*
INSHierarchyIntegrator::getStokesSpecifications() const
{
    return &d_problem_coefs;
} // getStokesSpecifications

void
INSHierarchyIntegrator::registerPhysicalBoundaryConditions(const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(bc_coefs.size() == NDIM);
#endif
    d_bc_coefs = bc_coefs;
    return;
} // registerPhysicalBoundaryConditions

const std::vector<RobinBcCoefStrategy<NDIM>*>&
INSHierarchyIntegrator::getVelocityBoundaryConditions() const
{
    return d_U_bc_coefs;
} // getVelocityBoundaryConditions

RobinBcCoefStrategy<NDIM>*
INSHierarchyIntegrator::getPressureBoundaryConditions() const
{
    return d_P_bc_coef;
} // getPressureBoundaryConditions

void
INSHierarchyIntegrator::registerVelocityInitialConditions(Pointer<CartGridFunction> U_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_U_init = U_init;
    return;
} // registerVelocityInitialConditions

void
INSHierarchyIntegrator::registerPressureInitialConditions(Pointer<CartGridFunction> P_init)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_P_init = P_init;
    return;
} // registerPressureInitialConditions

void
INSHierarchyIntegrator::registerBodyForceFunction(Pointer<CartGridFunction> F_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    if (d_F_fcn)
    {
        Pointer<CartGridFunctionSet> p_F_fcn = d_F_fcn;
        if (!p_F_fcn)
        {
            pout << d_object_name << "::registerBodyForceFunction(): WARNING:\n"
                 << "  body force function has already been set.\n"
                 << "  functions will be evaluated in the order in which they were registered "
                    "with "
                    "the solver\n"
                 << "  when evaluating the body force term value.\n";
            p_F_fcn = new CartGridFunctionSet(d_object_name + "::body_force_function_set");
            p_F_fcn->addFunction(d_F_fcn);
        }
        p_F_fcn->addFunction(F_fcn);
    }
    else
    {
        d_F_fcn = F_fcn;
    }
    return;
} // registerBodyForceFunction

void
INSHierarchyIntegrator::registerFluidSourceFunction(Pointer<CartGridFunction> Q_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    if (d_Q_fcn)
    {
        Pointer<CartGridFunctionSet> p_Q_fcn = d_Q_fcn;
        if (!p_Q_fcn)
        {
            pout << d_object_name << "::registerFluidSourceFunction(): WARNING:\n"
                 << "  fluid source function has already been set.\n"
                 << "  functions will be evaluated in the order in which they were registered "
                    "with "
                    "the solver\n"
                 << "  when evaluating the fluid source term value.\n";
            p_Q_fcn = new CartGridFunctionSet(d_object_name + "::fluid_source_function_set");
            p_Q_fcn->addFunction(d_Q_fcn);
        }
        p_Q_fcn->addFunction(Q_fcn);
    }
    else
    {
        d_Q_fcn = Q_fcn;
    }
    return;
} // registerFluidSourceFunction

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getVelocityVariable() const
{
    return d_U_var;
} // getVelocityVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getPressureVariable() const
{
    return d_P_var;
} // getPressureVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getBodyForceVariable() const
{
    return d_F_var;
} // getBodyForceVariable

Pointer<Variable<NDIM> >
INSHierarchyIntegrator::getFluidSourceVariable() const
{
    return d_Q_var;
} // getFluidSourceVariable

Pointer<FaceVariable<NDIM, double> >
INSHierarchyIntegrator::getAdvectionVelocityVariable() const
{
    return d_U_adv_diff_var;
} // getAdvectionVelocityVariable

std::vector<RobinBcCoefStrategy<NDIM>*>
INSHierarchyIntegrator::getIntermediateVelocityBoundaryConditions() const
{
    return d_U_star_bc_coefs;
} // getIntermediateVelocityBoundaryConditions

RobinBcCoefStrategy<NDIM>*
INSHierarchyIntegrator::getProjectionBoundaryConditions() const
{
    return d_Phi_bc_coef;
} // getProjectionBoundaryConditions

void
INSHierarchyIntegrator::registerMassDensityVariable(Pointer<Variable<NDIM> > rho_var)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_rho_var);
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_rho_var = rho_var;
    return;
} // registerMassDensityVariable

void
INSHierarchyIntegrator::setMassDensityFunction(Pointer<CartGridFunction> rho_fcn)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_rho_fcn = rho_fcn;
    return;
} // registerMassDensityFunction

Pointer<CartGridFunction>
INSHierarchyIntegrator::getMassDensityFunction() const
{
    return d_rho_fcn;
} // getMassDensityFunction

void
INSHierarchyIntegrator::setCreepingFlow(bool creeping_flow)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
#endif
    d_creeping_flow = creeping_flow;
    d_convective_op.setNull();
    d_convective_difference_form = UNKNOWN_CONVECTIVE_DIFFERENCING_TYPE;
    return;
} // setCreepingFlow

bool
INSHierarchyIntegrator::getCreepingFlow() const
{
    return d_creeping_flow;
} // getCreepingFlow

void
INSHierarchyIntegrator::setConvectiveOperatorType(const std::string& op_type)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(!d_convective_op);
    TBOX_ASSERT(!d_creeping_flow);
#endif
    d_convective_op_type = op_type;
    return;
} // setConvectiveOperatorType

const std::string&
INSHierarchyIntegrator::getConvectiveOperatorType() const
{
    return d_convective_op_type;
} // getConvectiveOperatorType

void
INSHierarchyIntegrator::setConvectiveDifferencingType(ConvectiveDifferencingType difference_form)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(!d_convective_op);
    TBOX_ASSERT(!d_creeping_flow);
#endif
    d_convective_difference_form = difference_form;
    return;
} // setConvectiveDifferencingType

ConvectiveDifferencingType
INSHierarchyIntegrator::getConvectiveDifferencingType() const
{
    return d_convective_difference_form;
} // getConvectiveDifferencingType

void
INSHierarchyIntegrator::setConvectiveOperator(Pointer<ConvectiveOperator> convective_op)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(!d_convective_op);
#endif
    d_convective_op = convective_op;
    d_creeping_flow = !d_convective_op;
    return;
} // setConvectiveOperator

void
INSHierarchyIntegrator::setConvectiveOperatorNeedsInit()
{
    d_convective_op_needs_init = true;
    return;
}

void
INSHierarchyIntegrator::setVelocitySubdomainSolver(Pointer<PoissonSolver> velocity_solver)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(!d_velocity_solver);
#endif
    d_velocity_solver = velocity_solver;
    return;
} // setVelocitySubdomainSolver

void
INSHierarchyIntegrator::setVelocitySubdomainSolverNeedsInit()
{
    d_velocity_solver_needs_init = true;
    return;
}

void
INSHierarchyIntegrator::setPressureSubdomainSolver(Pointer<PoissonSolver> pressure_solver)
{
#if !defined(NDEBUG)
    TBOX_ASSERT(!d_integrator_is_initialized);
    TBOX_ASSERT(!d_pressure_solver);
#endif
    d_pressure_solver = pressure_solver;
    return;
} // setPressureSubdomainSolver

void
INSHierarchyIntegrator::setPressureSubdomainSolverNeedsInit()
{
    d_pressure_solver_needs_init = true;
    return;
}

int
INSHierarchyIntegrator::getNumberOfCycles() const
{
    int num_cycles = d_num_cycles;
    if (!d_creeping_flow && MathUtilities<double>::equalEps(d_integrator_time, d_start_time) &&
        is_multistep_time_stepping_type(d_convective_time_stepping_type) &&
        d_init_convective_time_stepping_type != FORWARD_EULER)
    {
        num_cycles = std::max(2, num_cycles);
    }
    return num_cycles;
} // getNumberOfCycles

/////////////////////////////// PROTECTED ////////////////////////////////////

INSHierarchyIntegrator::INSHierarchyIntegrator(const std::string& object_name,
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
      d_default_bc_coefs(d_object_name + "::default_bc_coefs", Pointer<Database>(NULL)),
      d_bc_coefs(NDIM, static_cast<RobinBcCoefStrategy<NDIM>*>(NULL)),
      d_traction_bc_type(TRACTION),
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
    d_normalize_velocity = false;
    d_convective_op_type = "DEFAULT";
    d_convective_difference_form = ADVECTIVE;
    d_convective_op_input_db = new MemoryDatabase(d_object_name + "::convective_op_input_db");
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
        d_default_bc_coefs.setBoundaryValue(2 * d, 0.0);
        d_default_bc_coefs.setBoundaryValue(2 * d + 1, 0.0);
    }
    registerPhysicalBoundaryConditions(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, &d_default_bc_coefs));

    // Setup physical boundary conditions objects.
    d_U_star_bc_coefs.resize(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        d_U_star_bc_coefs[d] = new INSIntermediateVelocityBcCoef(d, d_bc_coefs);
    }
    d_Phi_bc_coef = new INSProjectionBcCoef(d_bc_coefs);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    if (input_db) getFromInput(input_db, from_restart);

    // Initialize an advection velocity variable.  NOTE: Patch data are
    // allocated for this variable only when an advection-diffusion solver is
    // registered with the INSHierarchyIntegrator.
    d_U_adv_diff_var = new FaceVariable<NDIM, double>(d_object_name + "::U_adv_diff");
    return;
} // INSHierarchyIntegrator

double
INSHierarchyIntegrator::getMaximumTimeStepSizeSpecialized()
{
    double dt = HierarchyIntegrator::getMaximumTimeStepSizeSpecialized();
    for (int ln = 0; ln <= d_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        dt = std::min(dt, d_cfl_max * getStableTimestep(level));
    }
    return dt;
} // getMaximumTimeStepSizeSpecialized

double
INSHierarchyIntegrator::getStableTimestep(Pointer<PatchLevel<NDIM> > level) const
{
    double stable_dt = std::numeric_limits<double>::max();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM> > patch = level->getPatch(p());
        stable_dt = std::min(stable_dt, getStableTimestep(patch));
    }
    stable_dt = SAMRAI_MPI::minReduction(stable_dt);
    return stable_dt;
} // getStableTimestep

void
INSHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    db->putInteger("INS_HIERARCHY_INTEGRATOR_VERSION", INS_HIERARCHY_INTEGRATOR_VERSION);
    db->putString("d_viscous_time_stepping_type", enum_to_string<TimeSteppingType>(d_viscous_time_stepping_type));
    db->putString("d_convective_time_stepping_type", enum_to_string<TimeSteppingType>(d_convective_time_stepping_type));
    db->putString("d_init_convective_time_stepping_type",
                  enum_to_string<TimeSteppingType>(d_init_convective_time_stepping_type));
    db->putDouble("d_rho", d_problem_coefs.getRho());
    db->putDouble("d_mu", d_problem_coefs.getMu());
    db->putDouble("d_lambda", d_problem_coefs.getLambda());
    db->putDouble("d_cfl_max", d_cfl_max);
    db->putBool("d_using_vorticity_tagging", d_using_vorticity_tagging);
    if (d_Omega_rel_thresh.size() > 0) db->putDoubleArray("d_Omega_rel_thresh", d_Omega_rel_thresh);
    if (d_Omega_abs_thresh.size() > 0) db->putDoubleArray("d_Omega_abs_thresh", d_Omega_abs_thresh);
    db->putDouble("d_Omega_max", d_Omega_max);
    db->putBool("d_normalize_pressure", d_normalize_pressure);
    db->putBool("d_normalize_velocity", d_normalize_velocity);
    db->putString("d_convective_op_type", d_convective_op_type);
    db->putString("d_convective_difference_form",
                  enum_to_string<ConvectiveDifferencingType>(d_convective_difference_form));
    db->putBool("d_creeping_flow", d_creeping_flow);
    db->putDouble("d_regrid_max_div_growth_factor", d_regrid_max_div_growth_factor);
    db->putDouble("d_U_scale", d_U_scale);
    db->putDouble("d_P_scale", d_P_scale);
    db->putDouble("d_F_scale", d_F_scale);
    db->putDouble("d_Q_scale", d_Q_scale);
    db->putDouble("d_Omega_scale", d_Omega_scale);
    db->putDouble("d_Div_U_scale", d_Div_U_scale);
    db->putBool("d_output_U", d_output_U);
    db->putBool("d_output_P", d_output_P);
    db->putBool("d_output_F", d_output_F);
    db->putBool("d_output_Q", d_output_Q);
    db->putBool("d_output_Omega", d_output_Omega);
    db->putBool("d_output_Div_U", d_output_Div_U);
    return;
} // putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
INSHierarchyIntegrator::getFromInput(Pointer<Database> db, const bool is_from_restart)
{
    if (!is_from_restart)
    {
        if (db->keyExists("viscous_time_stepping_type"))
            d_viscous_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("viscous_time_stepping_type"));
        else if (db->keyExists("viscous_timestepping_type"))
            d_viscous_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("viscous_timestepping_type"));
        if (db->keyExists("convective_time_stepping_type"))
            d_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("convective_time_stepping_type"));
        else if (db->keyExists("convective_timestepping_type"))
            d_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("convective_timestepping_type"));
        if (db->keyExists("init_convective_time_stepping_type"))
            d_init_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("init_convective_time_stepping_type"));
        else if (db->keyExists("init_convective_timestepping_type"))
            d_init_convective_time_stepping_type =
                string_to_enum<TimeSteppingType>(db->getString("init_convective_timestepping_type"));
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
    if (db->keyExists("cfl"))
        d_cfl_max = db->getDouble("cfl");
    else if (db->keyExists("cfl_max"))
        d_cfl_max = db->getDouble("cfl_max");
    else if (db->keyExists("CFL"))
        d_cfl_max = db->getDouble("CFL");
    else if (db->keyExists("CFL_max"))
        d_cfl_max = db->getDouble("CFL_max");
    if (db->keyExists("using_vorticity_tagging")) d_using_vorticity_tagging = db->getBool("using_vorticity_tagging");
    if (db->keyExists("Omega_rel_thresh"))
        d_Omega_rel_thresh = db->getDoubleArray("Omega_rel_thresh");
    else if (db->keyExists("omega_rel_thresh"))
        d_Omega_rel_thresh = db->getDoubleArray("omega_rel_thresh");
    else if (db->keyExists("vorticity_rel_thresh"))
        d_Omega_rel_thresh = db->getDoubleArray("vorticity_rel_thresh");
    if (db->keyExists("Omega_abs_thresh"))
        d_Omega_abs_thresh = db->getDoubleArray("Omega_abs_thresh");
    else if (db->keyExists("omega_abs_thresh"))
        d_Omega_abs_thresh = db->getDoubleArray("omega_abs_thresh");
    else if (db->keyExists("vorticity_abs_thresh"))
        d_Omega_abs_thresh = db->getDoubleArray("vorticity_abs_thresh");
    if (db->keyExists("normalize_pressure")) d_normalize_pressure = db->getBool("normalize_pressure");
    if (db->keyExists("normalize_velocity")) d_normalize_velocity = db->getBool("normalize_velocity");
    if (db->keyExists("convective_op_type"))
        d_convective_op_type = db->getString("convective_op_type");
    else if (db->keyExists("convective_operator_type"))
        d_convective_op_type = db->getString("convective_operator_type");
    else if (db->keyExists("default_convective_op_type"))
        d_convective_op_type = db->getString("default_convective_op_type");
    else if (db->keyExists("default_convective_operator_type"))
        d_convective_op_type = db->getString("default_convective_operator_type");
    if (db->keyExists("convective_difference_form"))
        d_convective_difference_form =
            string_to_enum<ConvectiveDifferencingType>(db->getString("convective_difference_form"));
    else if (db->keyExists("convective_difference_type"))
        d_convective_difference_form =
            string_to_enum<ConvectiveDifferencingType>(db->getString("convective_difference_type"));
    else if (db->keyExists("default_convective_difference_form"))
        d_convective_difference_form =
            string_to_enum<ConvectiveDifferencingType>(db->getString("default_convective_difference_form"));
    else if (db->keyExists("default_convective_difference_type"))
        d_convective_difference_form =
            string_to_enum<ConvectiveDifferencingType>(db->getString("default_convective_difference_type"));
    if (db->keyExists("convective_op_db"))
        d_convective_op_input_db = db->getDatabase("convective_op_db");
    else if (db->keyExists("default_convective_op_db"))
        d_convective_op_input_db = db->getDatabase("default_convective_op_db");
    if (db->keyExists("creeping_flow")) d_creeping_flow = db->getBool("creeping_flow");
    if (db->keyExists("regrid_max_div_growth_factor"))
        d_regrid_max_div_growth_factor = db->getDouble("regrid_max_div_growth_factor");
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
    if (db->keyExists("traction_bc_type"))
        d_traction_bc_type = string_to_enum<TractionBcType>(db->getString("traction_bc_type"));

    if (db->keyExists("velocity_solver_type"))
    {
        d_velocity_solver_type = db->getString("velocity_solver_type");
        if (db->keyExists("velocity_solver_db")) d_velocity_solver_db = db->getDatabase("velocity_solver_db");
    }
    if (!d_velocity_solver_db) d_velocity_solver_db = new MemoryDatabase("velocity_solver_db");

    if (db->keyExists("velocity_precond_type"))
    {
        d_velocity_precond_type = db->getString("velocity_precond_type");
        if (db->keyExists("velocity_precond_db")) d_velocity_precond_db = db->getDatabase("velocity_precond_db");
    }
    if (!d_velocity_precond_db) d_velocity_precond_db = new MemoryDatabase("velocity_precond_db");

    if (db->keyExists("pressure_solver_type"))
    {
        d_pressure_solver_type = db->getString("pressure_solver_type");
        if (db->keyExists("pressure_solver_db")) d_pressure_solver_db = db->getDatabase("pressure_solver_db");
    }
    if (!d_pressure_solver_db) d_pressure_solver_db = new MemoryDatabase("pressure_solver_db");

    if (db->keyExists("pressure_precond_type"))
    {
        d_pressure_precond_type = db->getString("pressure_precond_type");
        if (db->keyExists("pressure_precond_db")) d_pressure_precond_db = db->getDatabase("pressure_precond_db");
    }
    if (!d_pressure_precond_db) d_pressure_precond_db = new MemoryDatabase("pressure_precond_db");

    if (db->keyExists("regrid_projection_solver_type"))
    {
        d_regrid_projection_solver_type = db->getString("regrid_projection_solver_type");
        if (db->keyExists("regrid_projection_solver_db"))
            d_regrid_projection_solver_db = db->getDatabase("regrid_projection_solver_db");
    }
    if (!d_regrid_projection_solver_db)
        d_regrid_projection_solver_db = new MemoryDatabase("regrid_projection_solver_db");

    if (db->keyExists("regrid_projection_precond_type"))
    {
        d_regrid_projection_precond_type = db->getString("regrid_projection_precond_type");
        if (db->keyExists("regrid_projection_precond_db"))
            d_regrid_projection_precond_db = db->getDatabase("regrid_projection_precond_db");
    }
    if (!d_regrid_projection_precond_db)
        d_regrid_projection_precond_db = new MemoryDatabase("regrid_projection_precond_db");
    return;
} // getFromInput

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
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file."
                                 << std::endl);
    }
    int ver = db->getInteger("INS_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != INS_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    d_viscous_time_stepping_type = string_to_enum<TimeSteppingType>(db->getString("d_viscous_time_stepping_type"));
    d_convective_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_convective_time_stepping_type"));
    d_init_convective_time_stepping_type =
        string_to_enum<TimeSteppingType>(db->getString("d_init_convective_time_stepping_type"));
    d_problem_coefs.setRho(db->getDouble("d_rho"));
    d_problem_coefs.setMu(db->getDouble("d_mu"));
    d_problem_coefs.setLambda(db->getDouble("d_lambda"));
    d_num_cycles = db->getInteger("d_num_cycles");
    d_cfl_max = db->getDouble("d_cfl_max");
    d_using_vorticity_tagging = db->getBool("d_using_vorticity_tagging");
    if (db->keyExists("d_Omega_rel_thresh"))
        d_Omega_rel_thresh = db->getDoubleArray("d_Omega_rel_thresh");
    else
        d_Omega_rel_thresh.resizeArray(0);
    if (db->keyExists("d_Omega_abs_thresh"))
        d_Omega_abs_thresh = db->getDoubleArray("d_Omega_abs_thresh");
    else
        d_Omega_abs_thresh.resizeArray(0);
    d_Omega_max = db->getDouble("d_Omega_max");
    d_normalize_pressure = db->getBool("d_normalize_pressure");
    d_normalize_velocity = db->getBool("d_normalize_velocity");
    d_convective_op_type = db->getString("d_convective_op_type");
    d_convective_difference_form =
        string_to_enum<ConvectiveDifferencingType>(db->getString("d_convective_difference_form"));
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
} // getFromRestart

//////////////////////////////////////////////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
