// Filename: IBImplicitStaggeredHierarchyIntegrator.C
// Created on 07 Apr 2012 by Boyce Griffith
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

#include "IBImplicitStaggeredHierarchyIntegrator.h"

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
#include <ibamr/IBImplicitStaggeredPETScLevelSolver.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/namespaces.h>

// IBTK INCLUDES
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScMatUtilities.h>

// C++ STDLIB INCLUDES
#include <algorithm>

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBImplicitStaggeredHierarchyIntegrator restart file data.
static const int IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitStaggeredHierarchyIntegrator::IBImplicitStaggeredHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<IBStrategy> ib_method_ops,
    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
    bool register_for_restart)
    : IBHierarchyIntegrator(object_name, input_db, ib_method_ops, ins_hier_integrator, register_for_restart)
{
    // Setup IB ops object to use "fixed" Lagrangian-Eulerian coupling
    // operators.
    d_ib_method_ops->setUseFixedLEOperators(true);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;
}// IBImplicitStaggeredHierarchyIntegrator

IBImplicitStaggeredHierarchyIntegrator::~IBImplicitStaggeredHierarchyIntegrator()
{
    // intentionally blank
    return;
}// ~IBImplicitStaggeredHierarchyIntegrator

void
IBImplicitStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const int num_cycles)
{
    IBHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    d_current_time = current_time;
    d_new_time = new_time;
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    TBOX_ASSERT(d_time_stepping_type == MIDPOINT_RULE);

    // Allocate Eulerian scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_u_idx, current_time);
        level->allocatePatchData(d_f_idx, current_time);
        if (d_ib_method_ops->hasFluidSources())
        {
            level->allocatePatchData(d_p_idx, current_time);
            level->allocatePatchData(d_q_idx, current_time);
        }
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data    ,     new_time);
    }

    // Initialize IB data.
    d_ib_method_ops->preprocessIntegrateData(current_time, new_time, num_cycles);

    // Initialize the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    if (ins_num_cycles != d_current_num_cycles && d_current_num_cycles != 1)
    {
        TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                   << "  attempting to perform " << d_current_num_cycles << " cycles of fixed point iteration.\n"
                   << "  number of cycles required by Navier-Stokes solver = " << ins_num_cycles << ".\n"
                   << "  current implementation requires either that both solvers use the same number of cycles,\n"
                   << "  or that the IB solver use only a single cycle.\n");
    }
    d_ins_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, ins_num_cycles);

    // Compute an initial prediction of the updated positions of the Lagrangian
    // structure.
    //
    // NOTE: The velocity should already have been interpolated to the
    // curvilinear mesh and should not need to be re-interpolated.
    if (d_enable_logging) plog << d_object_name << "::preprocessIntegrateHierarchy(): performing Lagrangian forward Euler step\n";
    d_ib_method_ops->eulerStep(current_time, new_time);

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
}// preprocessIntegrateHierarchy

void
IBImplicitStaggeredHierarchyIntegrator::integrateHierarchy(
    const double current_time,
    const double new_time,
    const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    const double half_time = current_time+0.5*(new_time-current_time);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getCurrentContext());
    const int u_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(), d_ins_hier_integrator->getNewContext());
    const int p_new_idx     = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getPressureVariable(), d_ins_hier_integrator->getNewContext());

    // Fix the positions of the Lagrangian-Eulerian interaction operators.
    d_ib_method_ops->updateFixedLEOperators();
    d_ib_method_ops->getLEOperatorPositions(d_X_LE_vec, d_hierarchy->getFinestLevelNumber(), half_time);

    // Compute the Lagrangian source/sink strengths and spread them to the
    // Eulerian grid.
    if (d_ib_method_ops->hasFluidSources())
    {
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian fluid source strength\n";
        d_ib_method_ops->computeLagrangianFluidSource(half_time);
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian fluid source strength to the Eulerian grid\n";
        d_hier_pressure_data_ops->setToScalar(d_q_idx, 0.0);
        d_ib_method_ops->spreadFluidSource(d_q_idx, getProlongRefineSchedules(d_object_name+"::q"), half_time);
    }

    // Solve the incompressible Navier-Stokes equations.
    d_ib_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): solving the modified incompressible Navier-Stokes equations\n";
    if (d_current_num_cycles > 1)
    {
        d_ins_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);
    }
    else
    {
        const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
        for (int cycle = 0; cycle < ins_num_cycles; ++cycle)
        {
            d_ins_hier_integrator->integrateHierarchy(current_time, new_time, cycle_num);
        }
    }
    d_ib_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
    if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): interpolating Eulerian velocity to the Lagrangian mesh\n";
    d_ib_method_ops->interpolateVelocity(d_u_idx, getCoarsenSchedules(d_object_name+"::u::CONSERVATIVE_COARSEN"), getGhostfillRefineSchedules(d_object_name+"::u"), half_time);

    // Compute an updated prediction of the updated positions of the Lagrangian
    // structure.
    d_ib_method_ops->midpointStep(current_time, new_time);

    // Compute the pressure at the updated locations of any distributed internal
    // fluid sources or sinks.
    if (d_ib_method_ops->hasFluidSources())
    {
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): interpolating Eulerian fluid pressure to the Lagrangian mesh\n";
        d_hier_pressure_data_ops->copyData(d_p_idx, p_new_idx);
        d_ib_method_ops->interpolatePressure(d_p_idx, getCoarsenSchedules(d_object_name+"::p::CONSERVATIVE_COARSEN"), getGhostfillRefineSchedules(d_object_name+"::p"), half_time);
    }

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
}// integrateHierarchy

void
IBImplicitStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(
    const double current_time,
    const double new_time,
    const bool skip_synchronize_new_state_data,
    const int num_cycles)
{
    IBHierarchyIntegrator::postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time-current_time;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
    if (d_enable_logging) plog << d_object_name << "::postprocessIntegrateHierarchy(): interpolating Eulerian velocity to the Lagrangian mesh\n";
    d_ib_method_ops->interpolateVelocity(d_u_idx, getCoarsenSchedules(d_object_name+"::u::CONSERVATIVE_COARSEN"), getGhostfillRefineSchedules(d_object_name+"::u"), new_time);

    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_enable_logging) plog << d_object_name << "::postprocessIntegrateHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(NEW_DATA);
    }

    // Determine the CFL number.
    double cfl_max = 0.0;
    PatchCellDataOpsReal<NDIM,double> patch_cc_ops;
    PatchSideDataOpsReal<NDIM,double> patch_sc_ops;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double dx_min = *(std::min_element(dx,dx+NDIM));
            Pointer<CellData<NDIM,double> > u_cc_new_data = patch->getPatchData(u_new_idx);
            Pointer<SideData<NDIM,double> > u_sc_new_data = patch->getPatchData(u_new_idx);
            double u_max = 0.0;
            if (u_cc_new_data) u_max = patch_cc_ops.maxNorm(u_cc_new_data, patch_box);
            if (u_sc_new_data) u_max = patch_sc_ops.maxNorm(u_sc_new_data, patch_box);
            cfl_max = std::max(cfl_max, u_max*dt/dx_min);
        }
    }
    cfl_max = SAMRAI_MPI::maxReduction(cfl_max);
    d_regrid_cfl_estimate += cfl_max;
    if (d_enable_logging) plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl_max << "\n";
    if (d_enable_logging) plog << d_object_name << "::postprocessIntegrateHierarchy(): estimated upper bound on IB point displacement since last regrid = " << d_regrid_cfl_estimate << "\n";

    // Deallocate the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    d_ins_hier_integrator->postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data, ins_num_cycles);

    // Deallocate IB data.
    d_ib_method_ops->postprocessIntegrateData(current_time, new_time, num_cycles);

    // Deallocate Eulerian scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_u_idx);
        level->deallocatePatchData(d_f_idx);
        if (d_ib_method_ops->hasFluidSources())
        {
            level->deallocatePatchData(d_p_idx);
            level->deallocatePatchData(d_q_idx);
        }
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data    );
    }

    d_current_time = std::numeric_limits<double>::quiet_NaN();
    d_new_time = std::numeric_limits<double>::quiet_NaN();

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
}// postprocessIntegrateHierarchy

void
IBImplicitStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;
#if 0
    // Setup the fluid solver for implicit coupling.
    if (d_body_force_fcn)
    {
        d_ins_hier_integrator->registerBodyForceFunction(d_body_force_fcn);
    }
    const std::string stokes_prefix = "stokes_";
    Pointer<INSStaggeredHierarchyIntegrator> p_ins_hier_integrator = d_ins_hier_integrator;
    const StokesSpecifications* const problem_coefs = p_ins_hier_integrator->getStokesSpecifications();
    const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs = p_ins_hier_integrator->getVelocityBoundaryConditions();
    RobinBcCoefStrategy<NDIM>* const P_bc_coef = p_ins_hier_integrator->getPressureBoundaryConditions();
    TBOX_ERROR("need to add U bc helper . . . ?\n");
//  d_stokes_op = new StaggeredStokesOperator(d_object_name+"::StaggeredStokesOperator", problem_coefs, TRAPEZOIDAL_RULE, U_bc_coefs, NULL, P_bc_coef, buildHierarchyMathOps(hierarchy));
    Pointer<NewtonKrylovSolver> modified_stokes_solver = NULL; // new PETScNewtonKrylovSolver(d_object_name+"::stokes_solver", stokes_prefix);
    d_F_op = new IBImplicitStaggeredHierarchyIntegrator::Operator(this);
    modified_stokes_solver->setOperator(d_F_op);
    d_J_op = new IBImplicitStaggeredHierarchyIntegrator::Jacobian(this);
    modified_stokes_solver->setJacobian(d_J_op);
//  p_ins_hier_integrator->setStokesSolver(d_stokes_op, modified_stokes_solver);
//  d_modified_stokes_pc = new IBImplicitStaggeredPETScLevelSolver(d_object_name+"::PETScLevelSolver", *problem_coefs, &d_J_mat, PETScMatUtilities::ib_4_interp_fcn, PETScMatUtilities::ib_4_interp_stencil, &d_X_LE_vec, U_bc_coefs);
//  modified_stokes_solver->getLinearSolver()->setPreconditioner(d_modified_stokes_pc);
#endif
    TBOX_ERROR("not currently implemented!\n"); // XXXX

    // Finish initializing the hierarchy integrator.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);
}// initializeHierarchyIntegrator

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBImplicitStaggeredHierarchyIntegrator::putToDatabaseSpecialized(
    Pointer<Database> db)
{
    IBHierarchyIntegrator::putToDatabaseSpecialized(db);
    db->putInteger("IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION",IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION);
    return;
}// putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBImplicitStaggeredHierarchyIntegrator::getFromRestart()
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
    int ver = db->getInteger("IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
}// getFromRestart

IBImplicitStaggeredHierarchyIntegrator::Operator::Operator(
    const IBImplicitStaggeredHierarchyIntegrator* ib_solver)
    : GeneralOperator(ib_solver->d_object_name+"::Operator"),
      d_ib_solver(ib_solver)
{
    // intentionally blank
    return;
}// Operator

IBImplicitStaggeredHierarchyIntegrator::Operator::~Operator()
{
    // intentionally blank
    return;
}// ~Operator();

void
IBImplicitStaggeredHierarchyIntegrator::Operator::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    const double current_time = d_ib_solver->d_current_time;
    const double new_time = d_ib_solver->d_new_time;
    const double half_time = current_time+0.5*(new_time-current_time);
    IBStrategy* ib_method_ops = d_ib_solver->d_ib_method_ops;
    Pointer<HierarchyDataOpsReal<NDIM,double> > hier_velocity_data_ops = d_ib_solver->d_hier_velocity_data_ops;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<SideVariable<NDIM,double> > u_current_var = d_ib_solver->d_ins_hier_integrator->getVelocityVariable();
    Pointer<VariableContext> u_current_ctx = d_ib_solver->d_ins_hier_integrator->getCurrentContext();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(u_current_var, u_current_ctx);
    const int u_new_idx = x.getComponentDescriptorIndex(0);
    const int u_half_ib_idx = d_ib_solver->d_u_idx;
    hier_velocity_data_ops->linearSum(u_half_ib_idx, 0.5, u_current_idx, 0.5, u_new_idx);

    const int f_half_idx = y.getComponentDescriptorIndex(0);

    // Compute the "fluid" part of the operator.
    Pointer<StaggeredStokesOperator> stokes_op = d_ib_solver->d_stokes_op;
    stokes_op->setHomogeneousBc(false);
    stokes_op->setTimeInterval(current_time, new_time);
    stokes_op->apply(x, y);

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_coarsen_scheds = d_ib_solver->getCoarsenSchedules(d_ib_solver->d_object_name+"::u::CONSERVATIVE_COARSEN");
    const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghostfill_scheds = d_ib_solver->getGhostfillRefineSchedules(d_ib_solver->d_object_name+"::u");
    ib_method_ops->interpolateVelocity(u_half_ib_idx, u_coarsen_scheds, u_ghostfill_scheds, half_time);

    // Compute an updated prediction of the updated positions of the Lagrangian
    // structure.
    ib_method_ops->midpointStep(current_time, new_time);

    // Compute the "structure" part of the operator.
    ib_method_ops->computeLagrangianForce(half_time);
    hier_velocity_data_ops->scale(f_half_idx, -1.0, f_half_idx);
    const std::vector<Pointer<RefineSchedule<NDIM> > > f_prolong_scheds(u_ghostfill_scheds.size(), Pointer<RefineSchedule<NDIM> >(NULL));
    ib_method_ops->spreadForce(f_half_idx, f_prolong_scheds, half_time);
    hier_velocity_data_ops->scale(f_half_idx, -1.0, f_half_idx);
    return;
}// apply

void
IBImplicitStaggeredHierarchyIntegrator::Operator::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& in,
    const SAMRAIVectorReal<NDIM,double>& out)
{
    d_ib_solver->d_stokes_op->initializeOperatorState(in,out);
    return;
}// initializeOperatorState

void
IBImplicitStaggeredHierarchyIntegrator::Operator::deallocateOperatorState()
{
    d_ib_solver->d_stokes_op->deallocateOperatorState();
    return;
}// deallocateOperatorState

IBImplicitStaggeredHierarchyIntegrator::Jacobian::Jacobian(
    IBImplicitStaggeredHierarchyIntegrator* ib_solver)
    : JacobianOperator(ib_solver->getName()+"::Jacobian"),
      d_ib_solver(ib_solver),
      d_J_mat(NULL),
      d_J_is_set(false),
      d_x_base(NULL)
{
    // intentionally blank
    return;
}// Jacobian

IBImplicitStaggeredHierarchyIntegrator::Jacobian::~Jacobian()
{
    deallocateOperatorState();
    return;
}// ~Jacobian();

void
IBImplicitStaggeredHierarchyIntegrator::Jacobian::formJacobian(
    SAMRAIVectorReal<NDIM,double>& x)
{
    const double current_time = d_ib_solver->d_current_time;
    const double new_time = d_ib_solver->d_new_time;
    const double half_time = current_time+0.5*(new_time-current_time);
    const double dt = new_time-current_time;
    IBStrategy* ib_method_ops = d_ib_solver->d_ib_method_ops;

    d_x_base = Pointer<SAMRAIVectorReal<NDIM,double> >(&x, false);
    if (d_J_is_set)
    {
        int ierr = MatZeroEntries(d_J_mat); IBTK_CHKERRQ(ierr);
    }
    ib_method_ops->computeLagrangianForceJacobian(d_J_mat, MAT_FINAL_ASSEMBLY, /* X_coef */ -0.25*dt, /* U_coef */ -0.5, half_time);
    Pointer<IBImplicitStaggeredPETScLevelSolver> p_modified_stokes_pc = d_ib_solver->d_modified_stokes_pc;
    p_modified_stokes_pc->initializeOperator();
    d_J_is_set = true;
    return;
}// formJacobian

Pointer<SAMRAIVectorReal<NDIM,double> >
IBImplicitStaggeredHierarchyIntegrator::Jacobian::getBaseVector() const
{
    return d_x_base;
}// getBaseVector

void
IBImplicitStaggeredHierarchyIntegrator::Jacobian::apply(
    SAMRAIVectorReal<NDIM,double>& x,
    SAMRAIVectorReal<NDIM,double>& y)
{
    const double current_time = d_ib_solver->d_current_time;
    const double new_time = d_ib_solver->d_new_time;
    const double half_time = current_time+0.5*(new_time-current_time);
    IBStrategy* ib_method_ops = d_ib_solver->d_ib_method_ops;
    Pointer<HierarchyDataOpsReal<NDIM,double> > hier_velocity_data_ops = d_ib_solver->d_hier_velocity_data_ops;

    const int u_new_idx = x.getComponentDescriptorIndex(0);
    const int u_new_ib_idx = d_ib_solver->d_u_idx;
    hier_velocity_data_ops->copyData(u_new_ib_idx, u_new_idx);

    const int f_half_idx = y.getComponentDescriptorIndex(0);

    // Compute the "fluid" part of the operator.
    Pointer<StaggeredStokesOperator> stokes_op = d_ib_solver->d_stokes_op;
    stokes_op->setHomogeneousBc(true);
    stokes_op->setTimeInterval(current_time, new_time);
    stokes_op->apply(x, y);

    // Compute the "structure" part of the operator.
    const std::vector<Pointer<CoarsenSchedule<NDIM> > >& u_coarsen_scheds = d_ib_solver->getCoarsenSchedules(d_ib_solver->d_object_name+"::u::CONSERVATIVE_COARSEN");
    const std::vector<Pointer<RefineSchedule<NDIM> > >& u_ghostfill_scheds = d_ib_solver->getGhostfillRefineSchedules(d_ib_solver->d_object_name+"::u");
    const std::vector<Pointer<RefineSchedule<NDIM> > > f_prolong_scheds(u_ghostfill_scheds.size(), Pointer<RefineSchedule<NDIM> >(NULL));
    ib_method_ops->applyLagrangianForceJacobian(f_half_idx, f_prolong_scheds, u_new_ib_idx, u_coarsen_scheds, u_ghostfill_scheds, half_time, d_J_mat);
    return;
}// apply

void
IBImplicitStaggeredHierarchyIntegrator::Jacobian::initializeOperatorState(
    const SAMRAIVectorReal<NDIM,double>& /*in*/,
    const SAMRAIVectorReal<NDIM,double>& /*out*/)
{
    deallocateOperatorState();
    IBStrategy* ib_method_ops = d_ib_solver->d_ib_method_ops;
    std::vector<int> d_nnz, o_nnz;
    ib_method_ops->computeLagrangianForceJacobianNonzeroStructure(d_nnz, o_nnz);
    const int n_local = d_nnz.size();
    int ierr;
    ierr = MatCreateAIJ(PETSC_COMM_WORLD, n_local, n_local, PETSC_DETERMINE, PETSC_DETERMINE, 0, (n_local == 0 ? NULL : &d_nnz[0]), 0, (n_local == 0 ? NULL : &o_nnz[0]), &d_J_mat); IBTK_CHKERRQ(ierr);
    ierr = MatSetBlockSize(d_J_mat, NDIM);
    d_ib_solver->d_J_mat = d_J_mat;
    d_J_is_set = false;
    return;
}// initializeOperatorState

void
IBImplicitStaggeredHierarchyIntegrator::Jacobian::deallocateOperatorState()
{
    if (d_J_mat)
    {
        int ierr = MatDestroy(&d_J_mat); IBTK_CHKERRQ(ierr);
        d_J_mat = NULL;
    }
    return;
}// deallocateOperatorState

void
IBImplicitStaggeredHierarchyIntegrator::Jacobian::enableLogging(
    bool /*enabled*/)
{
    // intentionally blank
    return;
}// enableLogging

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
