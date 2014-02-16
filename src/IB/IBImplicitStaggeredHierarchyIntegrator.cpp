// Filename: IBImplicitStaggeredHierarchyIntegrator.cpp
// Created on 07 Apr 2012 by Boyce Griffith
//
// Copyright (c) 2002-2013, Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <algorithm>
#include <iosfwd>
#include <limits>
#include <memory>
#include <ostream>
#include <vector>

#include "ibtk/KrylovLinearSolver.h"
#include "petscsnes.h"
#include "ibtk/PETScMultiVec.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "GriddingAlgorithm.h"
#include "HierarchyDataOpsReal.h"
#include "IBImplicitStaggeredHierarchyIntegrator.h"
#include "Patch.h"
#include "PatchCellDataOpsReal.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataOpsReal.h"
#include "RefineSchedule.h"
#include "SideData.h"
#include "SideVariable.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/IBImplicitStaggeredPETScLevelSolver.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesOperator.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/ibtk_enums.h"
#include "petscsys.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

namespace IBAMR {
class INSStaggeredHierarchyIntegrator;
}  // namespace IBAMR
namespace SAMRAI {
namespace hier {
template <int DIM> class Box;
}  // namespace hier
namespace xfer {
template <int DIM> class CoarsenSchedule;
}  // namespace xfer
}  // namespace SAMRAI

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

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    TBOX_ASSERT(d_time_stepping_type == MIDPOINT_RULE);

    // Allocate Eulerian scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_u_idx, current_time);
        level->allocatePatchData(d_f_idx, current_time);
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

    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator = d_ins_hier_integrator;
    TBOX_ASSERT(ins_hier_integrator);

    PetscErrorCode ierr;
    
    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double half_time = current_time+0.5*(new_time-current_time);
    
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = ins_hier_integrator->getCurrentContext();
    Pointer<VariableContext> scratch_ctx = ins_hier_integrator->getScratchContext();
    Pointer<VariableContext> new_ctx     = ins_hier_integrator->getNewContext();

    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();  // XXXX: This should be made more generic.
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();
    
    Pointer<Variable<NDIM> > u_var = ins_hier_integrator->getVelocityVariable();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, current_ctx);
    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, scratch_ctx);
    const int u_new_idx     = var_db->mapVariableAndContextToIndex(u_var, new_ctx);

    Pointer<Variable<NDIM> > p_var = ins_hier_integrator->getPressureVariable();
//  const int p_current_idx = var_db->mapVariableAndContextToIndex(p_var, current_ctx);
    const int p_scratch_idx = var_db->mapVariableAndContextToIndex(p_var, scratch_ctx);
//  const int p_new_idx     = var_db->mapVariableAndContextToIndex(p_var, new_ctx);

    // Skip all cycles in the INS solver --- we advance the state data here.
    ins_hier_integrator->skipCycle(current_time, new_time, cycle_num);

    // Setup Eulerian vectors used in solving the implicit IB equations.
    Pointer<SAMRAIVectorReal<NDIM,double> > eul_sol_vec = new SAMRAIVectorReal<NDIM,double>(d_object_name+"::eulerian_sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    eul_sol_vec->addComponent(u_var,u_scratch_idx,wgt_sc_idx,d_hier_velocity_data_ops);
    eul_sol_vec->addComponent(p_var,p_scratch_idx,wgt_cc_idx,d_hier_pressure_data_ops);

    Pointer<SAMRAIVectorReal<NDIM,double> > eul_rhs_vec = eul_sol_vec->cloneVector(d_object_name+"::eulerian_rhs_vec");
    eul_rhs_vec->allocateVectorData(current_time);

    ins_hier_integrator->setupSolverVectors(eul_sol_vec, eul_rhs_vec, current_time, new_time, cycle_num);

    d_stokes_solver = ins_hier_integrator->getStokesSolver();
    Pointer<KrylovLinearSolver> p_stokes_solver = d_stokes_solver;
    TBOX_ASSERT(p_stokes_solver);
    d_stokes_op = p_stokes_solver->getOperator();
    TBOX_ASSERT(d_stokes_op);

    // Setup Lagrangian vectors used in solving the implicit IB equations.
    Vec lag_sol_petsc_vec, lag_rhs_petsc_vec;
    d_ib_method_ops->createSolverVecs(lag_sol_petsc_vec, lag_rhs_petsc_vec, d_hierarchy->getFinestLevelNumber());
    d_ib_method_ops->setupSolverVecs(lag_sol_petsc_vec, lag_rhs_petsc_vec, d_hierarchy->getFinestLevelNumber());
    
    d_ib_method_ops->updateFixedLEOperators();

    // Setup multi-vec objects to store the composite solution and
    // right-hand-side vectors.
    Vec eul_sol_petsc_vec = PETScSAMRAIVectorReal::createPETScVector(eul_sol_vec, PETSC_COMM_WORLD);
    Vec eul_rhs_petsc_vec = PETScSAMRAIVectorReal::createPETScVector(eul_rhs_vec, PETSC_COMM_WORLD);

    Vec sol_petsc_vecs[] = { eul_sol_petsc_vec , lag_sol_petsc_vec };
    Vec rhs_petsc_vecs[] = { eul_rhs_petsc_vec , lag_rhs_petsc_vec };

    Vec composite_sol_petsc_vec, composite_rhs_petsc_vec, composite_res_petsc_vec;
    ierr = VecCreateMultiVec(PETSC_COMM_WORLD, 2, sol_petsc_vecs, &composite_sol_petsc_vec); IBTK_CHKERRQ(ierr);
    ierr = VecCreateMultiVec(PETSC_COMM_WORLD, 2, rhs_petsc_vecs, &composite_rhs_petsc_vec); IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(composite_rhs_petsc_vec, &composite_res_petsc_vec); IBTK_CHKERRQ(ierr);
    
    // Solve the implicit IB equations.
    d_ib_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
       
    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes); IBTK_CHKERRQ(ierr);
    ierr = SNESSetFunction(snes, composite_res_petsc_vec, SNESFunction_SAMRAI, this); IBTK_CHKERRQ(ierr);
    ierr = SNESSetOptionsPrefix(snes, "ib_"); IBTK_CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); IBTK_CHKERRQ(ierr);
    ierr = SNESSolve(snes, composite_rhs_petsc_vec, composite_sol_petsc_vec); IBTK_CHKERRQ(ierr);
    ierr = SNESDestroy(&snes); IBTK_CHKERRQ(ierr);

    d_ib_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);
    
    // Reset Eulerian solver vectors and Eulerian state data.
    ins_hier_integrator->resetSolverVectors(eul_sol_vec, eul_rhs_vec, current_time, new_time, cycle_num);
    
    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
    if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): interpolating Eulerian velocity to the Lagrangian mesh\n";
    d_ib_method_ops->interpolateVelocity(d_u_idx, getCoarsenSchedules(d_object_name+"::u::CONSERVATIVE_COARSEN"), getGhostfillRefineSchedules(d_object_name+"::u"), half_time);

    // Compute an updated prediction of the updated positions of the Lagrangian
    // structure.
    d_ib_method_ops->midpointStep(current_time, new_time);

    // Deallocate temporary data.
    ierr = VecDestroy(&composite_sol_petsc_vec); IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&composite_rhs_petsc_vec); IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&composite_res_petsc_vec); IBTK_CHKERRQ(ierr);
    PETScSAMRAIVectorReal::destroyPETScVector(eul_sol_petsc_vec);
    PETScSAMRAIVectorReal::destroyPETScVector(eul_rhs_petsc_vec);
    eul_rhs_vec->freeVectorComponents();
    ierr = VecDestroy(&lag_sol_petsc_vec); IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&lag_rhs_petsc_vec); IBTK_CHKERRQ(ierr);
    
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
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data    );
    }

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

    // Finish initializing the hierarchy integrator.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);
}// initializeHierarchyIntegrator

int
IBImplicitStaggeredHierarchyIntegrator::getNumberOfCycles() const
{
    return d_ins_hier_integrator->getNumberOfCycles();
}// getNumberOfCycles

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

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::SNESFunction_SAMRAI(
    SNES snes,
    Vec x,
    Vec f,
    void* ctx)
{
    IBImplicitStaggeredHierarchyIntegrator* ib_integrator = static_cast<IBImplicitStaggeredHierarchyIntegrator*>(ctx);
    return ib_integrator->SNESFunction(snes, x, f);
}// SNESFunction_SAMRAI

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::SNESFunction(
    SNES /*snes*/,
    Vec x,
    Vec f)
{
    PetscErrorCode ierr;
    const int level_num = d_hierarchy->getFinestLevelNumber();

    const double half_time = d_integrator_time+0.5*d_current_dt;
    
    Vec* component_sol_vecs;
    Vec* component_rhs_vecs;
    ierr = VecMultiVecGetSubVecs(x, &component_sol_vecs); CHKERRQ(ierr);
    ierr = VecMultiVecGetSubVecs(f, &component_rhs_vecs); CHKERRQ(ierr);

    Pointer<SAMRAIVectorReal<NDIM,double> > u   = PETScSAMRAIVectorReal::getSAMRAIVector(component_sol_vecs[0]);
    Pointer<SAMRAIVectorReal<NDIM,double> > f_u = PETScSAMRAIVectorReal::getSAMRAIVector(component_rhs_vecs[0]);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = d_ins_hier_integrator->getCurrentContext();
//  Pointer<VariableContext> scratch_ctx = d_ins_hier_integrator->getScratchContext();
//  Pointer<VariableContext> new_ctx     = d_ins_hier_integrator->getNewContext();

    Pointer<Variable<NDIM> > u_var = d_ins_hier_integrator->getVelocityVariable();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, current_ctx);
    const int u_new_idx = u->getComponentDescriptorIndex(0);
    const int f_u_idx = f_u->getComponentDescriptorIndex(0);
    
    Vec X = component_sol_vecs[1];
    Vec F = component_rhs_vecs[1];
    d_ib_method_ops->setSolution(X, level_num);

    // Evaluate the Eulerian terms.
    d_stokes_op->apply(*u, *f_u);

    d_ib_method_ops->computeLagrangianForce(half_time);
    if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
    d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
    d_ib_method_ops->spreadForce(d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name+"::f"), half_time);
    d_hier_velocity_data_ops->subtract(f_u_idx, f_u_idx, d_f_idx);
    
    // Evaluate the Lagrangian terms.
    d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
    d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
    d_ib_method_ops->interpolateVelocity(d_u_idx, getCoarsenSchedules(d_object_name+"::u::CONSERVATIVE_COARSEN"), getGhostfillRefineSchedules(d_object_name+"::u"), half_time);
    d_ib_method_ops->computeResidual(F, level_num);    

    // Ensure that PETSc sees that the state of the RHS vector has changed.
    // This is a nasty hack.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f)); CHKERRQ(ierr);
    return 0;
}// SNESFunction

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
