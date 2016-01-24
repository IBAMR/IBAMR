// Filename: IBImplicitStaggeredHierarchyIntegrator.cpp
// Created on 07 Apr 2012 by Boyce Griffith
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

#include <algorithm>
#include <ostream>
#include <string>

#include "CartesianPatchGeometry.h"
#include "CellData.h"
#include "GriddingAlgorithm.h"
#include "HierarchyDataOpsReal.h"
#include "IntVector.h"
#include "Patch.h"
#include "PatchCellDataOpsReal.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchSideDataOpsReal.h"
#include "PoissonSpecifications.h"
#include "SAMRAIVectorReal.h"
#include "SideData.h"
#include "Variable.h"
#include "VariableContext.h"
#include "VariableDatabase.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBImplicitStaggeredHierarchyIntegrator.h"
#include "ibamr/IBImplicitStrategy.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesOperator.h"
#include "ibamr/StaggeredStokesPETScVecUtilities.h"
#include "ibamr/StaggeredStokesSolver.h"
#include "ibamr/ibamr_enums.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "ibtk/HierarchyMathOps.h"
#include "ibtk/IBTK_CHKERRQ.h"
#include "ibtk/KrylovLinearSolver.h"
#include "ibtk/PETScMultiVec.h"
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_enums.h"
#include "petscerror.h"
#include "petscksp.h"
#include "petscmat.h"
#include "petscpc.h"
#include "petscsnes.h"
#include "petscsys.h"
#include "petscvec.h"
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

namespace IBAMR
{
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{
template <int DIM>
class Box;
} // namespace hier
namespace xfer
{
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Version of IBImplicitStaggeredHierarchyIntegrator restart file data.
static const int IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION = 1;

// IB-4 interpolation function.
static void
ib_4_interp_fcn(const double r, double* const w)
{
    const double q = sqrt(-7.0 + 12.0 * r - 4.0 * r * r);
    w[0] = 0.125 * (5.0 - 2.0 * r - q);
    w[1] = 0.125 * (5.0 - 2.0 * r + q);
    w[2] = 0.125 * (-1.0 + 2.0 * r + q);
    w[3] = 0.125 * (-1.0 + 2.0 * r - q);
    return;
} // ib_4_interp_fcn
static const int ib_4_interp_stencil = 4;

// Piecewise linear interpolation function
static void
pwl_interp_fcn(const double r, double* const w)
{
    w[0] = 1.0 - r;
    w[1] = r;
    return;
} // pwl_interp_fcn
static const int pwl_interp_stencil = 2;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitStaggeredHierarchyIntegrator::IBImplicitStaggeredHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<IBImplicitStrategy> ib_implicit_ops,
    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
    bool register_for_restart)
    : IBHierarchyIntegrator(object_name, input_db, ib_implicit_ops, ins_hier_integrator, register_for_restart),
      d_ib_implicit_ops(ib_implicit_ops)
{
    // Setup IB ops object to use "fixed" Lagrangian-Eulerian coupling
    // operators.
    d_ib_implicit_ops->setUseFixedLEOperators(true);

    // Set default configuration options.
    d_solve_for_position = false;
    d_jac_delta_fcn = "IB_4";

    // Set options from input.
    if (input_db)
    {
        if (input_db->keyExists("eliminate_eulerian_vars"))
            d_solve_for_position = input_db->getBool("eliminate_eulerian_vars");
        else if (input_db->keyExists("eliminate_lagrangian_vars"))
            d_solve_for_position = !input_db->getBool("eliminate_lagrangian_vars");

        if (input_db->keyExists("jacobian_delta_fcn")) d_jac_delta_fcn = input_db->getString("jacobian_delta_fcn");
    }

    if (!d_solve_for_position)
    {
        Pointer<Database> stokes_ib_pc_db = input_db->getDatabase("stokes_ib_precond_db");
        d_stokes_solver = new IBImplicitStaggeredStokesSolver(object_name, stokes_ib_pc_db);
        ins_hier_integrator->setStokesSolver(d_stokes_solver);
    }

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;
} // IBImplicitStaggeredHierarchyIntegrator

IBImplicitStaggeredHierarchyIntegrator::~IBImplicitStaggeredHierarchyIntegrator()
{
    // intentionally blank
    return;
} // ~IBImplicitStaggeredHierarchyIntegrator

void
IBImplicitStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                                     const double new_time,
                                                                     const int num_cycles)
{
    IBHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    if (d_time_stepping_type != MIDPOINT_RULE)
    {
        TBOX_ERROR("IBImplicitStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(): time_stepping_type = "
                   << enum_to_string<TimeSteppingType>(d_time_stepping_type)
                   << "\n"
                   << "  currently only time_stepping_type = "
                   << enum_to_string<TimeSteppingType>(MIDPOINT_RULE)
                   << " is supported\n");
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate Eulerian scratch and new data.
    d_num_dofs_per_proc.resize(finest_ln + 1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_u_idx, current_time);
        level->allocatePatchData(d_f_idx, current_time);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
        if (!d_solve_for_position && ln == finest_ln)
        {
            level->allocatePatchData(d_u_dof_index_idx, current_time);
            level->allocatePatchData(d_p_dof_index_idx, current_time);
            StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
                d_num_dofs_per_proc[ln], d_u_dof_index_idx, d_p_dof_index_idx, level);
        }
    }

    // Initialize IB data.
    d_ib_implicit_ops->preprocessIntegrateData(current_time, new_time, num_cycles);

    // Initialize the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    if (ins_num_cycles != d_current_num_cycles && d_current_num_cycles != 1)
    {
        TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                 << "  attempting to perform "
                                 << d_current_num_cycles
                                 << " cycles of fixed point iteration.\n"
                                 << "  number of cycles required by Navier-Stokes solver = "
                                 << ins_num_cycles
                                 << ".\n"
                                 << "  current implementation requires either that both solvers "
                                    "use the same number of cycles,\n"
                                 << "  or that the IB solver use only a single cycle.\n");
    }
    d_ins_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, ins_num_cycles);

    // Compute an initial prediction of the updated positions of the Lagrangian
    // structure.
    //
    // NOTE: The velocity should already have been interpolated to the
    // curvilinear mesh and should not need to be re-interpolated.
    if (d_enable_logging)
        plog << d_object_name << "::preprocessIntegrateHierarchy(): performing Lagrangian forward Euler step\n";
    d_ib_implicit_ops->eulerStep(current_time, new_time);

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
IBImplicitStaggeredHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                           const double new_time,
                                                           const int cycle_num)
{
    if (d_solve_for_position)
    {
        integrateHierarchy_position(current_time, new_time, cycle_num);
    }
    else
    {
        integrateHierarchy_velocity(current_time, new_time, cycle_num);
    }

    return;
} // integrateHierarchy

void
IBImplicitStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                      const double new_time,
                                                                      const bool skip_synchronize_new_state_data,
                                                                      const int num_cycles)
{
    IBHierarchyIntegrator::postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): interpolating Eulerian "
                                 "velocity to the Lagrangian mesh\n";
    d_ib_implicit_ops->interpolateVelocity(d_u_idx,
                                           getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                           getGhostfillRefineSchedules(d_object_name + "::u"),
                                           new_time);

    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_enable_logging)
            plog << d_object_name << "::postprocessIntegrateHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(NEW_DATA);
    }

    // Determine the CFL number.
    double cfl_max = 0.0;
    PatchCellDataOpsReal<NDIM, double> patch_cc_ops;
    PatchSideDataOpsReal<NDIM, double> patch_sc_ops;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double dx_min = *(std::min_element(dx, dx + NDIM));
            Pointer<CellData<NDIM, double> > u_cc_new_data = patch->getPatchData(u_new_idx);
            Pointer<SideData<NDIM, double> > u_sc_new_data = patch->getPatchData(u_new_idx);
            double u_max = 0.0;
            if (u_cc_new_data) u_max = patch_cc_ops.maxNorm(u_cc_new_data, patch_box);
            if (u_sc_new_data) u_max = patch_sc_ops.maxNorm(u_sc_new_data, patch_box);
            cfl_max = std::max(cfl_max, u_max * dt / dx_min);
        }
    }
    cfl_max = SAMRAI_MPI::maxReduction(cfl_max);
    d_regrid_cfl_estimate += cfl_max;
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl_max << "\n";
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): estimated upper bound on IB "
                                 "point displacement since last regrid = "
             << d_regrid_cfl_estimate << "\n";

    // Deallocate the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    d_ins_hier_integrator->postprocessIntegrateHierarchy(
        current_time, new_time, skip_synchronize_new_state_data, ins_num_cycles);

    // Deallocate IB data.
    d_ib_implicit_ops->postprocessIntegrateData(current_time, new_time, num_cycles);

    // Deallocate Eulerian scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_u_idx);
        level->deallocatePatchData(d_f_idx);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
        if (!d_solve_for_position && ln == finest_ln)
        {
            level->deallocatePatchData(d_u_dof_index_idx);
            level->deallocatePatchData(d_p_dof_index_idx);
        }
    }

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
IBImplicitStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                                      Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Register u and p DOF variables.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_u_dof_index_var = new SideVariable<NDIM, int>(d_object_name + "::u_dof_index");
    d_p_dof_index_var = new CellVariable<NDIM, int>(d_object_name + "::p_dof_index");
    const IntVector<NDIM> ib_ghosts = d_ib_method_ops->getMinimumGhostCellWidth();
    const IntVector<NDIM> no_ghosts = 0;
    d_u_dof_index_idx = var_db->registerVariableAndContext(d_u_dof_index_var, getScratchContext(), ib_ghosts);
    d_p_dof_index_idx = var_db->registerVariableAndContext(d_p_dof_index_var, getScratchContext(), no_ghosts);

    // Register body force function with INSHierarchyIntegrator
    d_ins_hier_integrator->registerBodyForceFunction(d_body_force_fcn);

    // Finish initializing the hierarchy integrator.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);
} // initializeHierarchyIntegrator

int
IBImplicitStaggeredHierarchyIntegrator::getNumberOfCycles() const
{
    return d_ins_hier_integrator->getNumberOfCycles();
} // getNumberOfCycles

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBImplicitStaggeredHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    IBHierarchyIntegrator::putToDatabaseSpecialized(db);
    db->putInteger("IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION",
                   IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION);
    return;
} // putToDatabaseSpecialized

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
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file."
                                 << std::endl);
    }
    int ver = db->getInteger("IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
} // getFromRestart

void
IBImplicitStaggeredHierarchyIntegrator::integrateHierarchy_position(const double current_time,
                                                                    const double new_time,
                                                                    const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);

    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator = d_ins_hier_integrator;
    TBOX_ASSERT(ins_hier_integrator);

    PetscErrorCode ierr;
    int n_local;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = ins_hier_integrator->getCurrentContext();
    Pointer<VariableContext> scratch_ctx = ins_hier_integrator->getScratchContext();
    Pointer<VariableContext> new_ctx = ins_hier_integrator->getNewContext();

    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();

    Pointer<Variable<NDIM> > u_var = ins_hier_integrator->getVelocityVariable();
    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, scratch_ctx);

    Pointer<Variable<NDIM> > p_var = ins_hier_integrator->getPressureVariable();
    const int p_scratch_idx = var_db->mapVariableAndContextToIndex(p_var, scratch_ctx);

    // Skip all cycles in the INS solver --- we advance the state data here.
    ins_hier_integrator->skipCycle(current_time, new_time, cycle_num);

    // Setup Eulerian vectors used in solving the implicit IB equations.
    Pointer<SAMRAIVectorReal<NDIM, double> > eul_sol_vec =
        new SAMRAIVectorReal<NDIM, double>(d_object_name + "::eulerian_sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    eul_sol_vec->addComponent(u_var, u_scratch_idx, wgt_sc_idx, d_hier_velocity_data_ops);
    eul_sol_vec->addComponent(p_var, p_scratch_idx, wgt_cc_idx, d_hier_pressure_data_ops);

    Pointer<SAMRAIVectorReal<NDIM, double> > eul_rhs_vec =
        eul_sol_vec->cloneVector(d_object_name + "::eulerian_rhs_vec");
    eul_rhs_vec->allocateVectorData(current_time);

    d_u_scratch_vec = eul_sol_vec->cloneVector(d_object_name + "::u_scratch_vec");
    d_f_scratch_vec = eul_rhs_vec->cloneVector(d_object_name + "::f_scratch_vec");
    d_u_scratch_vec->allocateVectorData(current_time);
    d_f_scratch_vec->allocateVectorData(current_time);

    ins_hier_integrator->setupSolverVectors(eul_sol_vec, eul_rhs_vec, current_time, new_time, cycle_num);

    d_stokes_solver = ins_hier_integrator->getStokesSolver();
    Pointer<KrylovLinearSolver> p_stokes_solver = d_stokes_solver;
    TBOX_ASSERT(p_stokes_solver);
    d_stokes_op = p_stokes_solver->getOperator();
    TBOX_ASSERT(d_stokes_op);

    // Setup Lagrangian vectors used in solving the implicit IB equations.
    Vec lag_sol_petsc_vec, lag_rhs_petsc_vec;
    d_ib_implicit_ops->createSolverVecs(&lag_sol_petsc_vec, &lag_rhs_petsc_vec);
    d_ib_implicit_ops->setupSolverVecs(&lag_sol_petsc_vec, &lag_rhs_petsc_vec);

    // Indicate that the current approximation to position of the structure
    // should be used for Lagrangian-Eulerian coupling.
    d_ib_implicit_ops->updateFixedLEOperators();

    // Setup multi-vec objects to store the composite solution and
    // right-hand-side vectors.
    Vec eul_sol_petsc_vec = PETScSAMRAIVectorReal::createPETScVector(eul_sol_vec, PETSC_COMM_WORLD);
    Vec eul_rhs_petsc_vec = PETScSAMRAIVectorReal::createPETScVector(eul_rhs_vec, PETSC_COMM_WORLD);

    Vec sol_petsc_vecs[] = { eul_sol_petsc_vec, lag_sol_petsc_vec };
    Vec rhs_petsc_vecs[] = { eul_rhs_petsc_vec, lag_rhs_petsc_vec };

    Vec composite_sol_petsc_vec, composite_rhs_petsc_vec, composite_res_petsc_vec;
    ierr = VecCreateMultiVec(PETSC_COMM_WORLD, 2, sol_petsc_vecs, &composite_sol_petsc_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMultiVec(PETSC_COMM_WORLD, 2, rhs_petsc_vecs, &composite_rhs_petsc_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(composite_rhs_petsc_vec, &composite_res_petsc_vec);
    IBTK_CHKERRQ(ierr);

    // Solve the implicit IB equations.
    d_ib_implicit_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);

    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
    IBTK_CHKERRQ(ierr);
    ierr = SNESSetFunction(snes, composite_res_petsc_vec, IBFunction_SAMRAI, this);
    IBTK_CHKERRQ(ierr);
    ierr = SNESSetOptionsPrefix(snes, "ib_");
    IBTK_CHKERRQ(ierr);

    Mat jac;
    ierr = VecGetLocalSize(composite_sol_petsc_vec, &n_local);
    IBTK_CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD, n_local, n_local, PETSC_DETERMINE, PETSC_DETERMINE, this, &jac);
    IBTK_CHKERRQ(ierr);
    ierr = MatShellSetOperation(jac, MATOP_MULT, reinterpret_cast<void (*)(void)>(IBJacobianApply_SAMRAI));
    IBTK_CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes, jac, jac, IBJacobianSetup_SAMRAI, this);
    IBTK_CHKERRQ(ierr);

    Mat schur;
    ierr = VecGetLocalSize(lag_sol_petsc_vec, &n_local);
    IBTK_CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD, n_local, n_local, PETSC_DETERMINE, PETSC_DETERMINE, this, &schur);
    IBTK_CHKERRQ(ierr);
    ierr = MatShellSetOperation(schur, MATOP_MULT, reinterpret_cast<void (*)(void)>(lagrangianSchurApply_SAMRAI));
    IBTK_CHKERRQ(ierr);
    ierr = KSPCreate(PETSC_COMM_WORLD, &d_schur_solver);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetOptionsPrefix(d_schur_solver, "ib_schur_");
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetOperators(d_schur_solver, schur, schur);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetReusePreconditioner(d_schur_solver, PETSC_TRUE);
    IBTK_CHKERRQ(ierr);
    PC schur_pc;
    ierr = KSPGetPC(d_schur_solver, &schur_pc);
    IBTK_CHKERRQ(ierr);
    ierr = PCSetType(schur_pc, PCNONE);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetFromOptions(d_schur_solver);
    IBTK_CHKERRQ(ierr);

    KSP snes_ksp;
    ierr = SNESGetKSP(snes, &snes_ksp);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetType(snes_ksp, KSPFGMRES);
    IBTK_CHKERRQ(ierr);
    PC snes_pc;
    ierr = KSPGetPC(snes_ksp, &snes_pc);
    IBTK_CHKERRQ(ierr);
    ierr = PCSetType(snes_pc, PCSHELL);
    IBTK_CHKERRQ(ierr);
    ierr = PCShellSetContext(snes_pc, this);
    IBTK_CHKERRQ(ierr);
    ierr = PCShellSetApply(snes_pc, IBPCApply_SAMRAI);
    IBTK_CHKERRQ(ierr);

    ierr = SNESSetFromOptions(snes);
    IBTK_CHKERRQ(ierr);
    ierr = SNESSolve(snes, composite_rhs_petsc_vec, composite_sol_petsc_vec);
    IBTK_CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&jac);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&schur);
    IBTK_CHKERRQ(ierr);
    ierr = KSPDestroy(&d_schur_solver);
    IBTK_CHKERRQ(ierr);

    d_ib_implicit_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);

    // Reset Eulerian solver vectors and Eulerian state data.
    ins_hier_integrator->resetSolverVectors(eul_sol_vec, eul_rhs_vec, current_time, new_time, cycle_num);

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_ib_implicit_ops->setUpdatedPosition(lag_sol_petsc_vec);

    // Deallocate temporary data.
    ierr = VecDestroy(&composite_sol_petsc_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&composite_rhs_petsc_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&composite_res_petsc_vec);
    IBTK_CHKERRQ(ierr);
    PETScSAMRAIVectorReal::destroyPETScVector(eul_sol_petsc_vec);
    PETScSAMRAIVectorReal::destroyPETScVector(eul_rhs_petsc_vec);
    eul_rhs_vec->freeVectorComponents();
    d_u_scratch_vec->freeVectorComponents();
    d_f_scratch_vec->freeVectorComponents();
    ierr = VecDestroy(&lag_sol_petsc_vec);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&lag_rhs_petsc_vec);
    IBTK_CHKERRQ(ierr);

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy_position

void
IBImplicitStaggeredHierarchyIntegrator::integrateHierarchy_velocity(const double current_time,
                                                                    const double new_time,
                                                                    const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);

    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator = d_ins_hier_integrator;
    TBOX_ASSERT(ins_hier_integrator);

    PetscErrorCode ierr;
    int n_local;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double half_time = current_time + 0.5 * (new_time - current_time);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = ins_hier_integrator->getCurrentContext();
    Pointer<VariableContext> scratch_ctx = ins_hier_integrator->getScratchContext();
    Pointer<VariableContext> new_ctx = ins_hier_integrator->getNewContext();

    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();

    Pointer<Variable<NDIM> > u_var = ins_hier_integrator->getVelocityVariable();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, current_ctx);
    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, scratch_ctx);

    Pointer<Variable<NDIM> > p_var = ins_hier_integrator->getPressureVariable();
    const int p_scratch_idx = var_db->mapVariableAndContextToIndex(p_var, scratch_ctx);

    // Skip all cycles in the INS solver --- we advance the state data here.
    ins_hier_integrator->skipCycle(current_time, new_time, cycle_num);

    // Setup Eulerian vectors used in solving the implicit IB equations.
    Pointer<SAMRAIVectorReal<NDIM, double> > eul_sol_vec =
        new SAMRAIVectorReal<NDIM, double>(d_object_name + "::eulerian_sol_vec", d_hierarchy, coarsest_ln, finest_ln);
    eul_sol_vec->addComponent(u_var, u_scratch_idx, wgt_sc_idx, d_hier_velocity_data_ops);
    eul_sol_vec->addComponent(p_var, p_scratch_idx, wgt_cc_idx, d_hier_pressure_data_ops);

    Pointer<SAMRAIVectorReal<NDIM, double> > eul_rhs_vec =
        eul_sol_vec->cloneVector(d_object_name + "::eulerian_rhs_vec");
    eul_rhs_vec->allocateVectorData(current_time);

    d_f_scratch_vec = eul_rhs_vec->cloneVector(d_object_name + "::f_scratch_vec");
    d_f_scratch_vec->allocateVectorData(current_time);

    // Compute convective and previous time-step diffusive terms in the rhs vec.
    ins_hier_integrator->setupSolverVectors(eul_sol_vec, eul_rhs_vec, current_time, new_time, cycle_num);

    // Setup Stokes operator.
    Pointer<IBImplicitStaggeredStokesSolver> p_stokes_solver = d_stokes_solver;
#if !defined(NDEBUG)
    TBOX_ASSERT(p_stokes_solver);
#endif
    d_stokes_op = p_stokes_solver->getStaggeredStokesOperator();
    d_stokes_op->setSolutionTime(new_time);
    d_stokes_op->setTimeInterval(current_time, new_time);
    d_stokes_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);

    // Setup Stokes-IB preconditioner.
    Pointer<StaggeredStokesFACPreconditioner> stokes_fac_pc = p_stokes_solver->getStaggeredStokesFACPreconditioner();
    stokes_fac_pc->setSolutionTime(new_time);
    stokes_fac_pc->setTimeInterval(current_time, new_time);
    stokes_fac_pc->setPhysicalBcCoefs(d_ins_hier_integrator->getIntermediateVelocityBoundaryConditions(),
                                      d_ins_hier_integrator->getProjectionBoundaryConditions());
    Pointer<StaggeredStokesIBLevelRelaxationFACOperator> stokes_fac_op = stokes_fac_pc->getFACPreconditionerStrategy();
    Mat elastic_op = NULL;
    d_ib_implicit_ops->constructLagrangianForceJacobian(elastic_op, MATAIJ);
    stokes_fac_op->setIBForceJacobian(elastic_op);
    Mat interp_op = NULL;
    if (d_jac_delta_fcn == "IB_4")
    {
        d_ib_implicit_ops->constructInterpOp(
            interp_op, ib_4_interp_fcn, ib_4_interp_stencil, d_num_dofs_per_proc[finest_ln], d_u_dof_index_idx);
    }
    else if (d_jac_delta_fcn == "PIECEWISE_LINEAR")
    {
        d_ib_implicit_ops->constructInterpOp(
            interp_op, pwl_interp_fcn, pwl_interp_stencil, d_num_dofs_per_proc[finest_ln], d_u_dof_index_idx);
    }
    else
    {
        TBOX_ERROR("IBImplicitStaggeredHierarchyIntegrator::integrateHierarchy_velocity()."
                   << " Delta function "
                   << d_jac_delta_fcn
                   << " is not supported in creating Jacobian."
                   << std::endl);
    }
    stokes_fac_op->setIBInterpOp(interp_op);
    stokes_fac_pc->initializeSolverState(*eul_sol_vec, *eul_rhs_vec);

    // Indicate that the current approximation to position of the structure
    // should be used for Lagrangian-Eulerian coupling.
    d_ib_implicit_ops->updateFixedLEOperators();

    // Setup right-hand-side vectors.
    d_ib_implicit_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    d_stokes_op->setHomogeneousBc(false);
    d_stokes_op->modifyRhsForBcs(*eul_rhs_vec);
    d_stokes_solver->setHomogeneousBc(true);
    d_stokes_op->setHomogeneousBc(true);

    // Create PETSc Vecs to be used with PETSc solvers.
    d_ib_implicit_ops->createSolverVecs(&d_X_current, NULL);
    d_ib_implicit_ops->setupSolverVecs(&d_X_current, NULL);
    Vec eul_sol_petsc_vec = PETScSAMRAIVectorReal::createPETScVector(eul_sol_vec, PETSC_COMM_WORLD);
    Vec eul_rhs_petsc_vec = PETScSAMRAIVectorReal::createPETScVector(eul_rhs_vec, PETSC_COMM_WORLD);
    Vec eul_res_petsc_vec = PETScSAMRAIVectorReal::createPETScVector(d_f_scratch_vec, PETSC_COMM_WORLD);

    // Create the outer nonlinear solver.
    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
    IBTK_CHKERRQ(ierr);
    ierr = SNESSetFunction(snes, eul_res_petsc_vec, IBFunction_SAMRAI, this);
    IBTK_CHKERRQ(ierr);
    ierr = SNESSetOptionsPrefix(snes, "ib_");
    IBTK_CHKERRQ(ierr);

    // Create the Jacobian for Newton iterations.
    Mat jac;
    ierr = VecGetLocalSize(eul_sol_petsc_vec, &n_local);
    IBTK_CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD, n_local, n_local, PETSC_DETERMINE, PETSC_DETERMINE, this, &jac);
    IBTK_CHKERRQ(ierr);
    ierr = MatShellSetOperation(jac, MATOP_MULT, reinterpret_cast<void (*)(void)>(IBJacobianApply_SAMRAI));
    IBTK_CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes, jac, jac, IBJacobianSetup_SAMRAI, this);
    IBTK_CHKERRQ(ierr);

    // Create the KSP for Jacobian setup for SNES.
    KSP snes_ksp;
    ierr = SNESGetKSP(snes, &snes_ksp);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSetType(snes_ksp, KSPFGMRES);
    IBTK_CHKERRQ(ierr);
    PC snes_pc;
    ierr = KSPGetPC(snes_ksp, &snes_pc);
    IBTK_CHKERRQ(ierr);
    ierr = PCSetType(snes_pc, PCSHELL);
    IBTK_CHKERRQ(ierr);
    ierr = PCShellSetContext(snes_pc, this);
    IBTK_CHKERRQ(ierr);
    ierr = PCShellSetApply(snes_pc, IBPCApply_SAMRAI);
    IBTK_CHKERRQ(ierr);

    // Solve the system.
    ierr = SNESSetFromOptions(snes);
    IBTK_CHKERRQ(ierr);
    ierr = SNESSolve(snes, eul_rhs_petsc_vec, eul_sol_petsc_vec);
    IBTK_CHKERRQ(ierr);
    ierr = SNESDestroy(&snes);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&jac);
    IBTK_CHKERRQ(ierr);

    d_stokes_op->imposeSolBcs(*eul_sol_vec);
    d_ib_implicit_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);

    // Reset Eulerian solver vectors and Eulerian state data.
    ins_hier_integrator->resetSolverVectors(eul_sol_vec, eul_rhs_vec, current_time, new_time, cycle_num);

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_scratch_idx);
    if (d_enable_logging)
    {
        plog << d_object_name
             << "::integrateHierarchy_velocity(): interpolating Eulerian velocity to the Lagrangian mesh\n";
    }
    d_ib_implicit_ops->interpolateVelocity(d_u_idx,
                                           getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                           getGhostfillRefineSchedules(d_object_name + "::u"),
                                           half_time);

    // Compute the final value of the updated positions of the Lagrangian
    // structure.
    d_ib_implicit_ops->midpointStep(current_time, new_time);

    // Deallocate PETSc Vecs.
    ierr = VecDestroy(&d_X_current);
    IBTK_CHKERRQ(ierr);
    PETScSAMRAIVectorReal::destroyPETScVector(eul_sol_petsc_vec);
    PETScSAMRAIVectorReal::destroyPETScVector(eul_rhs_petsc_vec);
    PETScSAMRAIVectorReal::destroyPETScVector(eul_res_petsc_vec);

    // Deallocate Eulerian components.
    eul_rhs_vec->freeVectorComponents();
    d_f_scratch_vec->freeVectorComponents();

    // Deallocate solvers and operators.
    p_stokes_solver->deallocateSolverState();
    d_stokes_op->deallocateOperatorState();
    stokes_fac_pc->deallocateSolverState();
    stokes_fac_op->deallocateOperatorState();
    ierr = MatDestroy(&elastic_op);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&interp_op);
    IBTK_CHKERRQ(ierr);

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
    return;
} // integrateHierarchy_velocity

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBFunction_SAMRAI(SNES snes, Vec x, Vec f, void* ctx)
{
    IBImplicitStaggeredHierarchyIntegrator* ib_integrator = static_cast<IBImplicitStaggeredHierarchyIntegrator*>(ctx);

    PetscErrorCode ierr = 1;
    if (ib_integrator->d_solve_for_position)
    {
        ierr = ib_integrator->IBFunction_position(snes, x, f);
    }
    else
    {
        ierr = ib_integrator->IBFunction_velocity(snes, x, f);
    }

    return ierr;
} // IBFunction_SAMRAI

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBFunction_position(SNES /*snes*/, Vec x, Vec f)
{
    PetscErrorCode ierr;
    const double half_time = d_integrator_time + 0.5 * d_current_dt;

    Vec* component_sol_vecs;
    Vec* component_rhs_vecs;
    ierr = VecMultiVecGetSubVecs(x, &component_sol_vecs);
    IBTK_CHKERRQ(ierr);
    ierr = VecMultiVecGetSubVecs(f, &component_rhs_vecs);
    IBTK_CHKERRQ(ierr);

    Pointer<SAMRAIVectorReal<NDIM, double> > u = PETScSAMRAIVectorReal::getSAMRAIVector(component_sol_vecs[0]);
    Pointer<SAMRAIVectorReal<NDIM, double> > f_u = PETScSAMRAIVectorReal::getSAMRAIVector(component_rhs_vecs[0]);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = d_ins_hier_integrator->getCurrentContext();
    Pointer<Variable<NDIM> > u_var = d_ins_hier_integrator->getVelocityVariable();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, current_ctx);
    const int u_new_idx = u->getComponentDescriptorIndex(0);
    const int f_u_idx = f_u->getComponentDescriptorIndex(0);

    Vec X = component_sol_vecs[1];
    Vec R = component_rhs_vecs[1];

    // Evaluate the Eulerian terms.
    d_stokes_op->setHomogeneousBc(true);
    d_stokes_op->apply(*u, *f_u);

    d_ib_implicit_ops->setUpdatedPosition(X);
    d_ib_implicit_ops->computeLagrangianForce(half_time);
    if (d_enable_logging)
    {
        plog << d_object_name << "::integrateHierarchy_position(): spreading Lagrangian force to the Eulerian grid\n";
        plog << "Spreading being done from " << d_object_name << "::IBFunction_position().\n";
    }
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0, /*interior_only*/ false);
    d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
    d_u_phys_bdry_op->setHomogeneousBc(true); // use homogeneous BCs to define spreading at physical boundaries
    d_ib_implicit_ops->spreadForce(
        d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);
    d_hier_velocity_data_ops->subtract(f_u_idx, f_u_idx, d_f_idx);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(component_rhs_vecs[0]));
    IBTK_CHKERRQ(ierr);

    // Evaluate the Lagrangian terms.
    d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
    d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
    d_u_phys_bdry_op->setHomogeneousBc(false);
    d_ib_implicit_ops->interpolateVelocity(d_u_idx,
                                           getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                           getGhostfillRefineSchedules(d_object_name + "::u"),
                                           half_time);
    d_ib_implicit_ops->computeResidual(R);

    // Ensure that PETSc sees that the state of the RHS vector has changed.
    // This is a nasty hack.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f));
    IBTK_CHKERRQ(ierr);

    return ierr;
} // IBFunction_position

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBFunction_velocity(SNES /*snes*/, Vec x, Vec f)
{
    PetscErrorCode ierr;
    const double current_time = d_integrator_time;
    const double new_time = current_time + d_current_dt;
    const double half_time = d_integrator_time + 0.5 * d_current_dt;

    Pointer<SAMRAIVectorReal<NDIM, double> > u = PETScSAMRAIVectorReal::getSAMRAIVector(x);
    Pointer<SAMRAIVectorReal<NDIM, double> > f_u = PETScSAMRAIVectorReal::getSAMRAIVector(f);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = d_ins_hier_integrator->getCurrentContext();
    Pointer<Variable<NDIM> > u_var = d_ins_hier_integrator->getVelocityVariable();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, current_ctx);
    const int u_new_idx = u->getComponentDescriptorIndex(0);
    const int f_u_idx = f_u->getComponentDescriptorIndex(0);

    // Apply the Stokes part.
    d_stokes_op->setHomogeneousBc(true);
    d_stokes_op->apply(*u, *f_u);

    // Compute the new position of the structure.
    d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_new_idx, 0.5, u_current_idx);
    d_ib_implicit_ops->interpolateVelocity(d_u_idx,
                                           getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                           getGhostfillRefineSchedules(d_object_name + "::u"),
                                           half_time);
    d_ib_implicit_ops->midpointStep(current_time, new_time);

    // Subtract the spreaded F[X^n+1/2] from the RHS
    d_ib_implicit_ops->computeLagrangianForce(half_time);
    if (d_enable_logging)
    {
        plog << d_object_name << "::integrateHierarchy_velocity(): spreading Lagrangian force to the Eulerian grid\n";
        plog << "Spreading being done from " << d_object_name << "::IBFunction_velocity().\n";
    }
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0, /*interior_only*/ false);
    d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
    d_u_phys_bdry_op->setHomogeneousBc(true); // use homogeneous BCs to define spreading at physical boundaries
    d_ib_implicit_ops->spreadForce(
        d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);
    d_hier_velocity_data_ops->subtract(f_u_idx, f_u_idx, d_f_idx);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f));
    IBTK_CHKERRQ(ierr);

    return ierr;
} // IBFunction_velocity

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBJacobianSetup_SAMRAI(SNES snes, Vec x, Mat A, Mat B, void* ctx)
{
    IBImplicitStaggeredHierarchyIntegrator* ib_integrator = static_cast<IBImplicitStaggeredHierarchyIntegrator*>(ctx);

    PetscErrorCode ierr = 1;
    if (ib_integrator->d_solve_for_position)
    {
        ierr = ib_integrator->IBJacobianSetup_position(snes, x, A, B);
    }
    else
    {
        ierr = ib_integrator->IBJacobianSetup_velocity(snes, x, A, B);
    }

    return ierr;
} // IBJacobianSetup_SAMRAI

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBJacobianSetup_position(SNES /*snes*/, Vec x, Mat A, Mat /*B*/)
{
    PetscErrorCode ierr;
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    Vec* component_sol_vecs;
    ierr = VecMultiVecGetSubVecs(x, &component_sol_vecs);
    IBTK_CHKERRQ(ierr);
    Vec X = component_sol_vecs[1];
    d_ib_implicit_ops->setLinearizedPosition(X);

    return ierr;
} // IBJacobianSetup_position

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBJacobianSetup_velocity(SNES /*snes*/, Vec x, Mat A, Mat /*B*/)
{
    PetscErrorCode ierr;
    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);

    // Get the estimate of X^{n+1} from the current iterate U^{n+1} and set as it as
    // a base vector in matrix-free Lagrangian force Jacobian.
    const double half_time = d_integrator_time + 0.5 * d_current_dt;
    Vec X_new;
    ierr = VecDuplicate(d_X_current, &X_new);
    IBTK_CHKERRQ(ierr);

    Pointer<SAMRAIVectorReal<NDIM, double> > u = PETScSAMRAIVectorReal::getSAMRAIVector(x);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> current_ctx = d_ins_hier_integrator->getCurrentContext();
    Pointer<Variable<NDIM> > u_var = d_ins_hier_integrator->getVelocityVariable();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, current_ctx);
    const int u_new_idx = u->getComponentDescriptorIndex(0);
    d_hier_velocity_data_ops->linearSum(d_u_idx, -0.5, u_new_idx, -0.5, u_current_idx);

    d_ib_implicit_ops->interpolateLinearizedVelocity(d_u_idx,
                                                     getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                                     getGhostfillRefineSchedules(d_object_name + "::u"),
                                                     half_time);
    d_ib_implicit_ops->computeLinearizedResidual(d_X_current, X_new);
    d_ib_implicit_ops->setLinearizedPosition(X_new);

    ierr = VecDestroy(&X_new);
    IBTK_CHKERRQ(ierr);

    return ierr;
} // IBJacobianSetup_velocity

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBJacobianApply_SAMRAI(Mat A, Vec x, Vec f)
{
    PetscErrorCode ierr = 1;
    void* ctx;
    ierr = MatShellGetContext(A, &ctx);
    IBTK_CHKERRQ(ierr);
    IBImplicitStaggeredHierarchyIntegrator* ib_integrator = static_cast<IBImplicitStaggeredHierarchyIntegrator*>(ctx);
    if (ib_integrator->d_solve_for_position)
    {
        ierr = ib_integrator->IBJacobianApply_position(x, f);
    }
    else
    {
        ierr = ib_integrator->IBJacobianApply_velocity(x, f);
    }

    return ierr;
} // IBJacobianApply_SAMRAI

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBJacobianApply_position(Vec x, Vec f)
{
    PetscErrorCode ierr;
    const double half_time = d_integrator_time + 0.5 * d_current_dt;

    Vec* component_sol_vecs;
    Vec* component_rhs_vecs;
    ierr = VecMultiVecGetSubVecs(x, &component_sol_vecs);
    IBTK_CHKERRQ(ierr);
    ierr = VecMultiVecGetSubVecs(f, &component_rhs_vecs);
    IBTK_CHKERRQ(ierr);

    Pointer<SAMRAIVectorReal<NDIM, double> > u = PETScSAMRAIVectorReal::getSAMRAIVector(component_sol_vecs[0]);
    Pointer<SAMRAIVectorReal<NDIM, double> > f_u = PETScSAMRAIVectorReal::getSAMRAIVector(component_rhs_vecs[0]);

    Pointer<Variable<NDIM> > u_var = d_ins_hier_integrator->getVelocityVariable();
    const int u_idx = u->getComponentDescriptorIndex(0);
    const int f_u_idx = f_u->getComponentDescriptorIndex(0);

    Vec X = component_sol_vecs[1];
    Vec R = component_rhs_vecs[1];

    // Evaluate the Eulerian terms.
    d_stokes_op->setHomogeneousBc(true);
    d_stokes_op->apply(*u, *f_u);

    d_ib_implicit_ops->computeLinearizedLagrangianForce(X, half_time);
    if (d_enable_logging)
    {
        plog << d_object_name << "::integrateHierarchy_position(): spreading Lagrangian force to the Eulerian grid\n";
        plog << "Spreading being done from " << d_object_name << "::IBJacobianApply_position().\n";
    }
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0, /*interior_only*/ false);
    d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
    d_u_phys_bdry_op->setHomogeneousBc(true); // use homogeneous BCs to define spreading at physical boundaries
    d_ib_implicit_ops->spreadLinearizedForce(
        d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);
    d_hier_velocity_data_ops->subtract(f_u_idx, f_u_idx, d_f_idx);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(component_rhs_vecs[0]));
    IBTK_CHKERRQ(ierr);

    // Evaluate the Lagrangian terms.
    d_hier_velocity_data_ops->scale(d_u_idx, 0.5, u_idx);
    d_ib_implicit_ops->interpolateLinearizedVelocity(d_u_idx,
                                                     getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                                     getGhostfillRefineSchedules(d_object_name + "::u"),
                                                     half_time);
    d_ib_implicit_ops->computeLinearizedResidual(X, R);

    // Ensure that PETSc sees that the state of the RHS vector has changed.
    // This is a nasty hack.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f));
    IBTK_CHKERRQ(ierr);

    return ierr;
} // IBJacobianApply_position

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBJacobianApply_velocity(Vec x, Vec f)
{
    PetscErrorCode ierr;
    const double half_time = d_integrator_time + 0.5 * d_current_dt;

    Pointer<SAMRAIVectorReal<NDIM, double> > u = PETScSAMRAIVectorReal::getSAMRAIVector(x);
    Pointer<SAMRAIVectorReal<NDIM, double> > f_u = PETScSAMRAIVectorReal::getSAMRAIVector(f);

    const int u_idx = u->getComponentDescriptorIndex(0);
    const int f_u_idx = f_u->getComponentDescriptorIndex(0);

    // Evaluate the Eulerian terms.
    d_stokes_op->setHomogeneousBc(true);
    d_stokes_op->apply(*u, *f_u);

    // Compute position residual X = dt*J*[u/2] = 0 - dt*J*[-u/2]
    Vec X, X0;
    d_ib_implicit_ops->createSolverVecs(&X, &X0);
    d_ib_implicit_ops->setupSolverVecs(NULL, &X0);

    d_hier_velocity_data_ops->scale(d_u_idx, -0.5, u_idx);
    d_ib_implicit_ops->interpolateLinearizedVelocity(d_u_idx,
                                                     getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                                     getGhostfillRefineSchedules(d_object_name + "::u"),
                                                     half_time);
    d_ib_implicit_ops->computeLinearizedResidual(X0, X);

    // Compute linearized force F = A/2*[dt*J*[u/2]]
    d_ib_implicit_ops->computeLinearizedLagrangianForce(X, half_time);
    if (d_enable_logging)
    {
        plog << d_object_name << "::integrateHierarchy_velocity(): spreading Lagrangian force to the Eulerian grid\n";
        plog << "Spreading being done from " << d_object_name << "::IBJacobianApply_velocity().\n";
    }
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0, /*interior_only*/ false);
    d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
    d_u_phys_bdry_op->setHomogeneousBc(true); // use homogeneous BCs to define spreading at physical boundaries
    d_ib_implicit_ops->spreadLinearizedForce(
        d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);
    d_hier_velocity_data_ops->subtract(f_u_idx, f_u_idx, d_f_idx);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f));
    IBTK_CHKERRQ(ierr);

    return ierr;

} // IBJacobianApply_velocity

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBPCApply_SAMRAI(PC pc, Vec x, Vec y)
{
    PetscErrorCode ierr = 1;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx);
    IBTK_CHKERRQ(ierr);
    IBImplicitStaggeredHierarchyIntegrator* ib_integrator = static_cast<IBImplicitStaggeredHierarchyIntegrator*>(ctx);
    if (ib_integrator->d_solve_for_position)
    {
        ierr = ib_integrator->IBPCApply_position(x, y);
    }
    else
    {
        ierr = ib_integrator->IBPCApply_velocity(x, y);
    }
    IBTK_CHKERRQ(ierr);
    return ierr;
} // IBPCApply_SAMRAI

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBPCApply_position(Vec x, Vec y)
{
    PetscErrorCode ierr;
    const double half_time = d_integrator_time + 0.5 * d_current_dt;

    Vec* component_x_vecs;
    Vec* component_y_vecs;
    ierr = VecMultiVecGetSubVecs(x, &component_x_vecs);
    IBTK_CHKERRQ(ierr);
    ierr = VecMultiVecGetSubVecs(y, &component_y_vecs);
    IBTK_CHKERRQ(ierr);

    Pointer<SAMRAIVectorReal<NDIM, double> > eul_x = PETScSAMRAIVectorReal::getSAMRAIVector(component_x_vecs[0]);
    Pointer<SAMRAIVectorReal<NDIM, double> > eul_y = PETScSAMRAIVectorReal::getSAMRAIVector(component_y_vecs[0]);

    Vec lag_x = component_x_vecs[1];
    Vec lag_y = component_y_vecs[1];

    // The full (nonlinear) system is:
    //
    //   L u(n+1) = S*F[X(n+1/2)] + f
    //   X(n+1) - X(n) = dt*U(n+1/2)
    //
    // where:
    //
    //   L = Eulerian operator (i.e. Stokes)
    //   F = Lagrangian force operator (potentially nonlinear)
    //   S = spreading operator
    //   J = interpolation operator = S^*
    //   X(n+1/2) = (X(n+1)+X(n))/2
    //   f = explicit right-hand side term
    //
    // For simplicity, only "lagged" S and J are considered, i.e., S and J are
    // not functions of the unknown X(n+1/2), but rather of some lagged
    // approximation to X(n+1/2).  This does not affect the stability of the
    // time stepping scheme, and the lagged values can be chosen so that the
    // overall scheme is second-order accurate.
    //
    // The linearized system is:
    //
    //   [L         -S*A/2] [u]
    //   [-dt*J/2   I     ] [X]
    //
    // where:
    //
    //   L = Eulerian operator
    //   A = dF/dX = Lagrangian operator
    //   S = spreading operator
    //   J = interpolation operator = S^*
    //
    // The Lagrangian Schur complement preconditioner is P = (4)*(3)*(2)*(1),
    // which is the inverse of the linearized system.
    //
    // (1) = [inv(L)  0]  ==>  [I         -inv(L)*S*A/2]
    //       [0       I]       [-dt*J/2   I            ]
    //
    // (2) = [I        0]  ==>  [I   -inv(L)*S*A/2      ]
    //       [dt*J/2   I]       [0   I-dt*J*inv(L)*S*A/4]
    //
    // Sc = Schur complement = I-dt*J*inv(L)*S*A/4
    //
    // (3) = [I   0      ]  ==>  [I   -inv(L)*S*A/2]
    //       [0   inv(Sc)]       [0   I            ]
    //
    // (4) = [I   inv(L)*S*A/2]  ==>  [I   0]
    //       [0   I           ]       [0   I]

    // Step 1: eul_y := inv(L)*eul_x
    eul_y->setToScalar(0.0);
    d_stokes_solver->setHomogeneousBc(true);
    d_stokes_solver->solveSystem(*eul_y, *eul_x);

    // Step 2: lag_y := lag_x + dt*J*eul_y/2
    d_hier_velocity_data_ops->scale(d_u_idx, -0.5, eul_y->getComponentDescriptorIndex(0));
    d_ib_implicit_ops->interpolateLinearizedVelocity(d_u_idx,
                                                     getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                                     getGhostfillRefineSchedules(d_object_name + "::u"),
                                                     half_time);
    d_ib_implicit_ops->computeLinearizedResidual(lag_x, lag_y);

    // Step 3: lag_y := inv(Sc)*lag_y
    ierr = KSPSolve(d_schur_solver, lag_y, lag_y);
    IBTK_CHKERRQ(ierr);

    // Step 4: eul_y := eul_y + inv(L)*S*A*lag_y/2
    d_ib_implicit_ops->computeLinearizedLagrangianForce(lag_y, half_time);
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0, /*interior_only*/ false);
    d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
    d_u_phys_bdry_op->setHomogeneousBc(true); // use homogeneous BCs to define spreading at physical boundaries
    d_ib_implicit_ops->spreadLinearizedForce(
        d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);
    d_u_scratch_vec->setToScalar(0.0);
    d_f_scratch_vec->setToScalar(0.0);
    d_hier_velocity_data_ops->copyData(d_f_scratch_vec->getComponentDescriptorIndex(0), d_f_idx);
    d_stokes_solver->setHomogeneousBc(true);
    d_stokes_solver->solveSystem(*d_u_scratch_vec, *d_f_scratch_vec);
    eul_y->add(eul_y, d_u_scratch_vec);

    // Ensure that PETSc sees that the state of the RHS vector has changed.
    // This is a nasty hack.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(component_y_vecs[0]));
    IBTK_CHKERRQ(ierr);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
    IBTK_CHKERRQ(ierr);
    return ierr;
} // IBPCApply_position

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::IBPCApply_velocity(Vec x, Vec y)
{
    Pointer<SAMRAIVectorReal<NDIM, double> > f_g = PETScSAMRAIVectorReal::getSAMRAIVector(x);
    Pointer<SAMRAIVectorReal<NDIM, double> > u_p = PETScSAMRAIVectorReal::getSAMRAIVector(y);
    Pointer<IBImplicitStaggeredStokesSolver> p_stokes_solver = d_stokes_solver;
#if !defined(NDEBUG)
    TBOX_ASSERT(p_stokes_solver);
#endif
    bool converged = p_stokes_solver->getStaggeredStokesFACPreconditioner()->solveSystem(*u_p, *f_g);
    PetscErrorCode ierr = !converged;
    return ierr;
} // IBPCApply_velocity

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::lagrangianSchurApply_SAMRAI(Mat A, Vec x, Vec y)
{
    PetscErrorCode ierr;
    void* ctx;
    ierr = MatShellGetContext(A, &ctx);
    IBTK_CHKERRQ(ierr);
    IBImplicitStaggeredHierarchyIntegrator* ib_integrator = static_cast<IBImplicitStaggeredHierarchyIntegrator*>(ctx);
    ierr = ib_integrator->lagrangianSchurApply(x, y);
    return ierr;
} // lagrangianSchurApply_SAMRAI

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::lagrangianSchurApply(Vec X, Vec Y)
{
    const double half_time = d_integrator_time + 0.5 * d_current_dt;

    // The Schur complement is: I-dt*J*inv(L)*S*A/4
    d_ib_implicit_ops->computeLinearizedLagrangianForce(X, half_time);
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0, /*interior_only*/ false);
    d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
    d_u_phys_bdry_op->setHomogeneousBc(true); // use homogeneous BCs to define spreading at physical boundaries
    d_ib_implicit_ops->spreadLinearizedForce(
        d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);
    d_u_scratch_vec->setToScalar(0.0);
    d_hier_velocity_data_ops->copyData(d_f_scratch_vec->getComponentDescriptorIndex(0), d_f_idx);
    d_stokes_solver->setHomogeneousBc(true);
    d_stokes_solver->solveSystem(*d_u_scratch_vec, *d_f_scratch_vec);
    d_hier_velocity_data_ops->scale(d_u_idx, 0.5, d_u_scratch_vec->getComponentDescriptorIndex(0));
    d_ib_implicit_ops->interpolateLinearizedVelocity(d_u_idx,
                                                     getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                                     getGhostfillRefineSchedules(d_object_name + "::u"),
                                                     half_time);
    d_ib_implicit_ops->computeLinearizedResidual(X, Y);
    return 0;
} // lagrangianSchurApply

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
