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

#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/pdat/CellData.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/math/PatchCellDataOpsReal.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/math/PatchSideDataOpsReal.h"
#include "SAMRAI/solv/PoissonSpecifications.h"
#include "SAMRAI/solv/SAMRAIVectorReal.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBImplicitStaggeredHierarchyIntegrator.h"
#include "ibamr/IBImplicitStrategy.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/INSStaggeredHierarchyIntegrator.h"
#include "ibamr/StaggeredStokesOperator.h"
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
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/PIO.h"

#include "SAMRAI/tbox/RestartManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

namespace IBAMR
{
} // namespace IBAMR
namespace SAMRAI
{
namespace hier
{

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
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitStaggeredHierarchyIntegrator::IBImplicitStaggeredHierarchyIntegrator(
    const std::string& object_name,
    const boost::shared_ptr<Database>& input_db,
    const boost::shared_ptr<IBImplicitStrategy>& ib_implicit_ops,
    const boost::shared_ptr<INSStaggeredHierarchyIntegrator>& ins_hier_integrator,
    bool register_for_restart)
    : IBHierarchyIntegrator(object_name, input_db, ib_implicit_ops, ins_hier_integrator, register_for_restart),
      d_ib_implicit_ops(ib_implicit_ops)
{
    // Setup IB ops object to use "fixed" Lagrangian-Eulerian coupling
    // operators.
    d_ib_implicit_ops->setUseFixedLEOperators(true);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;
}

IBImplicitStaggeredHierarchyIntegrator::~IBImplicitStaggeredHierarchyIntegrator()
{
    // intentionally blank
    return;
}

void IBImplicitStaggeredHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
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
        auto level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_u_idx, current_time);
        level->allocatePatchData(d_f_idx, current_time);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
    }

    // Initialize IB data.
    d_ib_implicit_ops->preprocessIntegrateData(current_time, new_time, num_cycles);

    // Initialize the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    if (ins_num_cycles != d_current_num_cycles && d_current_num_cycles != 1)
    {
        TBOX_ERROR(d_object_name << "::preprocessIntegrateHierarchy():\n"
                                 << "  attempting to perform " << d_current_num_cycles
                                 << " cycles of fixed point iteration.\n"
                                 << "  number of cycles required by Navier-Stokes solver = " << ins_num_cycles << ".\n"
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
}

void IBImplicitStaggeredHierarchyIntegrator::integrateHierarchy(const double current_time,
                                                                const double new_time,
                                                                const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);

    auto p_ins_hier_integrator = BOOST_CAST<INSStaggeredHierarchyIntegrator>(d_ins_hier_integrator);

    PetscErrorCode ierr;
    int n_local;

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    //  const double half_time = current_time+0.5*(new_time-current_time);

    auto var_db = VariableDatabase::getDatabase();
    auto current_ctx = p_ins_hier_integrator->getCurrentContext();
    auto scratch_ctx = p_ins_hier_integrator->getScratchContext();
    auto new_ctx = p_ins_hier_integrator->getNewContext();

    const int wgt_cc_idx = d_hier_math_ops->getCellWeightPatchDescriptorIndex();
    const int wgt_sc_idx = d_hier_math_ops->getSideWeightPatchDescriptorIndex();

    auto u_var = p_ins_hier_integrator->getVelocityVariable();
    //  const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, current_ctx);
    const int u_scratch_idx = var_db->mapVariableAndContextToIndex(u_var, scratch_ctx);
    //  const int u_new_idx     = var_db->mapVariableAndContextToIndex(u_var, new_ctx);

    auto p_var = p_ins_hier_integrator->getPressureVariable();
    //  const int p_current_idx = var_db->mapVariableAndContextToIndex(p_var, current_ctx);
    const int p_scratch_idx = var_db->mapVariableAndContextToIndex(p_var, scratch_ctx);
    //  const int p_new_idx     = var_db->mapVariableAndContextToIndex(p_var, new_ctx);

    // Skip all cycles in the INS solver --- we advance the state data here.
    p_ins_hier_integrator->skipCycle(current_time, new_time, cycle_num);

    // Setup Eulerian vectors used in solving the implicit IB equations.
    auto eul_sol_vec = boost::make_shared<SAMRAIVectorReal<double> >(d_object_name + "::eulerian_sol_vec", d_hierarchy,
                                                                     coarsest_ln, finest_ln);
    eul_sol_vec->addComponent(u_var, u_scratch_idx, wgt_sc_idx, d_hier_velocity_data_ops);
    eul_sol_vec->addComponent(p_var, p_scratch_idx, wgt_cc_idx, d_hier_pressure_data_ops);

    auto eul_rhs_vec = eul_sol_vec->cloneVector(d_object_name + "::eulerian_rhs_vec");
    eul_rhs_vec->allocateVectorData(current_time);

    d_u_scratch_vec = eul_sol_vec->cloneVector(d_object_name + "::u_scratch_vec");
    d_f_scratch_vec = eul_rhs_vec->cloneVector(d_object_name + "::f_scratch_vec");
    d_u_scratch_vec->allocateVectorData(current_time);
    d_f_scratch_vec->allocateVectorData(current_time);

    p_ins_hier_integrator->setupSolverVectors(eul_sol_vec, eul_rhs_vec, current_time, new_time, cycle_num);

    d_stokes_solver = p_ins_hier_integrator->getStokesSolver();
    auto p_stokes_solver = BOOST_CAST<KrylovLinearSolver>(d_stokes_solver);
    d_stokes_op = BOOST_CAST<StaggeredStokesOperator>(p_stokes_solver->getOperator());

    // Setup Lagrangian vectors used in solving the implicit IB equations.
    Vec lag_sol_petsc_vec, lag_rhs_petsc_vec;
    d_ib_implicit_ops->createSolverVecs(lag_sol_petsc_vec, lag_rhs_petsc_vec);
    d_ib_implicit_ops->setupSolverVecs(lag_sol_petsc_vec, lag_rhs_petsc_vec);

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
    ierr = SNESSetFunction(snes, composite_res_petsc_vec, compositeIBFunction_SAMRAI, this);
    IBTK_CHKERRQ(ierr);
    ierr = SNESSetOptionsPrefix(snes, "ib_");
    IBTK_CHKERRQ(ierr);

    Mat jac;
    ierr = VecGetLocalSize(composite_sol_petsc_vec, &n_local);
    IBTK_CHKERRQ(ierr);
    ierr = MatCreateShell(PETSC_COMM_WORLD, n_local, n_local, PETSC_DETERMINE, PETSC_DETERMINE, this, &jac);
    IBTK_CHKERRQ(ierr);
    ierr = MatShellSetOperation(jac, MATOP_MULT, reinterpret_cast<void (*)(void)>(compositeIBJacobianApply_SAMRAI));
    IBTK_CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes, jac, jac, compositeIBJacobianSetup_SAMRAI, this);
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
    ierr = PCShellSetApply(snes_pc, compositeIBPCApply_SAMRAI);
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
    p_ins_hier_integrator->resetSolverVectors(eul_sol_vec, eul_rhs_vec, current_time, new_time, cycle_num);

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_ib_implicit_ops->setUpdatedPosition(lag_sol_petsc_vec);
#if 0
    d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
    if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): interpolating Eulerian velocity to the Lagrangian mesh\n";
    d_ib_implicit_ops->interpolateVelocity(d_u_idx, getCoarsenSchedules(d_object_name+"::u::CONSERVATIVE_COARSEN"), getGhostfillRefineSchedules(d_object_name+"::u"), half_time);

    // Compute the final value of the updated positions of the Lagrangian
    // structure.
    d_ib_implicit_ops->midpointStep(current_time, new_time);
#endif

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
}

void IBImplicitStaggeredHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
                                                                           const double new_time,
                                                                           const bool skip_synchronize_new_state_data,
                                                                           const int num_cycles)
{
    IBHierarchyIntegrator::postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data,
                                                         num_cycles);

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();
    const double dt = new_time - current_time;
    auto var_db = VariableDatabase::getDatabase();
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): interpolating Eulerian "
                                 "velocity to the Lagrangian mesh\n";
    d_ib_implicit_ops->interpolateVelocity(d_u_idx, getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                           getGhostfillRefineSchedules(d_object_name + "::u"), new_time);

    // Synchronize new state data.
    if (!skip_synchronize_new_state_data)
    {
        if (d_enable_logging)
            plog << d_object_name << "::postprocessIntegrateHierarchy(): synchronizing updated data\n";
        synchronizeHierarchyData(NEW_DATA);
    }

    // Determine the CFL number.
    double cfl_max = 0.0;
    PatchCellDataOpsReal<double> patch_cc_ops;
    PatchSideDataOpsReal<double> patch_sc_ops;
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            auto patch = *p;
            const Box& patch_box = patch->getBox();
            const auto pgeom = BOOST_CAST<CartesianPatchGeometry>(patch->getPatchGeometry());
            const double* const dx = pgeom->getDx();
            const double dx_min = *(std::min_element(dx, dx + NDIM));
            auto u_cc_new_data = boost::dynamic_pointer_cast<CellData<double> >(patch->getPatchData(u_new_idx));
            auto u_sc_new_data = boost::dynamic_pointer_cast<SideData<double> >(patch->getPatchData(u_new_idx));
            double u_max = 0.0;
            if (u_cc_new_data) u_max = patch_cc_ops.maxNorm(u_cc_new_data, patch_box);
            if (u_sc_new_data) u_max = patch_sc_ops.maxNorm(u_sc_new_data, patch_box);
            cfl_max = std::max(cfl_max, u_max * dt / dx_min);
        }
    }
    tbox::SAMRAI_MPI comm(MPI_COMM_WORLD);
    comm.AllReduce(&cfl_max, 1, MPI_MAX);
    d_regrid_cfl_estimate += cfl_max;
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): CFL number = " << cfl_max << "\n";
    if (d_enable_logging)
        plog << d_object_name << "::postprocessIntegrateHierarchy(): estimated upper bound on IB "
                                 "point displacement since last regrid = " << d_regrid_cfl_estimate << "\n";

    // Deallocate the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    d_ins_hier_integrator->postprocessIntegrateHierarchy(current_time, new_time, skip_synchronize_new_state_data,
                                                         ins_num_cycles);

    // Deallocate IB data.
    d_ib_implicit_ops->postprocessIntegrateData(current_time, new_time, num_cycles);

    // Deallocate Eulerian scratch data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        auto level = d_hierarchy->getPatchLevel(ln);
        level->deallocatePatchData(d_u_idx);
        level->deallocatePatchData(d_f_idx);
        level->deallocatePatchData(d_scratch_data);
        level->deallocatePatchData(d_new_data);
    }

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(current_time, new_time, skip_synchronize_new_state_data,
                                                     num_cycles);
    return;
}

void IBImplicitStaggeredHierarchyIntegrator::initializeHierarchyIntegrator(
    const boost::shared_ptr<PatchHierarchy>& hierarchy,
    const boost::shared_ptr<GriddingAlgorithm>& gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Finish initializing the hierarchy integrator.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);
}

int IBImplicitStaggeredHierarchyIntegrator::getNumberOfCycles() const
{
    return d_ins_hier_integrator->getNumberOfCycles();
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void IBImplicitStaggeredHierarchyIntegrator::putToRestartSpecialized(const boost::shared_ptr<Database>& db) const
{
    IBHierarchyIntegrator::putToRestartSpecialized(db);
    db->putInteger("IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION",
                   IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION);
    return;
}

/////////////////////////////// PRIVATE //////////////////////////////////////

void IBImplicitStaggeredHierarchyIntegrator::getFromRestart()
{
    auto restart_db = RestartManager::getManager()->getRootDatabase();
    boost::shared_ptr<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
}

PetscErrorCode IBImplicitStaggeredHierarchyIntegrator::compositeIBFunction_SAMRAI(SNES snes, Vec x, Vec f, void* ctx)
{
    IBImplicitStaggeredHierarchyIntegrator* ib_integrator = static_cast<IBImplicitStaggeredHierarchyIntegrator*>(ctx);
    return ib_integrator->compositeIBFunction(snes, x, f);
}

PetscErrorCode IBImplicitStaggeredHierarchyIntegrator::compositeIBFunction(SNES /*snes*/, Vec x, Vec f)
{
    PetscErrorCode ierr;
    const double half_time = d_integrator_time + 0.5 * d_current_dt;

    Vec* component_sol_vecs;
    Vec* component_rhs_vecs;
    ierr = VecMultiVecGetSubVecs(x, &component_sol_vecs);
    CHKERRQ(ierr);
    ierr = VecMultiVecGetSubVecs(f, &component_rhs_vecs);
    CHKERRQ(ierr);

    auto u = PETScSAMRAIVectorReal::getSAMRAIVector(component_sol_vecs[0]);
    auto f_u = PETScSAMRAIVectorReal::getSAMRAIVector(component_rhs_vecs[0]);

    auto var_db = VariableDatabase::getDatabase();
    auto current_ctx = d_ins_hier_integrator->getCurrentContext();
    //  auto scratch_ctx = d_ins_hier_integrator->getScratchContext();
    //  auto new_ctx     = d_ins_hier_integrator->getNewContext();

    auto u_var = d_ins_hier_integrator->getVelocityVariable();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(u_var, current_ctx);
    const int u_new_idx = u->getComponentDescriptorIndex(0);
    const int f_u_idx = f_u->getComponentDescriptorIndex(0);

    Vec X = component_sol_vecs[1];
    Vec R = component_rhs_vecs[1];

    // Evaluate the Eulerian terms.
    d_stokes_op->setHomogeneousBc(false);
    d_stokes_op->apply(*u, *f_u);

    d_ib_implicit_ops->setUpdatedPosition(X);
    d_ib_implicit_ops->computeLagrangianForce(half_time);
    if (d_enable_logging)
        plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
    d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
    d_ib_implicit_ops->spreadForce(d_f_idx, d_u_phys_bdry_op.get(), getProlongRefineSchedules(d_object_name + "::f"),
                                   half_time);
    d_hier_velocity_data_ops->subtract(f_u_idx, f_u_idx, d_f_idx);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(component_rhs_vecs[0]));
    CHKERRQ(ierr);

    // Evaluate the Lagrangian terms.
    d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
    d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
    d_ib_implicit_ops->interpolateVelocity(d_u_idx, getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                           getGhostfillRefineSchedules(d_object_name + "::u"), half_time);
    d_ib_implicit_ops->computeResidual(R);

    // Ensure that PETSc sees that the state of the RHS vector has changed.
    // This is a nasty hack.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f));
    CHKERRQ(ierr);
    return ierr;
}

PetscErrorCode
IBImplicitStaggeredHierarchyIntegrator::compositeIBJacobianSetup_SAMRAI(SNES snes, Vec x, Mat A, Mat B, void* ctx)
{
    IBImplicitStaggeredHierarchyIntegrator* ib_integrator = static_cast<IBImplicitStaggeredHierarchyIntegrator*>(ctx);
    return ib_integrator->compositeIBJacobianSetup(snes, x, A, B);
}

PetscErrorCode IBImplicitStaggeredHierarchyIntegrator::compositeIBJacobianSetup(SNES /*snes*/, Vec x, Mat A, Mat /*B*/)
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
    return 0;
}

PetscErrorCode IBImplicitStaggeredHierarchyIntegrator::compositeIBJacobianApply_SAMRAI(Mat A, Vec x, Vec f)
{
    PetscErrorCode ierr;
    void* ctx;
    ierr = MatShellGetContext(A, &ctx);
    IBTK_CHKERRQ(ierr);
    IBImplicitStaggeredHierarchyIntegrator* ib_integrator = static_cast<IBImplicitStaggeredHierarchyIntegrator*>(ctx);
    return ib_integrator->compositeIBJacobianApply(x, f);
}

PetscErrorCode IBImplicitStaggeredHierarchyIntegrator::compositeIBJacobianApply(Vec x, Vec f)
{
    PetscErrorCode ierr;
    const double half_time = d_integrator_time + 0.5 * d_current_dt;

    Vec* component_sol_vecs;
    Vec* component_rhs_vecs;
    ierr = VecMultiVecGetSubVecs(x, &component_sol_vecs);
    IBTK_CHKERRQ(ierr);
    ierr = VecMultiVecGetSubVecs(f, &component_rhs_vecs);
    IBTK_CHKERRQ(ierr);

    auto u = PETScSAMRAIVectorReal::getSAMRAIVector(component_sol_vecs[0]);
    auto f_u = PETScSAMRAIVectorReal::getSAMRAIVector(component_rhs_vecs[0]);

    auto u_var = d_ins_hier_integrator->getVelocityVariable();
    const int u_idx = u->getComponentDescriptorIndex(0);
    const int f_u_idx = f_u->getComponentDescriptorIndex(0);

    Vec X = component_sol_vecs[1];
    Vec R = component_rhs_vecs[1];

    // Evaluate the Eulerian terms.
    d_stokes_op->setHomogeneousBc(true);
    d_stokes_op->apply(*u, *f_u);

    d_ib_implicit_ops->computeLinearizedLagrangianForce(X, half_time);
    if (d_enable_logging)
        plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
    d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
    d_ib_implicit_ops->spreadLinearizedForce(d_f_idx, d_u_phys_bdry_op.get(),
                                             getProlongRefineSchedules(d_object_name + "::f"), half_time);
    d_hier_velocity_data_ops->subtract(f_u_idx, f_u_idx, d_f_idx);
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(component_rhs_vecs[0]));
    IBTK_CHKERRQ(ierr);

    // Evaluate the Lagrangian terms.
    d_hier_velocity_data_ops->scale(d_u_idx, 0.5, u_idx);
    d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
    d_ib_implicit_ops->interpolateLinearizedVelocity(d_u_idx,
                                                     getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                                     getGhostfillRefineSchedules(d_object_name + "::u"), half_time);
    d_ib_implicit_ops->computeLinearizedResidual(X, R);

    // Ensure that PETSc sees that the state of the RHS vector has changed.
    // This is a nasty hack.
    ierr = PetscObjectStateIncrease(reinterpret_cast<PetscObject>(f));
    IBTK_CHKERRQ(ierr);
    return ierr;
}

PetscErrorCode IBImplicitStaggeredHierarchyIntegrator::compositeIBPCApply_SAMRAI(PC pc, Vec x, Vec y)
{
    PetscErrorCode ierr;
    void* ctx;
    ierr = PCShellGetContext(pc, &ctx);
    IBTK_CHKERRQ(ierr);
    IBImplicitStaggeredHierarchyIntegrator* ib_integrator = static_cast<IBImplicitStaggeredHierarchyIntegrator*>(ctx);
    ierr = ib_integrator->compositeIBPCApply(x, y);
    IBTK_CHKERRQ(ierr);
    return ierr;
}

PetscErrorCode IBImplicitStaggeredHierarchyIntegrator::compositeIBPCApply(Vec x, Vec y)
{
    PetscErrorCode ierr;
    const double half_time = d_integrator_time + 0.5 * d_current_dt;

    Vec* component_x_vecs;
    Vec* component_y_vecs;
    ierr = VecMultiVecGetSubVecs(x, &component_x_vecs);
    IBTK_CHKERRQ(ierr);
    ierr = VecMultiVecGetSubVecs(y, &component_y_vecs);
    IBTK_CHKERRQ(ierr);

    auto eul_x = PETScSAMRAIVectorReal::getSAMRAIVector(component_x_vecs[0]);
    auto eul_y = PETScSAMRAIVectorReal::getSAMRAIVector(component_y_vecs[0]);

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
                                                     getGhostfillRefineSchedules(d_object_name + "::u"), half_time);
    d_ib_implicit_ops->computeLinearizedResidual(lag_x, lag_y);

    // Step 3: lag_y := inv(Sc)*lag_y
    ierr = KSPSolve(d_schur_solver, lag_y, lag_y);
    IBTK_CHKERRQ(ierr);

    // Step 4: eul_y := eul_y + inv(L)*S*A*lag_y/2
    d_ib_implicit_ops->computeLinearizedLagrangianForce(lag_y, half_time);
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
    d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
    d_ib_implicit_ops->spreadLinearizedForce(d_f_idx, d_u_phys_bdry_op.get(),
                                             getProlongRefineSchedules(d_object_name + "::f"), half_time);
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
}

PetscErrorCode IBImplicitStaggeredHierarchyIntegrator::lagrangianSchurApply_SAMRAI(Mat A, Vec x, Vec y)
{
    PetscErrorCode ierr;
    void* ctx;
    ierr = MatShellGetContext(A, &ctx);
    IBTK_CHKERRQ(ierr);
    IBImplicitStaggeredHierarchyIntegrator* ib_integrator = static_cast<IBImplicitStaggeredHierarchyIntegrator*>(ctx);
    ierr = ib_integrator->lagrangianSchurApply(x, y);
    return ierr;
}

PetscErrorCode IBImplicitStaggeredHierarchyIntegrator::lagrangianSchurApply(Vec X, Vec Y)
{
    const double half_time = d_integrator_time + 0.5 * d_current_dt;

    // The Schur complement is: I-dt*J*inv(L)*S*A/4
    d_ib_implicit_ops->computeLinearizedLagrangianForce(X, half_time);
    d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
    d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
    d_ib_implicit_ops->spreadLinearizedForce(d_f_idx, d_u_phys_bdry_op.get(),
                                             getProlongRefineSchedules(d_object_name + "::f"), half_time);
    d_u_scratch_vec->setToScalar(0.0);
    d_hier_velocity_data_ops->copyData(d_f_scratch_vec->getComponentDescriptorIndex(0), d_f_idx);
#if 0
    d_stokes_solver->setHomogeneousBc(true);
    d_stokes_solver->solveSystem(*d_u_scratch_vec, *d_f_scratch_vec);
#else
    d_u_scratch_vec->scale(1.0 / d_stokes_op->getVelocityPoissonSpecifications().getCConstant(), d_f_scratch_vec);
#endif
    d_hier_velocity_data_ops->scale(d_u_idx, 0.5, d_u_scratch_vec->getComponentDescriptorIndex(0));
    d_ib_implicit_ops->interpolateLinearizedVelocity(d_u_idx,
                                                     getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                                     getGhostfillRefineSchedules(d_object_name + "::u"), half_time);
    d_ib_implicit_ops->computeLinearizedResidual(X, Y);
    return 0;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
