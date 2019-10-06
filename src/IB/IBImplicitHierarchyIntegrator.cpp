// Filename: IBImplicitHierarchyIntegrator.cpp
// Created on 07 Apr 2012 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
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

#include "ibamr/IBHierarchyIntegrator.h"
#include "ibamr/IBImplicitHierarchyIntegrator.h"
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
#include "ibtk/PETScSAMRAIVectorReal.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_enums.h"

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
#include "tbox/Database.h"
#include "tbox/PIO.h"
#include "tbox/Pointer.h"
#include "tbox/RestartManager.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"

#include "petscsnes.h"
#include "petscvec.h"

#include <algorithm>
#include <ostream>
#include <string>

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
// Version of IBImplicitHierarchyIntegrator restart file data.
static const int IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION = 1;
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBImplicitHierarchyIntegrator::IBImplicitHierarchyIntegrator(
    const std::string& object_name,
    Pointer<Database> input_db,
    Pointer<IBImplicitStrategy> ib_implicit_ops,
    Pointer<INSStaggeredHierarchyIntegrator> ins_hier_integrator,
    bool register_for_restart)
    : IBHierarchyIntegrator(object_name, input_db, ib_implicit_ops, ins_hier_integrator, register_for_restart),
      d_ib_implicit_ops(ib_implicit_ops)
{
    // Set default configuration options.
    d_use_structure_predictor = true;
    d_use_fixed_LE_operators = false;

    // Set options from input.
    if (input_db)
    {
        if (input_db->keyExists("use_structure_predictor"))
            d_use_structure_predictor = input_db->getBool("use_structure_predictor");
        if (input_db->keyExists("use_fixed_LE_operators"))
            d_use_fixed_LE_operators = input_db->getBool("use_fixed_LE_operators");
        if (input_db->keyExists("implicit_algorithm")) d_implicit_algorithm = input_db->getString("implicit_algorithm");
    }
    d_ib_implicit_ops->setUseFixedLEOperators(d_use_fixed_LE_operators);

    // Initialize object with data read from the input and restart databases.
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;
} // IBImplicitHierarchyIntegrator

void
IBImplicitHierarchyIntegrator::preprocessIntegrateHierarchy(const double current_time,
                                                            const double new_time,
                                                            const int num_cycles)
{
    IBHierarchyIntegrator::preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    switch (d_time_stepping_type)
    {
    case BACKWARD_EULER:
    case TRAPEZOIDAL_RULE:
    case MIDPOINT_RULE:
        break;
    default:
        TBOX_ERROR("IBImplicitHierarchyIntegrator::preprocessIntegrateHierarchy(): time_stepping_type = "
                   << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                   << "  only supported time_stepping_types are:\n"
                   << "    " << enum_to_string<TimeSteppingType>(BACKWARD_EULER) << "\n"
                   << "    " << enum_to_string<TimeSteppingType>(TRAPEZOIDAL_RULE) << "\n"
                   << "    " << enum_to_string<TimeSteppingType>(MIDPOINT_RULE) << "\n");
    }

    const int coarsest_ln = 0;
    const int finest_ln = d_hierarchy->getFinestLevelNumber();

    // Allocate Eulerian scratch and new data.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(d_u_idx, current_time);
        level->allocatePatchData(d_f_idx, current_time);
        if (d_f_current_idx != -1) level->allocatePatchData(d_f_current_idx, current_time);
        level->allocatePatchData(d_scratch_data, current_time);
        level->allocatePatchData(d_new_data, new_time);
    }

    // Initialize IB data.
    d_ib_implicit_ops->preprocessIntegrateData(current_time, new_time, num_cycles);

    // Initialize the fluid solver.
    const int ins_num_cycles = d_ins_hier_integrator->getNumberOfCycles();
    d_ins_hier_integrator->preprocessIntegrateHierarchy(current_time, new_time, ins_num_cycles);
    d_ins_hier_integrator->setSkipEnforceNumCycles();

    // Compute the Lagrangian forces and spread them to the Eulerian grid.
    switch (d_time_stepping_type)
    {
    case TRAPEZOIDAL_RULE:
        if (d_enable_logging) plog << d_object_name << "::preprocessIntegrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(current_time);
        if (d_enable_logging)
            plog << d_object_name
                 << "::preprocessIntegrateHierarchy(): spreading Lagrangian force "
                    "to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), current_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_hier_velocity_data_ops->copyData(d_f_current_idx, d_f_idx);
        break;
    case BACKWARD_EULER:
    case MIDPOINT_RULE:
        // intentionally blank
        break;
    default:
        TBOX_ERROR(
            d_object_name << "::preprocessIntegrateHierarchy():\n"
                          << "  unsupported time stepping type: "
                          << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                          << "  supported time stepping types are: BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute an initial prediction of the updated positions of the Lagrangian
    // structure.
    //
    // NOTE: The velocity should already have been interpolated to the
    // curvilinear mesh and should not need to be re-interpolated.
    if (d_use_structure_predictor)
    {
        if (d_enable_logging)
            plog << d_object_name << "::preprocessIntegrateHierarchy(): performing Lagrangian forward Euler step\n";
        d_ib_implicit_ops->forwardEulerStep(current_time, new_time);
    }

    // Execute any registered callbacks.
    executePreprocessIntegrateHierarchyCallbackFcns(current_time, new_time, num_cycles);
    return;
} // preprocessIntegrateHierarchy

void
IBImplicitHierarchyIntegrator::integrateHierarchy(const double current_time, const double new_time, const int cycle_num)
{
    IBHierarchyIntegrator::integrateHierarchy(current_time, new_time, cycle_num);
    d_current_time = current_time;
    d_new_time = new_time;
    d_cycle_num = cycle_num;

    if (d_use_fixed_LE_operators) d_ib_implicit_ops->updateFixedLEOperators();

    if (d_implicit_algorithm == "AITKEN" || d_implicit_algorithm == "ANDERSON")
    {
        // Use the crossed secant method (Aitken extrapolation) or the alternate secant method (Anderson extrapolation
        // with M = 1).
        //
        // Ref: Isabelle Ramière, Thomas Helfer. Iterative residual-based vector methods to accelerate fixed point
        // iterations. Computers and Mathematics with Applications, Elsevier, 2015, 70, pp. 2210 - 2226.
        //
        // See Pg. 21 of https://hal-cea.archives-ouvertes.fr/cea-01403292/document.
        PetscErrorCode ierr;
        Vec X_nm1, X_n, X_np1, GX_nm1, GX_n, DX_nm1, DX_n, DDX_n, DGX_n;
        d_ib_implicit_ops->createSolutionVec(&X_nm1);
        d_ib_implicit_ops->createSolutionVec(&X_n);
        d_ib_implicit_ops->createSolutionVec(&X_np1);
        d_ib_implicit_ops->createSolutionVec(&GX_nm1);
        d_ib_implicit_ops->createSolutionVec(&GX_n);
        d_ib_implicit_ops->createSolutionVec(&DX_nm1);
        d_ib_implicit_ops->createSolutionVec(&DX_n);
        d_ib_implicit_ops->createSolutionVec(&DDX_n);
        d_ib_implicit_ops->createSolutionVec(&DGX_n);

        // Initialize X_{n-1}:
        d_ib_implicit_ops->getUpdatedPosition(X_nm1);

        double D0;
        ierr = VecNorm(X_nm1, NORM_2, &D0);
        IBTK_CHKERRQ(ierr);

        int k = 0, n_solves = 0;
        double tol = 1.0e-5;
        int max_its = 100;
        bool converged = false;

        // Initialize X_{n} = G(X_{n-1}):
        iterateSolution(X_nm1);
        d_ib_implicit_ops->getUpdatedPosition(GX_nm1);
        ++n_solves;
        ierr = VecCopy(GX_nm1, X_n);
        IBTK_CHKERRQ(ierr);

        // Initialize DX_{n-1} := G(X_{n-1}) - X_{n-1}:
        ierr = VecWAXPY(DX_nm1, -1.0, X_nm1, GX_nm1);
        IBTK_CHKERRQ(ierr);

        for (k = 0; k < max_its && !converged; ++k)
        {
            // G(X_{n}):
            iterateSolution(X_n);
            d_ib_implicit_ops->getUpdatedPosition(GX_n);
            ++n_solves;

            // DX_{n} := G(X_{n}) - X_{n}:
            ierr = VecWAXPY(DX_n, -1.0, X_n, GX_n);
            IBTK_CHKERRQ(ierr);

            // DGX_{n} := G(X_{n}) - G(X_{n-1}):
            ierr = VecWAXPY(DGX_n, -1.0, GX_nm1, GX_n);
            IBTK_CHKERRQ(ierr);

            // DDX_{n} := DX_{n} - DX_{n-1}:
            ierr = VecWAXPY(DDX_n, -1.0, DX_nm1, DX_n);
            IBTK_CHKERRQ(ierr);

            if (d_implicit_algorithm == "AITKEN")
            {
                // X_{n+1} := G(X_{n}) - Q DX_{n}
                //
                // with Q = dot(G(X_{n}) - G(X_{n-1}, DX_{n} - DX_{n-1})/||DX_{n} - DX_{n-1}||^2
                //        = dot(DGX_{n}, DDX_{n})/||DDX_{n}||^2
                double Q1;
                ierr = VecDot(DGX_n, DDX_n, &Q1);
                IBTK_CHKERRQ(ierr);
                double Q2;
                ierr = VecDot(DDX_n, DDX_n, &Q2);
                IBTK_CHKERRQ(ierr);

                const double Q = Q1 / Q2;
                ierr = VecWAXPY(X_np1, -Q, DX_n, GX_n);
                IBTK_CHKERRQ(ierr);
            }
            else if (d_implicit_algorithm == "ANDERSON")
            {
                // X_{n+1} := G(X_{n}) - Q (G(X_{n}) - G(X_{n-1}))
                //          = G(X_{n}) - Q DGX_{n}
                //
                // with Q = dot(DX_{n} - DX_{n-1}, DX_{n})/||DX_{n} - DX_{n-1}||^2
                //        = dot(DDX_{n}, DX_{n})/||DDX_{n}||^2
                double Q1;
                ierr = VecDot(DDX_n, DX_n, &Q1);
                IBTK_CHKERRQ(ierr);
                double Q2;
                ierr = VecDot(DDX_n, DDX_n, &Q2);
                IBTK_CHKERRQ(ierr);

                const double Q = Q1 / Q2;
                ierr = VecWAXPY(X_np1, -Q, DGX_n, GX_n);
                IBTK_CHKERRQ(ierr);
            }

            // Check to see if ||G(X_{n}) - X_{n}||_2 is sufficiently small to declare convergence.
            double D;
            ierr = VecNorm(DX_n, NORM_2, &D);
            IBTK_CHKERRQ(ierr);
            if (D < tol * D0)
            {
                converged = true;
            }
            else
            {
                // If we haven't converged, get an updated solution, reset values and repeat.
                ierr = VecCopy(X_n, X_nm1);
                IBTK_CHKERRQ(ierr);
                ierr = VecCopy(X_np1, X_n);
                IBTK_CHKERRQ(ierr);
                ierr = VecCopy(GX_n, GX_nm1);
                IBTK_CHKERRQ(ierr);
                ierr = VecCopy(DX_n, DX_nm1);
                IBTK_CHKERRQ(ierr);
            }
        }

        if (!converged)
        {
            pout << "WARNING: " << d_implicit_algorithm << " extrapolation DIVERGED\n";
        }
        plog << d_object_name << "::integrateHierarchy():\n"
             << "  " << d_implicit_algorithm << " extrapolation " << (converged ? "CONVERGED" : "DIVERGED") << "\n"
             << "  requiring in " << k << " iterations and " << n_solves << " total Navier-Stokes solves\n";

        ierr = VecDestroy(&X_nm1);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&X_n);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&X_np1);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&GX_nm1);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&GX_n);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&DX_nm1);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&DX_n);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&DDX_n);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&DGX_n);
        IBTK_CHKERRQ(ierr);
    }
    else if (d_implicit_algorithm == "SNES")
    {
        // Setup Lagrangian vectors used in solving the implicit IB equations.
        PetscErrorCode ierr;
        Vec X, R, R_work;
        d_ib_implicit_ops->createSolverVecs(&X, &R);
        d_ib_implicit_ops->setupSolverVecs(X, R);
        ierr = VecDuplicate(R, &R_work);
        IBTK_CHKERRQ(ierr);

        // Solve the implicit IB equations.
        d_ins_cycle_num = 0;

        SNES snes;
        ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
        IBTK_CHKERRQ(ierr);
        ierr = SNESSetFunction(snes, R_work, IBFunction_SAMRAI, this);
        IBTK_CHKERRQ(ierr);
        ierr = SNESSetOptionsPrefix(snes, "ib_");
        IBTK_CHKERRQ(ierr);
        ierr = SNESSetFromOptions(snes);
        IBTK_CHKERRQ(ierr);
        ierr = SNESSolve(snes, R, X);
        IBTK_CHKERRQ(ierr);
        ierr = SNESDestroy(&snes);
        IBTK_CHKERRQ(ierr);

        // Ensure that the INS variables are consistent with the final structure configuration.
        iterateSolution(X);

        // Deallocate temporary data.
        ierr = VecDestroy(&X);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&R);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&R_work);
        IBTK_CHKERRQ(ierr);
    }
    else
    {
        TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                 << "  unsupported solver algorithm: " << d_implicit_algorithm << "\n"
                                 << "  supported time stepping types are: AITKEN, ANDERSON, SNES\n");
    }

    // Execute any registered callbacks.
    executeIntegrateHierarchyCallbackFcns(current_time, new_time, cycle_num);
} // integrateHierarchy

void
IBImplicitHierarchyIntegrator::postprocessIntegrateHierarchy(const double current_time,
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
        plog << d_object_name
             << "::postprocessIntegrateHierarchy(): interpolating Eulerian "
                "velocity to the Lagrangian mesh\n";
    d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
    d_u_phys_bdry_op->setHomogeneousBc(false);
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
        plog << d_object_name
             << "::postprocessIntegrateHierarchy(): estimated upper bound on IB "
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
        if (d_f_current_idx != -1) level->deallocatePatchData(d_f_current_idx);
    }

    // Execute any registered callbacks.
    executePostprocessIntegrateHierarchyCallbackFcns(
        current_time, new_time, skip_synchronize_new_state_data, num_cycles);
    return;
} // postprocessIntegrateHierarchy

void
IBImplicitHierarchyIntegrator::initializeHierarchyIntegrator(Pointer<PatchHierarchy<NDIM> > hierarchy,
                                                             Pointer<GriddingAlgorithm<NDIM> > gridding_alg)
{
    if (d_integrator_is_initialized) return;

    // Register body force function with INSHierarchyIntegrator
    d_ins_hier_integrator->registerBodyForceFunction(new IBEulerianForceFunction(this));

    // Finish initializing the hierarchy integrator.
    IBHierarchyIntegrator::initializeHierarchyIntegrator(hierarchy, gridding_alg);
    return;
} // initializeHierarchyIntegrator

int
IBImplicitHierarchyIntegrator::getNumberOfCycles() const
{
    return d_ins_hier_integrator->getNumberOfCycles();
} // getNumberOfCycles

/////////////////////////////// PROTECTED ////////////////////////////////////

void
IBImplicitHierarchyIntegrator::putToDatabaseSpecialized(Pointer<Database> db)
{
    IBHierarchyIntegrator::putToDatabaseSpecialized(db);
    db->putInteger("IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION",
                   IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION);
    return;
} // putToDatabaseSpecialized

/////////////////////////////// PRIVATE //////////////////////////////////////

void
IBImplicitHierarchyIntegrator::getFromRestart()
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
                                 << " not found in restart file." << std::endl);
    }
    int ver = db->getInteger("IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION");
    if (ver != IB_IMPLICIT_STAGGERED_HIERARCHY_INTEGRATOR_VERSION)
    {
        TBOX_ERROR(d_object_name << ":  Restart file version different than class version." << std::endl);
    }
    return;
} // getFromRestart

void
IBImplicitHierarchyIntegrator::iterateSolution(Vec X_new)
{
    const double current_time = d_current_time;
    const double new_time = d_new_time;
    const double half_time = current_time + 0.5 * (new_time - current_time);
    const int cycle_num = d_cycle_num;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                                   d_ins_hier_integrator->getCurrentContext());
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    // Set the current position data.
    d_ib_implicit_ops->setUpdatedPosition(X_new);

    // Compute the Lagrangian forces and spread them to the Eulerian grid.
    switch (d_time_stepping_type)
    {
    case FORWARD_EULER:
        // intentionally blank
        break;
    case BACKWARD_EULER:
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(d_new_time);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), new_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        break;
    case MIDPOINT_RULE:
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(half_time);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        break;
    case TRAPEZOIDAL_RULE:
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(new_time);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), new_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_hier_velocity_data_ops->linearSum(d_f_idx, 0.5, d_f_current_idx, 0.5, d_f_idx);
        break;
    default:
        TBOX_ERROR(
            d_object_name << "::integrateHierarchy():\n"
                          << "  unsupported time stepping type: "
                          << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                          << "  supported time stepping types are: BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Solve the incompressible Navier-Stokes equations.
    if (d_enable_logging)
        plog << d_object_name << "::integrateHierarchy(): solving the incompressible Navier-Stokes equations\n";
    d_ib_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    d_ins_hier_integrator->integrateHierarchy(current_time, new_time, d_ins_cycle_num);
    d_ins_cycle_num++;
    d_ib_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    switch (d_time_stepping_type)
    {
    case BACKWARD_EULER:
        d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             new_time);
        break;
    case MIDPOINT_RULE:
        d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             half_time);
        break;
    case TRAPEZOIDAL_RULE:
        d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             new_time);
        break;
    default:
        TBOX_ERROR(
            d_object_name << "::integrateHierarchy():\n"
                          << "  unsupported time stepping type: "
                          << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                          << "  supported time stepping types are: BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute the residual.
    switch (d_time_stepping_type)
    {
    case BACKWARD_EULER:
        d_ib_implicit_ops->backwardEulerStep(current_time, new_time);
        break;
    case MIDPOINT_RULE:
        d_ib_implicit_ops->midpointStep(current_time, new_time);
        break;
    case TRAPEZOIDAL_RULE:
        d_ib_implicit_ops->trapezoidalStep(current_time, new_time);
        break;
    default:
        TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                 << "  unsupported time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                                 << "  supported time stepping types are: BACKWARD_EULER, "
                                    "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }
    return;
}

PetscErrorCode
IBImplicitHierarchyIntegrator::IBFunction_SAMRAI(SNES snes, Vec X, Vec R, void* ctx)
{
    auto integrator = static_cast<IBImplicitHierarchyIntegrator*>(ctx);
    return integrator->IBFunction(snes, X, R);
}

PetscErrorCode
IBImplicitHierarchyIntegrator::IBFunction(SNES /*snes*/, Vec X, Vec R)
{
    const double current_time = d_current_time;
    const double new_time = d_new_time;
    const double half_time = current_time + 0.5 * (new_time - current_time);
    const int cycle_num = d_cycle_num;

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int u_current_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                                   d_ins_hier_integrator->getCurrentContext());
    const int u_new_idx = var_db->mapVariableAndContextToIndex(d_ins_hier_integrator->getVelocityVariable(),
                                                               d_ins_hier_integrator->getNewContext());

    // Set the current position data.
    d_ib_implicit_ops->setUpdatedPosition(X);

    // Compute the Lagrangian forces and spread them to the Eulerian grid.
    switch (d_time_stepping_type)
    {
    case FORWARD_EULER:
        // intentionally blank
        break;
    case BACKWARD_EULER:
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(d_new_time);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), new_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        break;
    case MIDPOINT_RULE:
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(half_time);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), half_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        break;
    case TRAPEZOIDAL_RULE:
        if (d_enable_logging) plog << d_object_name << "::integrateHierarchy(): computing Lagrangian force\n";
        d_ib_method_ops->computeLagrangianForce(new_time);
        if (d_enable_logging)
            plog << d_object_name << "::integrateHierarchy(): spreading Lagrangian force to the Eulerian grid\n";
        d_hier_velocity_data_ops->setToScalar(d_f_idx, 0.0);
        d_u_phys_bdry_op->setPatchDataIndex(d_f_idx);
        d_u_phys_bdry_op->setHomogeneousBc(true);
        d_ib_method_ops->spreadForce(
            d_f_idx, d_u_phys_bdry_op, getProlongRefineSchedules(d_object_name + "::f"), new_time);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_hier_velocity_data_ops->linearSum(d_f_idx, 0.5, d_f_current_idx, 0.5, d_f_idx);
        break;
    default:
        TBOX_ERROR(
            d_object_name << "::integrateHierarchy():\n"
                          << "  unsupported time stepping type: "
                          << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                          << "  supported time stepping types are: BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Solve the incompressible Navier-Stokes equations.
    if (d_enable_logging)
        plog << d_object_name << "::integrateHierarchy(): solving the incompressible Navier-Stokes equations\n";
    d_ib_method_ops->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    d_ins_hier_integrator->integrateHierarchy(current_time, new_time, d_ins_cycle_num);
    d_ins_cycle_num++;
    d_ib_method_ops->postprocessSolveFluidEquations(current_time, new_time, cycle_num);

    // Interpolate the Eulerian velocity to the curvilinear mesh.
    switch (d_time_stepping_type)
    {
    case BACKWARD_EULER:
        d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             new_time);
        break;
    case MIDPOINT_RULE:
        d_hier_velocity_data_ops->linearSum(d_u_idx, 0.5, u_current_idx, 0.5, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             half_time);
        break;
    case TRAPEZOIDAL_RULE:
        d_hier_velocity_data_ops->copyData(d_u_idx, u_new_idx);
        if (d_enable_logging)
            plog << d_object_name
                 << "::integrateHierarchy(): interpolating Eulerian velocity to "
                    "the Lagrangian mesh\n";
        d_u_phys_bdry_op->setPatchDataIndex(d_u_idx);
        d_u_phys_bdry_op->setHomogeneousBc(false);
        d_ib_method_ops->interpolateVelocity(d_u_idx,
                                             getCoarsenSchedules(d_object_name + "::u::CONSERVATIVE_COARSEN"),
                                             getGhostfillRefineSchedules(d_object_name + "::u"),
                                             new_time);
        break;
    default:
        TBOX_ERROR(
            d_object_name << "::integrateHierarchy():\n"
                          << "  unsupported time stepping type: "
                          << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                          << "  supported time stepping types are: BACKWARD_EULER, MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }

    // Compute the residual.
    switch (d_time_stepping_type)
    {
    case BACKWARD_EULER:
        d_ib_implicit_ops->computeResidualBackwardEuler(R);
        break;
    case MIDPOINT_RULE:
        d_ib_implicit_ops->computeResidualMidpointRule(R);
        break;
    case TRAPEZOIDAL_RULE:
        d_ib_implicit_ops->computeResidualTrapezoidalRule(R);
        break;
    default:
        TBOX_ERROR(d_object_name << "::integrateHierarchy():\n"
                                 << "  unsupported time stepping type: "
                                 << enum_to_string<TimeSteppingType>(d_time_stepping_type) << "\n"
                                 << "  supported time stepping types are: BACKWARD_EULER, "
                                    "MIDPOINT_RULE, TRAPEZOIDAL_RULE\n");
    }
    return 0;
}

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
