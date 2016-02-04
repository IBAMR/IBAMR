// Filename: main.cpp
// Created on 15 May 2015 by Amneet Bhalla

// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesIBBoxRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesFACPreconditioner.h>
#include <ibamr/StaggeredStokesSolver.h>
#include <ibamr/StaggeredStokesOperator.h>
#include <ibamr/app_namespaces.h>
#include <ibamr/ibamr_utilities.h>
#include <ibtk/CartSideDoubleRT0Coarsen.h>
#include <ibtk/CartSideDoubleSpecializedConstantRefine.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <petscksp.h>

// StokesIBSolver class solves linearized StokesIB equations.
// The class can use two versions of Mat-Vec multiply routine,
// in the Krylov solver (prefix is ib_). One version uses
// matrix-free [matApply()] while the other uses matrix version
// of SAJ operator [matApply2(): by default], to evaluate the IB
// part of the equation on the finest grid level. This can be
// switched in the initializeSolver() routine. The SAJ PETSc
// Mat passed to this class should be for the finest grid level.

class StokesIBSolver : public IBAMR::StaggeredStokesSolver
{
public:
    StokesIBSolver(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                   Pointer<StaggeredStokesIBBoxRelaxationFACOperator> fac_op,
                   Pointer<StaggeredStokesFACPreconditioner> fac_pc,
                   Pointer<IBMethod> ib_method_ops)
    {
        d_hierarchy = patch_hierarchy;
        d_stokes_op = new StaggeredStokesOperator("StokesIBSolver::stokes_op", false);
        d_fac_op = fac_op;
        d_fac_pc = fac_pc;
        d_ib_ops = ib_method_ops;

        // Some default options
        d_options_prefix = "ib_";
        d_ksp_type = "fgmres";
        d_initial_guess_nonzero = true;
        d_rel_residual_tol = 1e-10;
        d_abs_residual_tol = 1e-50;
        d_max_iterations = 1000;

        d_finest_ln = d_hierarchy->getFinestLevelNumber();
        d_hier_velocity_data_ops = new HierarchySideDataOpsReal<NDIM, double>(d_hierarchy, 0, d_finest_ln);
        return;

    } // StokesIBSolver

    void setScratchPatchDataIndices(const int u_idx, const int f_idx)
    {
        d_u_idx = u_idx;
        d_f_idx = f_idx;

        return;

    } // setScratchPatchDataIndices

    void setTimeInterval(double current_time, double new_time)
    {
        StaggeredStokesSolver::setTimeInterval(current_time, new_time);
        d_stokes_op->setTimeInterval(current_time, new_time);
        d_fac_pc->setTimeInterval(current_time, new_time);
        d_fac_op->setTimeInterval(current_time, new_time);

        return;
    } // setTimeInterval

    void setSolutionTime(double solution_time)
    {
        StaggeredStokesSolver::setSolutionTime(solution_time);
        d_stokes_op->setSolutionTime(solution_time);
        d_fac_pc->setSolutionTime(solution_time);
        d_fac_op->setSolutionTime(solution_time);

        return;
    }

    void setVelocityPoissonSpecifications(const PoissonSpecifications& U_problem_coefs)
    {
        StaggeredStokesSolver::setVelocityPoissonSpecifications(U_problem_coefs);
        d_stokes_op->setVelocityPoissonSpecifications(U_problem_coefs);
        d_fac_pc->setVelocityPoissonSpecifications(U_problem_coefs);
        d_fac_op->setVelocityPoissonSpecifications(U_problem_coefs);

        return;
    } // setVelocityPoissonSpecifications

    void setPhysicalBoundaryHelper(Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper)
    {
        StaggeredStokesSolver::setPhysicalBoundaryHelper(bc_helper);
        d_stokes_op->setPhysicalBoundaryHelper(bc_helper);
        d_fac_pc->setPhysicalBoundaryHelper(bc_helper);
        d_fac_op->setPhysicalBoundaryHelper(bc_helper);

        return;
    } // setPhysicalBoundaryHelper

    void setPhysicalBcCoefs(const std::vector<RobinBcCoefStrategy<NDIM>*>& U_bc_coefs,
                            RobinBcCoefStrategy<NDIM>* P_bc_coef)
    {
        StaggeredStokesSolver::setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);
        d_stokes_op->setPhysicalBcCoefs(U_bc_coefs, P_bc_coef);

        return;
    } // setPhysicalBcCoefs

    void setComponentsHaveNullspace(const bool has_velocity_nullspace, const bool has_pressure_nullspace)
    {
        StaggeredStokesSolver::setComponentsHaveNullspace(has_velocity_nullspace, has_pressure_nullspace);
        d_fac_pc->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);
        d_fac_op->setComponentsHaveNullspace(d_has_velocity_nullspace, d_has_pressure_nullspace);

        return;
    } // setComponentsHaveNullspace

    void initializeSolver(Vec x, Vec b)
    {
        d_finest_ln = d_hierarchy->getFinestLevelNumber();
        d_hier_velocity_data_ops->resetLevels(0, d_finest_ln);

        // Allocate space for scratch idx.
        {
            for (int ln = 0; ln <= d_finest_ln; ++ln)
            {
                Pointer<PatchLevel<NDIM> > level = d_hierarchy->getPatchLevel(ln);
                if (!level->checkAllocated(d_u_idx)) level->allocatePatchData(d_u_idx);
                if (!level->checkAllocated(d_f_idx)) level->allocatePatchData(d_f_idx);
            }
        }

        // Create ghost fill schedules
        {
            d_ghost_fill_alg.registerRefine(d_u_idx, d_u_idx, d_u_idx, NULL);
            d_ghost_fill_schd = d_ghost_fill_alg.createSchedule(d_hierarchy->getPatchLevel(d_finest_ln));
        }

        // Create KSP with default options
        {
            KSPCreate(PETSC_COMM_WORLD, &d_petsc_ksp);
            KSPSetType(d_petsc_ksp, d_ksp_type.c_str());

            PetscBool initial_guess_nonzero = d_initial_guess_nonzero ? PETSC_TRUE : PETSC_FALSE;
            KSPSetInitialGuessNonzero(d_petsc_ksp, initial_guess_nonzero);
            KSPSetTolerances(d_petsc_ksp, d_rel_residual_tol, d_abs_residual_tol, PETSC_DEFAULT, d_max_iterations);

            // Get command line options for the KSP.
            KSPSetOptionsPrefix(d_petsc_ksp, d_options_prefix.c_str());
            KSPSetFromOptions(d_petsc_ksp);

            // Reset the member state variables to correspond to the values used by the
            // KSP object.  (Command-line options always take precedence.)
            KSPType ksp_type;
            KSPGetType(d_petsc_ksp, (const char**)&ksp_type);
            d_ksp_type = ksp_type;
            PetscBool nz_init_guess;
            KSPGetInitialGuessNonzero(d_petsc_ksp, &nz_init_guess);
            d_initial_guess_nonzero = (nz_init_guess == PETSC_TRUE);
            KSPGetTolerances(d_petsc_ksp, &d_rel_residual_tol, &d_abs_residual_tol, PETSC_NULL, &d_max_iterations);

            // Set the KSP operator.
            if (d_ksp_mat)
            {
                MatDestroy(&d_ksp_mat);
                d_ksp_mat = PETSC_NULL;
            }
            if (!d_ksp_mat)
            {
                int n;
                VecGetLocalSize(b, &n);
                MatCreateShell(
                    PETSC_COMM_WORLD, n, n, PETSC_DETERMINE, PETSC_DETERMINE, static_cast<void*>(this), &d_ksp_mat);
            }
            MatShellSetOperation(d_ksp_mat, MATOP_MULT, reinterpret_cast<void (*)(void)>(StokesIBSolver::matApply2));

            // Set the shell matrix with KSP.
            if (d_petsc_ksp)
            {
                KSPSetOperators(d_petsc_ksp, d_ksp_mat, d_ksp_mat);
                KSPSetReusePreconditioner(d_petsc_ksp, PETSC_TRUE);
            }

            // Determine the preconditioner type to use.
            static const size_t len = 255;
            char pc_type_str[len];
            PetscBool flg;
#if (!PETSC_VERSION_RELEASE)
            PetscOptionsGetString(NULL, d_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg);
#else
            PetscOptionsGetString(d_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg);
#endif

            std::string pc_type = "shell";
            if (flg)
            {
                pc_type = std::string(pc_type_str);
            }

            PC petsc_pc;
            KSPGetPC(d_petsc_ksp, &petsc_pc);
            if (pc_type == "none")
            {
                PCSetType(petsc_pc, PCNONE);
            }
            else if (pc_type == "shell")
            {
                const std::string pc_name = pc_type;
                PCSetType(petsc_pc, PCSHELL);
                PCShellSetContext(petsc_pc, static_cast<void*>(this));
                PCShellSetApply(petsc_pc, StokesIBSolver::pcApply);
            }
        }

        // Initialize Stokes op and IB FAC pc
        {
            Pointer<SAMRAIVectorReal<NDIM, double> > u_p = PETScSAMRAIVectorReal::getSAMRAIVector(x);
            Pointer<SAMRAIVectorReal<NDIM, double> > f_g = PETScSAMRAIVectorReal::getSAMRAIVector(b);
            d_stokes_op->initializeOperatorState(*u_p, *f_g);
            d_fac_pc->initializeSolverState(*u_p, *f_g);
        }

        return;
    } // initializeSolver

    void registerSAJ(Mat& SAJ, int u_dof_idx, int p_dof_idx)
    {
        d_SAJ = SAJ;
        d_u_dof_idx = u_dof_idx;
        d_p_dof_idx = p_dof_idx;

        return;
    } // registerSAJ

    void deallocateSolver()
    {
        KSPDestroy(&d_petsc_ksp);
        d_stokes_op->deallocateOperatorState();
        d_fac_pc->deallocateSolverState();
        if (d_ksp_mat)
        {
            MatDestroy(&d_ksp_mat);
            d_ksp_mat = PETSC_NULL;
        }

        return;
    } // deallocateSolverState

    bool solveSystem(SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& /*x*/,
                     SAMRAI::solv::SAMRAIVectorReal<NDIM, double>& /*b*/)
    {
        // intentionally left blank.
        return false;
    } // solveSystem

    void solveSystem(Vec x, Vec b)
    {
        Vec petsc_x = x;
        Vec petsc_b;
        VecDuplicate(b, &petsc_b);
        VecCopy(b, petsc_b);
        Pointer<SAMRAIVectorReal<NDIM, double> > f_g = PETScSAMRAIVectorReal::getSAMRAIVector(petsc_b);

        d_stokes_op->setHomogeneousBc(false);
        d_stokes_op->modifyRhsForBcs(*f_g);
        d_stokes_op->setHomogeneousBc(true);

        KSPSolve(d_petsc_ksp, petsc_b, petsc_x);

        Pointer<SAMRAIVectorReal<NDIM, double> > u_p = PETScSAMRAIVectorReal::getSAMRAIVector(petsc_x);
        d_stokes_op->imposeSolBcs(*u_p);

        VecDestroy(&petsc_b);

        return;
    } // solveSystem

private:
    Pointer<PatchHierarchy<NDIM> > d_hierarchy;
    Pointer<StaggeredStokesOperator> d_stokes_op;
    Pointer<StaggeredStokesIBBoxRelaxationFACOperator> d_fac_op;
    Pointer<StaggeredStokesFACPreconditioner> d_fac_pc;
    Pointer<IBMethod> d_ib_ops;

    KSP d_petsc_ksp;
    std::string d_options_prefix, d_ksp_type;
    bool d_initial_guess_nonzero;
    double d_rel_residual_tol, d_abs_residual_tol;
    int d_max_iterations;
    Mat d_ksp_mat;
    Mat d_SAJ;

    double d_current_time, d_new_time;
    int d_u_idx, d_f_idx;
    int d_u_dof_idx, d_p_dof_idx;
    int d_finest_ln;
    RefineAlgorithm<NDIM> d_ghost_fill_alg;
    Pointer<RefineSchedule<NDIM> > d_ghost_fill_schd;
    Pointer<HierarchySideDataOpsReal<NDIM, double> > d_hier_velocity_data_ops;

    static PetscErrorCode matApply(Mat A, Vec x, Vec y)
    {
        void* p_ctx;
        MatShellGetContext(A, &p_ctx);
        StokesIBSolver* solver = static_cast<StokesIBSolver*>(p_ctx);

        Pointer<SAMRAIVectorReal<NDIM, double> > u_p = PETScSAMRAIVectorReal::getSAMRAIVector(x);
        Pointer<SAMRAIVectorReal<NDIM, double> > f_g = PETScSAMRAIVectorReal::getSAMRAIVector(y);

        const int u_idx = u_p->getComponentDescriptorIndex(0);
        const int f_u_idx = f_g->getComponentDescriptorIndex(0);
        const double half_time = 0.5 * (solver->d_current_time + solver->d_new_time);

        // Evaluate the Eulerian terms.
        solver->d_stokes_op->setHomogeneousBc(true);
        solver->d_stokes_op->apply(*u_p, *f_g);

        // Compute position residual X = dt*J*[u/2] = 0 - dt*J*[-u/2]
        Vec X, X0;
        solver->d_ib_ops->createSolverVecs(&X, &X0);
        solver->d_ib_ops->setupSolverVecs(PETSC_NULL, &X0);

        solver->d_hier_velocity_data_ops->scale(solver->d_u_idx, -0.5, u_idx);
        solver->d_ghost_fill_schd->fillData(half_time);
        solver->d_ib_ops->interpolateLinearizedVelocity(solver->d_u_idx,
                                                        std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                                                        std::vector<Pointer<RefineSchedule<NDIM> > >(),
                                                        half_time);
        solver->d_ib_ops->computeLinearizedResidual(X0, X);

        // Compute linearized force F = A/2*[dt*J*[u/2]]
        solver->d_ib_ops->computeLinearizedLagrangianForce(X, half_time);
        solver->d_hier_velocity_data_ops->setToScalar(solver->d_f_idx, 0.0);
        solver->d_ib_ops->spreadLinearizedForce(
            solver->d_f_idx, NULL, std::vector<Pointer<RefineSchedule<NDIM> > >(), half_time);
        solver->d_hier_velocity_data_ops->subtract(f_u_idx, f_u_idx, solver->d_f_idx);
        PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
        PetscFunctionReturn(0);
    } // matApply

    static PetscErrorCode matApply2(Mat A, Vec x, Vec y)
    {
        void* p_ctx;
        MatShellGetContext(A, &p_ctx);
        StokesIBSolver* solver = static_cast<StokesIBSolver*>(p_ctx);
        Pointer<PatchLevel<NDIM> > finest_level = solver->d_hierarchy->getPatchLevel(solver->d_finest_ln);

        Pointer<SAMRAIVectorReal<NDIM, double> > u_p = PETScSAMRAIVectorReal::getSAMRAIVector(x);
        Pointer<SAMRAIVectorReal<NDIM, double> > f_g = PETScSAMRAIVectorReal::getSAMRAIVector(y);
        const int u_idx = u_p->getComponentDescriptorIndex(0);
        const int p_idx = u_p->getComponentDescriptorIndex(1);
        const int f_u_idx = f_g->getComponentDescriptorIndex(0);

        Pointer<SAMRAIVectorReal<NDIM, double> > f_g_duplicate = f_g->cloneVector("");
        f_g_duplicate->allocateVectorData();
        f_g_duplicate->setToScalar(0.0);
        const int f_u_dup_idx = f_g_duplicate->getComponentDescriptorIndex(0);
        const int g_p_dup_idx = f_g_duplicate->getComponentDescriptorIndex(1);

        Vec left, right;
        MatCreateVecs(solver->d_SAJ, &right, &left);
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            right, u_idx, solver->d_u_dof_idx, p_idx, solver->d_p_dof_idx, finest_level);
        MatMult(solver->d_SAJ, right, left);
        Pointer<RefineSchedule<NDIM> > ghost_fill_sched =
            StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(f_u_dup_idx, g_p_dup_idx, finest_level);
        Pointer<RefineSchedule<NDIM> > data_synch_sched =
            StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(f_u_dup_idx, g_p_dup_idx, finest_level);
        StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(left,
                                                                f_u_dup_idx,
                                                                solver->d_u_dof_idx,
                                                                g_p_dup_idx,
                                                                solver->d_p_dof_idx,
                                                                finest_level,
                                                                data_synch_sched,
                                                                ghost_fill_sched);

        // Evaluate the Eulerian terms.
        solver->d_stokes_op->setHomogeneousBc(true);
        solver->d_stokes_op->apply(*u_p, *f_g);
        solver->d_hier_velocity_data_ops->add(f_u_idx, f_u_idx, f_u_dup_idx);
        VecDestroy(&left);
        VecDestroy(&right);
        f_g_duplicate->deallocateVectorData();
        f_g_duplicate->freeVectorComponents();
        PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
        PetscFunctionReturn(0);
    } // matApply2

    static PetscErrorCode pcApply(PC pc, Vec x, Vec y)
    {
        // Here we are solving the equation of the type : Py = x
        // in which P is the preconditioner.
        void* ctx;
        PCShellGetContext(pc, &ctx);
        StokesIBSolver* solver = static_cast<StokesIBSolver*>(ctx);
        Pointer<SAMRAIVectorReal<NDIM, double> > f_g = PETScSAMRAIVectorReal::getSAMRAIVector(x);
        Pointer<SAMRAIVectorReal<NDIM, double> > u_p = PETScSAMRAIVectorReal::getSAMRAIVector(y);
        solver->d_fac_pc->solveSystem(*u_p, *f_g);
        PetscObjectStateIncrease(reinterpret_cast<PetscObject>(y));
        PetscFunctionReturn(0);
    } // pcApply
};

// Function prototype to build coarse level SAJ from basis vectors.
void buildSAJCoarsestFromSAMRAIOperators(Mat& SAJ_coarse,
                                         Mat& SAJ_fine,
                                         std::vector<std::vector<int> > num_dofs_per_proc,
                                         Pointer<SideVariable<NDIM, double> > u_var,
                                         Pointer<CellVariable<NDIM, double> > p_var,
                                         const int u_idx,
                                         const int p_idx,
                                         const int u_dof_index_idx,
                                         const int p_dof_index_idx,
                                         Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                         Pointer<HierarchySideDataOpsReal<NDIM, double> > hier_velocity_data_ops,
                                         Pointer<HierarchyCellDataOpsReal<NDIM, double> > hier_pressure_data_ops,
                                         IntVector<NDIM> gcw);

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));

        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBImplicitStaggeredHierarchyIntegrator("IBHierarchyIntegrator",
                                                       app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                                       ib_method_ops,
                                                       navier_stokes_integrator);

        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

        // Create the IB FAC op/pc and StokesIBSolver.
        Pointer<Database> stokes_ib_precond_db = input_db->getDatabase("stokes_ib_precond_db");
        Pointer<StaggeredStokesIBBoxRelaxationFACOperator> fac_op = new StaggeredStokesIBBoxRelaxationFACOperator(
            "StaggeredStokesIBBoxRelaxationFACOperator", stokes_ib_precond_db, "stokes_ib_pc_");
        Pointer<StaggeredStokesFACPreconditioner> fac_pc =
            new StaggeredStokesFACPreconditioner("StaggeredStokesFACPC", fac_op, stokes_ib_precond_db, "stokes_ib_pc_");
        Pointer<StokesIBSolver> stokes_ib_solver = new StokesIBSolver(patch_hierarchy, fac_op, fac_pc, ib_method_ops);
        navier_stokes_integrator->setStokesSolver(stokes_ib_solver);

        // Create Eulerian initial condition specification objects.  These
        // objects also are used to specify exact solution values for error
        // analysis.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        // navier_stokes_integrator->registerVelocityInitialConditions(u_init);

        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        // navier_stokes_integrator->registerPressureInitialConditions(p_init);

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        LDataManager* lag_data_manager = ib_method_ops->getLDataManager();
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
        }

        //====================================================================
        //    START BUILDING THE DATA STRUCTURES OF INTEGRATORS AND SOLVERS

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
        }

        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        const double current_time = 0.0;
        const double dt = time_integrator->getMaximumTimeStepSize();
        const double new_time = current_time + dt;
        const double half_time = 0.5 * (current_time + new_time);
        const int current_num_cycles = time_integrator->getNumberOfCycles();

        // Initialize IB data and INS data.
        time_integrator->preprocessIntegrateHierarchy(current_time, new_time, current_num_cycles);
        ib_method_ops->updateFixedLEOperators();

        // Compute u and p DOFs per processor.
        std::vector<std::vector<int> > num_dofs_per_proc;
        Pointer<SideVariable<NDIM, int> > u_dof_index_var = new SideVariable<NDIM, int>("u_dof_index");
        ;
        Pointer<CellVariable<NDIM, int> > p_dof_index_var = new CellVariable<NDIM, int>("p_dof_index");
        const IntVector<NDIM> ib_ghosts = ib_method_ops->getMinimumGhostCellWidth();
        const IntVector<NDIM> no_ghosts = 0;
        const int u_dof_index_idx =
            var_db->registerVariableAndContext(u_dof_index_var, time_integrator->getScratchContext(), ib_ghosts);
        const int p_dof_index_idx =
            var_db->registerVariableAndContext(p_dof_index_var, time_integrator->getScratchContext(), no_ghosts);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        num_dofs_per_proc.resize(finest_ln + 1);
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);

            level->allocatePatchData(u_dof_index_idx, current_time);
            level->allocatePatchData(p_dof_index_idx, current_time);
            StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
                num_dofs_per_proc[ln], u_dof_index_idx, p_dof_index_idx, level);
        }

        //====================================================================
        // Get a sense of Lagrangian nodes distribution among processors
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,
                                "Local Lagrangian nodes on proc [%d] are : %d\n",
                                SAMRAI_MPI::getRank(),
                                lag_data_manager->getNumberOfLocalNodes(finest_ln));
        PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

        // Get the matrix/matrix-free representation of force Jacobian (A).
        Mat A = PETSC_NULL, A_MFFD = PETSC_NULL;
        ib_method_ops->constructLagrangianForceJacobian(A, MATAIJ);
        ib_method_ops->constructLagrangianForceJacobian(A_MFFD, MATSHELL);

        // ===================================================================
        // TEST matrix (or matrix-free) version of elasticity op.
        Vec dX, AxdX;
        ib_method_ops->createSolverVecs(&dX, &AxdX);
        VecSet(dX, 1.0);
        VecSet(AxdX, 0.0);
        MatMult(A, dX, AxdX);
        VecView(AxdX, PETSC_VIEWER_STDOUT_WORLD);
        //====================================================================

        // Get the matrix representation of J at the finest level
        Mat J = PETSC_NULL;
        ib_method_ops->constructInterpOp(J,
                                         PETScMatUtilities::ib_4_interp_fcn,
                                         PETScMatUtilities::ib_4_interp_stencil,
                                         num_dofs_per_proc[finest_ln],
                                         u_dof_index_idx);

        // Configure the fac pc/op
        fac_pc->setPhysicalBcCoefs(navier_stokes_integrator->getIntermediateVelocityBoundaryConditions(),
                                   navier_stokes_integrator->getProjectionBoundaryConditions());
        fac_op->setPhysicalBcCoefs(navier_stokes_integrator->getIntermediateVelocityBoundaryConditions(),
                                   navier_stokes_integrator->getProjectionBoundaryConditions());
        fac_op->setIBForceJacobian(A);
        fac_op->setIBInterpOp(J);

        // ===================  SAJ at the finest level ======================

        // Note that FAC op builds SAJ at the finest level in a similar way. We
        // are creating it here, as we want to pass it to StokesIBSolver class.

        Mat SAJ;
        MatPtAP(A, J, MAT_INITIAL_MATRIX, 1.0, &SAJ);

        // Compute the scale for the spreading operator.
        Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);
        Pointer<CartesianGridGeometry<NDIM> > grid_geom = patch_hierarchy->getGridGeometry();
        const double* const dx0 = grid_geom->getDx();
        IntVector<NDIM> ratio = finest_level->getRatio();
        double spread_scale = -0.25 * (dt);
        for (unsigned d = 0; d < NDIM; ++d) spread_scale *= ratio(d) / dx0[d];
        MatScale(SAJ, spread_scale);
        // ===================================================================

        // Create variables for velocity, force, pressure and incompressibility.
        Pointer<VariableContext> ib_ctx = var_db->getContext("ib_ctx");
        Pointer<VariableContext> ins_ctx = var_db->getContext("ins_ctx");
        Pointer<SideVariable<NDIM, double> > u_var = new SideVariable<NDIM, double>("u_var");
        Pointer<SideVariable<NDIM, double> > f_var = new SideVariable<NDIM, double>("f_var");
        Pointer<CellVariable<NDIM, double> > p_var = new CellVariable<NDIM, double>("p_var");
        Pointer<CellVariable<NDIM, double> > g_var = new CellVariable<NDIM, double>("g_var");
        const int u_ib_idx = var_db->registerVariableAndContext(u_var, ib_ctx, lag_data_manager->getGhostCellWidth());
        const int f_ib_idx = var_db->registerVariableAndContext(f_var, ib_ctx, lag_data_manager->getGhostCellWidth());
        const int u_ins_idx = var_db->registerVariableAndContext(u_var, ins_ctx, 1);
        const int f_ins_idx = var_db->registerVariableAndContext(f_var, ins_ctx, 1);
        const int p_ins_idx = var_db->registerVariableAndContext(p_var, ins_ctx, 1);
        const int g_ins_idx = var_db->registerVariableAndContext(g_var, ins_ctx, 1);

        Pointer<HierarchySideDataOpsReal<NDIM, double> > hier_velocity_data_ops =
            new HierarchySideDataOpsReal<NDIM, double>(patch_hierarchy, coarsest_ln, finest_ln);
        Pointer<HierarchyCellDataOpsReal<NDIM, double> > hier_pressure_data_ops =
            new HierarchyCellDataOpsReal<NDIM, double>(patch_hierarchy, coarsest_ln, finest_ln);

        Pointer<HierarchyMathOps> hier_math_ops = time_integrator->getHierarchyMathOps();
        const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops->getSideWeightPatchDescriptorIndex();

        // Setup Eulerian vectors used in solving the linear implicit IB equations.
        Pointer<SAMRAIVectorReal<NDIM, double> > eul_sol_vec =
            new SAMRAIVectorReal<NDIM, double>("eul_sol_vec", patch_hierarchy, coarsest_ln, finest_ln);
        eul_sol_vec->addComponent(u_var, u_ins_idx, wgt_sc_idx, hier_velocity_data_ops);
        eul_sol_vec->addComponent(p_var, p_ins_idx, wgt_cc_idx, hier_pressure_data_ops);
        eul_sol_vec->allocateVectorData();

        Pointer<SAMRAIVectorReal<NDIM, double> > eul_rhs_vec =
            new SAMRAIVectorReal<NDIM, double>("eul_rhs_vec", patch_hierarchy, coarsest_ln, finest_ln);
        eul_rhs_vec->addComponent(f_var, f_ins_idx, wgt_sc_idx, hier_velocity_data_ops);
        eul_rhs_vec->addComponent(g_var, g_ins_idx, wgt_cc_idx, hier_pressure_data_ops);
        eul_rhs_vec->allocateVectorData();
        //====================================================================

        // Compute convective and previous time-step diffusive terms in the rhs
        // vec and set an initial guess for the solution vec.
        for (int cycle_num = 0; cycle_num < current_num_cycles; ++cycle_num)
        {
            time_integrator->skipCycle(current_time, new_time, cycle_num);
            navier_stokes_integrator->skipCycle(current_time, new_time, cycle_num);
        }
        navier_stokes_integrator->setupSolverVectors(eul_sol_vec, eul_rhs_vec, current_time, new_time, 0);

        //====================================================================
        // Configure the StokesIBSolver and solve the linearized Stokes-IB eqns.
        Vec eul_sol_petsc_vec = PETScSAMRAIVectorReal::createPETScVector(eul_sol_vec, PETSC_COMM_WORLD);
        Vec eul_rhs_petsc_vec = PETScSAMRAIVectorReal::createPETScVector(eul_rhs_vec, PETSC_COMM_WORLD);

        stokes_ib_solver->setScratchPatchDataIndices(u_ib_idx, f_ib_idx);
        stokes_ib_solver->setTimeInterval(current_time, new_time);
        stokes_ib_solver->setSolutionTime(new_time);
        stokes_ib_solver->registerSAJ(SAJ, u_dof_index_idx, p_dof_index_idx);
        stokes_ib_solver->setComponentsHaveNullspace(false, true);
        stokes_ib_solver->initializeSolver(eul_sol_petsc_vec, eul_rhs_petsc_vec);
        stokes_ib_solver->solveSystem(eul_sol_petsc_vec, eul_rhs_petsc_vec);

        // Reset Eulerian solver vectors and Eulerian state data.
        navier_stokes_integrator->resetSolverVectors(eul_sol_vec, eul_rhs_vec, current_time, new_time, 0);

        // Interpolate the Eulerian velocity to the curvilinear mesh.
        RefineAlgorithm<NDIM> ghost_fill_alg;
        ghost_fill_alg.registerRefine(u_ib_idx, u_ib_idx, u_ib_idx, NULL);
        Pointer<RefineSchedule<NDIM> > ghost_fill_schd =
            ghost_fill_alg.createSchedule(patch_hierarchy->getPatchLevel(finest_ln));

        Pointer<SideVariable<NDIM, double> > vel_var = navier_stokes_integrator->getVelocityVariable();
        Pointer<VariableContext> current_ctx = navier_stokes_integrator->getCurrentContext();
        const int u_current_idx = var_db->mapVariableAndContextToIndex(vel_var, current_ctx);

        hier_velocity_data_ops->linearSum(u_ib_idx, 0.5, u_current_idx, 0.5, u_ins_idx);
        ghost_fill_schd->fillData(new_time);
        ib_method_ops->interpolateVelocity(u_ib_idx,
                                           std::vector<Pointer<CoarsenSchedule<NDIM> > >(),
                                           std::vector<Pointer<RefineSchedule<NDIM> > >(),
                                           half_time);

        // Compute the final value of the updated positions of the Lagrangian
        // structure.
        ib_method_ops->midpointStep(current_time, new_time);

        //====================================================================
        //========================= TEST SAJ Operator ========================

        // Get SAJ mat for the coarsest level
        Mat SAJ_coarsest_petsc = fac_op->getEulerianElasticityLevelOp(coarsest_ln);

        // Build SAJ mat for the coarsest level using SAMRAI operators
        Mat SAJ_coarsest_samrai;
        MatDuplicate(SAJ_coarsest_petsc, MAT_SHARE_NONZERO_PATTERN, &SAJ_coarsest_samrai);
        buildSAJCoarsestFromSAMRAIOperators(SAJ_coarsest_samrai,
                                            SAJ,
                                            num_dofs_per_proc,
                                            u_var,
                                            p_var,
                                            u_ib_idx,
                                            p_ins_idx,
                                            u_dof_index_idx,
                                            p_dof_index_idx,
                                            patch_hierarchy,
                                            hier_velocity_data_ops,
                                            hier_pressure_data_ops,
                                            lag_data_manager->getGhostCellWidth());

        // Print both versions of the matrices
        PetscViewer matlab_viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "PETSC_SAJ.dat", FILE_MODE_WRITE, &matlab_viewer);
        PetscViewerSetFormat(matlab_viewer, PETSC_VIEWER_NATIVE);
        MatView(SAJ_coarsest_petsc, matlab_viewer);
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, "SAMRAI_SAJ.dat", FILE_MODE_WRITE, &matlab_viewer);
        PetscViewerSetFormat(matlab_viewer, PETSC_VIEWER_NATIVE);
        MatView(SAJ_coarsest_samrai, matlab_viewer);
        //====================================================================

        // post process hierarchy
        time_integrator->postprocessIntegrateHierarchy(current_time, new_time, true, current_num_cycles);
        time_integrator->synchronizeHierarchyData(NEW_DATA);
        time_integrator->resetTimeDependentHierarchyData(new_time);

        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, 1, new_time);
            silo_data_writer->writePlotData(1, new_time);
        }

        // Deallocate solver and vector components
        stokes_ib_solver->deallocateSolver();
        eul_sol_vec->deallocateVectorData();
        eul_rhs_vec->deallocateVectorData();

        eul_sol_vec->freeVectorComponents();
        eul_rhs_vec->freeVectorComponents();

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

        // Cleanup PETSc objects
        MatDestroy(&SAJ);
        MatDestroy(&SAJ_coarsest_samrai);
        PetscViewerDestroy(&matlab_viewer);

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
} // main

void
buildSAJCoarsestFromSAMRAIOperators(Mat& SAJ_coarse,
                                    Mat& SAJ_fine,
                                    std::vector<std::vector<int> > num_dofs_per_proc,
                                    Pointer<SideVariable<NDIM, double> > u_var,
                                    Pointer<CellVariable<NDIM, double> > /*p_var*/,
                                    const int u_idx,
                                    const int p_idx,
                                    const int u_dof_index_idx,
                                    const int p_dof_index_idx,
                                    Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                    Pointer<HierarchySideDataOpsReal<NDIM, double> > hier_velocity_data_ops,
                                    Pointer<HierarchyCellDataOpsReal<NDIM, double> > hier_pressure_data_ops,
                                    IntVector<NDIM> gcw)
{
    // Level info.
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM> > coarsest_level = patch_hierarchy->getPatchLevel(coarsest_ln);
    Pointer<PatchLevel<NDIM> > finest_level = patch_hierarchy->getPatchLevel(finest_ln);

    // Get the transfer operators.
    Pointer<CartesianGridGeometry<NDIM> > geometry = patch_hierarchy->getGridGeometry();
    IBAMR_DO_ONCE(geometry->addSpatialCoarsenOperator(new CartSideDoubleRT0Coarsen(gcw));
                  geometry->addSpatialRefineOperator(new CartSideDoubleSpecializedConstantRefine()));

    Pointer<RefineOperator<NDIM> > prolongation_op =
        geometry->lookupRefineOperator(u_var, "SPECIALIZED_CONSTANT_REFINE");
    Pointer<CoarsenOperator<NDIM> > restriction_op = geometry->lookupCoarsenOperator(u_var, "RT0_COARSEN");

    // Define the prolongation and refine algorithms
    Pointer<RefineAlgorithm<NDIM> > prolongation_refine_algorithm = new RefineAlgorithm<NDIM>();
    Pointer<CoarsenAlgorithm<NDIM> > restriction_coarsen_algorithm = new CoarsenAlgorithm<NDIM>();
    prolongation_refine_algorithm->registerRefine(u_idx, u_idx, u_idx, prolongation_op, NULL);
    restriction_coarsen_algorithm->registerCoarsen(u_idx, u_idx, restriction_op, NULL);
    Pointer<RefineSchedule<NDIM> > prolongation_schedule = prolongation_refine_algorithm->createSchedule(
        finest_level, Pointer<PatchLevel<NDIM> >(), coarsest_ln, patch_hierarchy, NULL);
    Pointer<CoarsenSchedule<NDIM> > restriction_schedule =
        restriction_coarsen_algorithm->createSchedule(coarsest_level, finest_level);

    // Get DOFs info at the coarse and fine levels.
    const int mpi_rank = SAMRAI_MPI::getRank();
    const int n_local_coarsest = num_dofs_per_proc[coarsest_ln][mpi_rank];
    const int i_lower_coarsest =
        std::accumulate(num_dofs_per_proc[coarsest_ln].begin(), num_dofs_per_proc[coarsest_ln].begin() + mpi_rank, 0);
    const int n_total_coarsest =
        std::accumulate(num_dofs_per_proc[coarsest_ln].begin(), num_dofs_per_proc[coarsest_ln].end(), 0);

    const int n_local_finest = num_dofs_per_proc[finest_ln][mpi_rank];
    const int n_total_finest =
        std::accumulate(num_dofs_per_proc[finest_ln].begin(), num_dofs_per_proc[finest_ln].end(), 0);

    // Construct the coarse and fine level Vecs.
    Vec x, y;
    VecCreateMPI(PETSC_COMM_WORLD, n_local_coarsest, n_total_coarsest, &x);
    VecCreateMPI(PETSC_COMM_WORLD, n_local_coarsest, n_total_coarsest, &y);
    Vec X, Y;
    VecCreateMPI(PETSC_COMM_WORLD, n_local_finest, n_total_finest, &X);
    VecCreateMPI(PETSC_COMM_WORLD, n_local_finest, n_total_finest, &Y);

    // Utility schedules for copying to and from PETSc and SAMRAI vectors.
    Pointer<RefineSchedule<NDIM> > ghost_fill_sched_coarse =
        StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(u_idx, p_idx, coarsest_level);
    Pointer<RefineSchedule<NDIM> > ghost_fill_sched_fine =
        StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(u_idx, p_idx, finest_level);
    Pointer<RefineSchedule<NDIM> > data_synch_sched_coarse =
        StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(u_idx, p_idx, coarsest_level);
    Pointer<RefineSchedule<NDIM> > data_synch_sched_fine =
        StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(u_idx, p_idx, finest_level);

    // Construct the basis vecs and do matrix-free operations on them to
    // form an explicit representation of the matrix (column-by-column).
    for (int col = 0; col < n_total_coarsest; ++col)
    {
        pout << "+++++++++++++++ Filling up column " << col << " ++++++++++++++ \n\n" << std::endl;

        VecZeroEntries(x);
        VecSetValue(x, col, 1.0, INSERT_VALUES);
        VecAssemblyBegin(x);
        VecAssemblyEnd(x);

        StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(x,
                                                                u_idx,
                                                                u_dof_index_idx,
                                                                p_idx,
                                                                p_dof_index_idx,
                                                                coarsest_level,
                                                                data_synch_sched_coarse,
                                                                ghost_fill_sched_coarse);

        // Zero-out fine data
        hier_pressure_data_ops->resetLevels(finest_ln, finest_ln);
        hier_pressure_data_ops->setToScalar(p_idx, 0.0, /*interior_only*/ false);
        hier_pressure_data_ops->resetLevels(coarsest_ln, finest_ln);
        hier_velocity_data_ops->resetLevels(finest_ln, finest_ln);
        hier_velocity_data_ops->setToScalar(u_idx, 0.0, /*interior_only*/ false);
        hier_velocity_data_ops->resetLevels(coarsest_ln, finest_ln);

        prolongation_schedule->fillData(0.0);

        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            X, u_idx, u_dof_index_idx, p_idx, p_dof_index_idx, finest_level);
        MatMult(SAJ_fine, X, Y);

        // Zero-out fine data
        hier_pressure_data_ops->resetLevels(finest_ln, finest_ln);
        hier_pressure_data_ops->setToScalar(p_idx, 0.0, /*interior_only*/ false);
        hier_pressure_data_ops->resetLevels(coarsest_ln, finest_ln);
        hier_velocity_data_ops->resetLevels(finest_ln, finest_ln);
        hier_velocity_data_ops->setToScalar(u_idx, 0.0, /*interior_only*/ false);
        hier_velocity_data_ops->resetLevels(coarsest_ln, finest_ln);

        StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(Y,
                                                                u_idx,
                                                                u_dof_index_idx,
                                                                p_idx,
                                                                p_dof_index_idx,
                                                                finest_level,
                                                                data_synch_sched_fine,
                                                                ghost_fill_sched_fine);
        // Zero-out coarse data
        hier_pressure_data_ops->resetLevels(coarsest_ln, coarsest_ln);
        hier_pressure_data_ops->setToScalar(p_idx, 0.0);
        hier_pressure_data_ops->resetLevels(coarsest_ln, finest_ln);
        hier_velocity_data_ops->resetLevels(coarsest_ln, coarsest_ln);
        hier_velocity_data_ops->setToScalar(u_idx, 0.0);
        hier_velocity_data_ops->resetLevels(coarsest_ln, finest_ln);

        restriction_schedule->coarsenData();
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            y, u_idx, u_dof_index_idx, p_idx, p_dof_index_idx, coarsest_level);

        PetscScalar* y_array;
        VecGetArray(y, &y_array);
        std::vector<double> nnz_values;
        std::vector<int> nnz_indices;
        nnz_values.reserve(n_local_coarsest);
        nnz_indices.reserve(n_local_coarsest);
        for (int j = 0; j < n_local_coarsest; ++j)
        {
            if (!MathUtilities<double>::equalEps(y_array[j], 0.0))
            {
                const int global_idx = i_lower_coarsest + j;
                nnz_indices.push_back(global_idx);
                nnz_values.push_back(y_array[j]);
            }
        }
        int nnz_size = static_cast<int>(nnz_indices.size());
        MatSetValues(SAJ_coarse, nnz_size, &nnz_indices[0], 1, &col, &nnz_values[0], INSERT_VALUES);
        VecRestoreArray(y, &y_array);
    }

    MatAssemblyBegin(SAJ_coarse, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(SAJ_coarse, MAT_FINAL_ASSEMBLY);

    return;
} // buildSAJCoarsestFromSAMRAIOperators
