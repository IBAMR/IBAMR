// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

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
#include <VariableDatabase.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/CIBMethod.h>
#include <ibamr/CIBMobilitySolver.h>
#include <ibamr/CIBSaddlePointSolver.h>
#include <ibamr/CIBStaggeredStokesSolver.h>
#include <ibamr/DirectMobilitySolver.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/IBStrategySet.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/KrylovMobilitySolver.h>
#include <ibamr/app_namespaces.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Headers for complex fluids
#include "ibamr/CFINSForcing.h"
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>

//////////////////////////////////////////////////////////////////////////////

// Center of mass velocity
void
ConstrainedUpperVel(double /*data_time*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    return;
} // ConstrainedCOMVel

void
ConstrainedLowerVel(double /*data_time*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    return;
} // ConstrainedCOMVel
void
ConstrainedLeftVel(double /*data_time*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    return;
} // ConstrainedCOMVel

void
ConstrainedRightVel(double /*data_time*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    return;
} // ConstrainedCOMVel

void
ConstrainedNodalVel(Vec /*U_k*/, const RigidDOFVector& /*U*/, const Eigen::Vector3d& /*X_com*/, void* /*ctx*/)
{
    // intentionally left blank.
    return;
} // ConstrainedNodalVel

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 LDataManager* l_data_manager,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

struct MaskData
{
public:
    MaskData(const int cc_mask_idx,
             const int sc_mask_idx,
             Pointer<CartGridFunction> cc_mask_fcn,
             Pointer<CartGridFunction> sc_mask_fcn)
        : d_cc_mask_idx(cc_mask_idx), d_sc_mask_idx(sc_mask_idx), d_cc_mask_fcn(cc_mask_fcn), d_sc_mask_fcn(sc_mask_fcn)
    {
    }

    void fillMaskData(double time, Pointer<PatchHierarchy<NDIM> > hierarchy);

private:
    int d_cc_mask_idx = IBTK::invalid_index, d_sc_mask_idx = IBTK::invalid_index;
    Pointer<CartGridFunction> d_cc_mask_fcn, d_sc_mask_fcn;
};

static void
fill_mask_data(Pointer<BasePatchHierarchy<NDIM> > hierarchy, double data_time, bool /*initial_time*/, void* ctx)
{
    auto maskData = static_cast<MaskData*>(ctx);
    maskData->fillMaskData(data_time, hierarchy);
}

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
    SAMRAIManager::setMaxNumberPatchDataEntries(2054);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "CIB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Read default Petsc options
        if (input_db->keyExists("petsc_options_file"))
        {
            std::string petsc_options_file = input_db->getString("petsc_options_file");
            PetscOptionsInsertFile(PETSC_COMM_WORLD, NULL, petsc_options_file.c_str(), PETSC_TRUE);
        }

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.

        // INS integrator
        Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        // CIB method
        const unsigned int num_structures = input_db->getIntegerWithDefault("num_structures", 2);
        Pointer<CIBMethod> cib_method_ops =
            new CIBMethod("CIBMethod", app_initializer->getComponentDatabase("CIBMethod"), num_structures);
        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        std::vector<Pointer<IBStrategy> > ib_strategies{ cib_method_ops, ib_method_ops };
        Pointer<IBStrategy> ib_strategy_set = new IBStrategySet(ib_strategies.begin(), ib_strategies.end());

        // Krylov solver for INS integrator that solves for [u,p,U,L]
        Pointer<CIBStaggeredStokesSolver> CIBSolver =
            new CIBStaggeredStokesSolver("CIBStaggeredStokesSolver",
                                         input_db->getDatabase("CIBStaggeredStokesSolver"),
                                         navier_stokes_integrator,
                                         cib_method_ops,
                                         "SP_");

        // Register the Krylov solver with INS integrator
        navier_stokes_integrator->setStokesSolver(CIBSolver);

        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_strategy_set,
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
        // Configure the CIB solver.
        Pointer<IBStandardInitializer> cib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("CIBStandardInitializer"));
        cib_method_ops->registerLInitStrategy(cib_initializer);

        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

        /*
                // Adv diff
                Pointer<AdvDiffSemiImplicitHierarchyIntegrator> adv_diff_integrator;
                adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
                    "AdvDiffSemiImplicitHierarchyIntegrator",
                    app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
                navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);
        */
        // Specify structure kinematics
        FreeRigidDOFVector upper_dofs, lower_dofs, left_dofs, right_dofs;
        upper_dofs.setZero();
        lower_dofs.setZero();
        left_dofs.setZero();
        right_dofs.setZero();
        cib_method_ops->setSolveRigidBodyVelocity(0, upper_dofs);
        cib_method_ops->registerConstrainedVelocityFunction(NULL, &ConstrainedUpperVel, NULL, 0);
        cib_method_ops->setSolveRigidBodyVelocity(1, lower_dofs);
        cib_method_ops->registerConstrainedVelocityFunction(NULL, &ConstrainedLowerVel, NULL, 1);
        cib_method_ops->setSolveRigidBodyVelocity(2, left_dofs);
        cib_method_ops->registerConstrainedVelocityFunction(NULL, &ConstrainedLeftVel, NULL, 2);
        cib_method_ops->setSolveRigidBodyVelocity(3, right_dofs);
        cib_method_ops->registerConstrainedVelocityFunction(NULL, &ConstrainedRightVel, NULL, 3);

        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

        Pointer<CartGridFunction> sc_mask = new muParserCartGridFunction(
            "sc_mask", app_initializer->getComponentDatabase("VelocityMask"), grid_geometry);
        Pointer<CartGridFunction> cc_mask = new muParserCartGridFunction(
            "cc_mask", app_initializer->getComponentDatabase("PressureMask"), grid_geometry);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> cib_silo_data_writer = app_initializer->getLSiloDataWriter();
        const std::string ib_dirname = app_initializer->getVizDumpDirectory() + "/boundary_ib";
        Utilities::recursiveMkdir(ib_dirname);
        Pointer<LSiloDataWriter> ib_silo_data_writer = new LSiloDataWriter("IBLSiloDataWriter", ib_dirname);

        if (uses_visit)
        {
            cib_initializer->registerLSiloDataWriter(cib_silo_data_writer);
            ib_initializer->registerLSiloDataWriter(ib_silo_data_writer);

            cib_method_ops->registerLSiloDataWriter(cib_silo_data_writer);
            ib_method_ops->registerLSiloDataWriter(ib_silo_data_writer);

            cib_method_ops->registerVisItDataWriter(visit_data_writer);

            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Create boundary condition specification objects (when necessary).
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
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);

                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);

                Pointer<Database> bc_coefs_db = app_initializer->getComponentDatabase(bc_coefs_db_name);
                u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, bc_coefs_db, grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }
        /*
                // Complex fluid forcing
                Pointer<CFINSForcing> complex_fluid;
                bool using_exact_u = input_db->getBool("USING_EXACT_U");
                if (input_db->keyExists("ComplexFluid"))
                {
                    complex_fluid = new CFINSForcing("ComplexFluidForcing",
                        app_initializer->getComponentDatabase("ComplexFluid"),
                        (Pointer<INSHierarchyIntegrator>)navier_stokes_integrator,
                        grid_geometry,
                        adv_diff_integrator,
                        visit_data_writer);
                    time_integrator->registerBodyForceFunction(complex_fluid);
                }
        */
        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        Pointer<muParserCartGridFunction> s_init = new muParserCartGridFunction(
            "S_exact",
            app_initializer->getComponentDatabase("ComplexFluid")->getDatabase("InitialConditions"),
            grid_geometry);

        // Set physical boundary operator used in spreading.
        cib_method_ops->setVelocityPhysBdryOp(time_integrator->getVelocityPhysBdryOp());

        // Register mobility matrices (if needed)
        std::string mobility_solver_type = input_db->getString("MOBILITY_SOLVER_TYPE");
        if (mobility_solver_type == "DIRECT")
        {
            std::string mat_name1 = "MOB-MATRIX-1";
            std::string mat_name2 = "MOB-MATRIX-2";
            std::string mat_name3 = "MOB-MATRIX-3";
            std::string mat_name4 = "MOB-MATRIX-4";
            int managing_proc = 0;
            std::vector<std::vector<unsigned> > struct_ids1;
            std::vector<std::vector<unsigned> > struct_ids2;
            std::vector<std::vector<unsigned> > struct_ids3;
            std::vector<std::vector<unsigned> > struct_ids4;
            std::vector<unsigned> prototype_structs1;
            std::vector<unsigned> prototype_structs2;
            std::vector<unsigned> prototype_structs3;
            std::vector<unsigned> prototype_structs4;
            prototype_structs1.push_back(0);
            struct_ids1.push_back(prototype_structs1);
            prototype_structs2.push_back(1);
            struct_ids2.push_back(prototype_structs2);
            prototype_structs3.push_back(2);
            struct_ids1.push_back(prototype_structs3);
            prototype_structs4.push_back(3);
            struct_ids2.push_back(prototype_structs4);

            DirectMobilitySolver* direct_solvers = NULL;
            CIBSolver->getSaddlePointSolver()->getCIBMobilitySolver()->getMobilitySolvers(NULL, &direct_solvers, NULL);

            direct_solvers->registerMobilityMat(
                mat_name1, prototype_structs1, RPY, std::make_pair(LAPACK_LU, LAPACK_LU), managing_proc);
            direct_solvers->registerStructIDsWithMobilityMat(mat_name1, struct_ids1);

            if (SAMRAI_MPI::getNodes() > 1)
            {
                managing_proc += 1;
            }
            direct_solvers->registerMobilityMat(
                mat_name2, prototype_structs2, RPY, std::make_pair(LAPACK_LU, LAPACK_LU), 0); // managing_proc);
            direct_solvers->registerStructIDsWithMobilityMat(mat_name2, struct_ids2);
            direct_solvers->registerMobilityMat(
                mat_name3, prototype_structs3, RPY, std::make_pair(LAPACK_LU, LAPACK_LU), 0); // managing_proc);
            direct_solvers->registerStructIDsWithMobilityMat(mat_name3, struct_ids3);
            direct_solvers->registerMobilityMat(
                mat_name4, prototype_structs4, RPY, std::make_pair(LAPACK_LU, LAPACK_LU), 0); // managing_proc);
            direct_solvers->registerStructIDsWithMobilityMat(mat_name4, struct_ids4);
        }
        navier_stokes_integrator->setStokesSolverNeedsInit();

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Setup data used to determine the accuracy of the computed solution.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> main_ctx = var_db->getContext("main_context");

        const std::string weighting_fcn = input_db->getString("DELTA_FUNCTION");
        const int n_ghosts = LEInteractor::getMinimumGhostWidth(weighting_fcn);

        const Pointer<Variable<NDIM> > u_var = navier_stokes_integrator->getVelocityVariable();
        const int sc_masked_idx = var_db->registerVariableAndContext(u_var, main_ctx, IntVector<NDIM>(n_ghosts + 2));

        const Pointer<Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
        const int cc_masked_idx = var_db->registerVariableAndContext(p_var, main_ctx, IntVector<NDIM>(n_ghosts + 2));
        visit_data_writer->registerPlotQuantity("cc_mask", "SCALAR", cc_masked_idx);

        MaskData maskData(cc_masked_idx, sc_masked_idx, cc_mask, sc_mask);
        time_integrator->registerRegridHierarchyCallback(&fill_mask_data, static_cast<void*>(&maskData));

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();

        maskData.fillMaskData(loop_time, patch_hierarchy);

        // Register sc mask patch data index with CIBMethod
        cib_method_ops->registerMaskingPatchDataIndex(sc_masked_idx);

        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            cib_silo_data_writer->writePlotData(iteration_num, loop_time);
            ib_silo_data_writer->writePlotData(iteration_num, loop_time);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();

            pout << "Advancing hierarchy by timestep size dt = " << dt << "\n";

            if (time_integrator->atRegridPoint()) navier_stokes_integrator->setStokesSolverNeedsInit();
            if (cib_method_ops->flagRegrid())
            {
                time_integrator->regridHierarchy();
                navier_stokes_integrator->setStokesSolverNeedsInit();
            }
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            RigidDOFVector U0;
            cib_method_ops->getNewRigidBodyVelocity(0, U0);
            pout << "\nRigid body velocity of the structure is : \n" << U0 << "\n";

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                cib_silo_data_writer->writePlotData(iteration_num, loop_time);
                ib_silo_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

        // Cleanup boundary condition specification objects (when necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
} // main

void
MaskData::fillMaskData(const double time, Pointer<PatchHierarchy<NDIM> > hierarchy)
{
    plog << "Filling mask data at time " << time << ".\n";
    const int coarsest_ln = 0;
    const int finest_ln = hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(d_sc_mask_idx)) level->allocatePatchData(d_sc_mask_idx);
        if (!level->checkAllocated(d_cc_mask_idx)) level->allocatePatchData(d_cc_mask_idx);
    }
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<Variable<NDIM> > sc_var;
    Pointer<Variable<NDIM> > cc_var;
    var_db->mapIndexToVariable(d_sc_mask_idx, sc_var);
    var_db->mapIndexToVariable(d_cc_mask_idx, cc_var);
    d_sc_mask_fcn->setDataOnPatchHierarchy(d_sc_mask_idx, sc_var, hierarchy, time);
    d_cc_mask_fcn->setDataOnPatchHierarchy(d_cc_mask_idx, cc_var, hierarchy, time);

    // Fill ghost data.
    using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<ITC> ghost_cell_components(2);
    ghost_cell_components[0] =
        ITC(d_sc_mask_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, {}, nullptr);
    ghost_cell_components[1] =
        ITC(d_cc_mask_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false, {}, nullptr);
    HierarchyGhostCellInterpolation ghost_fill_op;
    ghost_fill_op.initializeOperatorState(ghost_cell_components, hierarchy);
    ghost_fill_op.fillData(time);
}
