// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2026 by the IBAMR developers
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

#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <LoadBalancer.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/CIBMethod.h>
#include <ibamr/CIBMobilitySolver.h>
#include <ibamr/CIBSaddlePointSolver.h>
#include <ibamr/CIBStaggeredStokesOperator.h>
#include <ibamr/CIBStaggeredStokesSolver.h>
#include <ibamr/DirectMobilitySolver.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/KrylovMobilitySolver.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LMesh.h>
#include <ibtk/LNode.h>
#include <ibtk/PETScSAMRAIVectorReal.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

#include <cmath>
#include <iomanip>

#include <ibamr/app_namespaces.h>

//////////////////////////////////////////////////////////////////////////////

// Center of mass velocity
void
ConstrainedCOMVel(double /*data_time*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();
    U_com[0] = 1.0;

    return;
} // ConstrainedCOMOuterVel

void
NetExternalForceTorque(double /*data_time*/, Eigen::Vector3d& F_ext, Eigen::Vector3d& T_ext, void* /*ctx*/)
{
    F_ext << 0.0, 0.0, 0.0;
    T_ext << 0.0, 0.0, 0.0;

    return;
} // NetExternalForceTorque

void
ConstrainedNodalVel(Vec /*U_k*/, const RigidDOFVector& /*U*/, const Eigen::Vector3d& /*X_com*/, void* /*ctx*/)
{
    // intentionally left blank
    return;
} // ConstrainedNodalVel

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                 LDataManager* l_data_manager,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

// Check that the coupling operators of the homogeneous CIB Stokes operator are
// adjoint: with constraint force L, velocity u, interpolation J, and spreading
// S, (L, J u) = h^d (S L, u), where the Eulerian inner product counts each
// face once. The constraint row of A [u; 0; 0] is -beta J u and the momentum
// row of A [0; L; 0] is -gamma S L, with beta and gamma the interpolation and
// spreading scale factors.
void
check_coupling_adjoint(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                       Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator,
                       Pointer<IBHierarchyIntegrator> time_integrator,
                       Pointer<CIBMethod> ib_method_ops,
                       const vector<RobinBcCoefStrategy<NDIM>*>& u_bc_coefs,
                       ofstream& output_file)
{
    // The fixtures use a single level, so the Eulerian inner product is taken
    // over the finest level.
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(finest_ln);
    const double current_time = time_integrator->getIntegratorTime();
    const double new_time = current_time + time_integrator->getMaximumTimeStepSize();
    const double dt = new_time - current_time;

    // Configure the operator as CIBStaggeredStokesSolver does, but with
    // homogeneous boundary conditions and non-unit scale factors.
    const double scale_interp = 2.0;
    const double scale_spread = 0.5;
    Pointer<CIBStaggeredStokesOperator> A =
        new CIBStaggeredStokesOperator("CIBStaggeredStokesOperator", ib_method_ops, /*homogeneous_bc*/ true);
    A->setInterpScaleFactor(scale_interp);
    A->setSpreadScaleFactor(scale_spread);
    A->setRegularizeMobilityFactor(0.0);
    A->setNormalizeSpreadForce(false);
    const StokesSpecifications* problem_coefs = navier_stokes_integrator->getStokesSpecifications();
    PoissonSpecifications U_problem_coefs("U_problem_coefs");
    U_problem_coefs.setCConstant(problem_coefs->getRho() / dt + problem_coefs->getLambda());
    U_problem_coefs.setDConstant(-problem_coefs->getMu());
    A->setVelocityPoissonSpecifications(U_problem_coefs);
    A->setPhysicalBcCoefs(navier_stokes_integrator->getVelocityBoundaryConditions(),
                          navier_stokes_integrator->getPressureBoundaryConditions());
    Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();
    bc_helper->cacheBcCoefData(u_bc_coefs, new_time, patch_hierarchy);
    A->setPhysicalBoundaryHelper(bc_helper);
    A->setSolutionTime(new_time);
    A->setTimeInterval(current_time, new_time);

    // Eulerian vectors with the ghost cell width that the IB operators need.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("check_coupling_adjoint");
    Pointer<Variable<NDIM>> u_var = navier_stokes_integrator->getVelocityVariable();
    Pointer<Variable<NDIM>> p_var = navier_stokes_integrator->getPressureVariable();
    const int u_idx = var_db->registerVariableAndContext(u_var, ctx, ib_method_ops->getMinimumGhostCellWidth());
    const int p_idx = var_db->registerVariableAndContext(p_var, ctx, IntVector<NDIM>(1));
    const int A_u_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
    const int A_p_idx = var_db->registerClonedPatchDataIndex(p_var, p_idx);
    const int u_save_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
    Pointer<HierarchyMathOps> hier_math_ops = navier_stokes_integrator->getHierarchyMathOps();
    const int wgt_sc_idx = hier_math_ops->getSideWeightPatchDescriptorIndex();
    const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();
    Pointer<SAMRAIVectorReal<NDIM, double>> x =
        new SAMRAIVectorReal<NDIM, double>("x", patch_hierarchy, coarsest_ln, finest_ln);
    x->addComponent(u_var, u_idx, wgt_sc_idx);
    x->addComponent(p_var, p_idx, wgt_cc_idx);
    Pointer<SAMRAIVectorReal<NDIM, double>> y =
        new SAMRAIVectorReal<NDIM, double>("y", patch_hierarchy, coarsest_ln, finest_ln);
    y->addComponent(u_var, A_u_idx, wgt_sc_idx);
    y->addComponent(p_var, A_p_idx, wgt_cc_idx);
    x->allocateVectorData(current_time);
    y->allocateVectorData(current_time);
    x->setToScalar(0.0);
    y->setToScalar(0.0);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_save_idx, current_time);
    }
    A->initializeOperatorState(*x, *y);

    // Rigid body vectors, combined as the solver combines them.
    Vec L_ref, U_ref;
    ib_method_ops->getConstraintForce(&L_ref, current_time);
    ib_method_ops->getFreeRigidVelocities(&U_ref, current_time);
    Vec L, U, V, F;
    VecDuplicate(L_ref, &L);
    VecDuplicate(U_ref, &U);
    VecDuplicate(L, &V);
    VecDuplicate(U, &F);
    VecSet(L, 0.0);
    VecSet(U, 0.0);
    Vec u_p = PETScSAMRAIVectorReal::createPETScVector(x);
    Vec g_f = PETScSAMRAIVectorReal::createPETScVector(y);
    std::vector<Vec> vx = { u_p, L, U };
    std::vector<Vec> vy = { g_f, V, F };
    Vec mv_x, mv_y;
    VecCreateNest(PETSC_COMM_WORLD, 3, nullptr, vx.data(), &mv_x);
    VecCreateNest(PETSC_COMM_WORLD, 3, nullptr, vy.data(), &mv_y);

    // Set an arbitrary velocity that depends only on position, with zero
    // values on the Dirichlet boundary faces.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const x_lower = pgeom->getXLower();
        const double* const dx = pgeom->getDx();
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_save_idx);
        u_data->fillAll(0.0);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                double phase = 1.7 * axis;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const double X = x_lower[d] + dx[d] * (b()(d) - patch_box.lower(d) + (d == axis ? 0.0 : 0.5));
                    phase += (3.1 + 2.2 * d) * X;
                }
                u_data->getArrayData(axis)(b(), 0) = std::sin(phase);
            }
        }
    }
    bc_helper->copyDataAtDirichletBoundaries(u_save_idx, A_u_idx);

    // Constraint row of A [u; 0; 0]: V = -beta J u.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double>> u_save_data = patch->getPatchData(u_save_idx);
        u_data->copy(*u_save_data);
    }
    A->apply(mv_x, mv_y);
    Vec Ju;
    VecDuplicate(V, &Ju);
    VecCopy(V, Ju);
    VecScale(Ju, -1.0 / scale_interp);

    // Momentum row of A [0; L; 0]: A_u = -gamma S L, with an arbitrary
    // constraint force that depends only on the Lagrangian index.
    x->setToScalar(0.0);
    LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    {
        double* L_array;
        VecGetArray(L, &L_array);
        for (const auto& node : l_data_manager->getLMesh(finest_ln)->getLocalNodes())
        {
            const int lag_idx = node->getLagrangianIndex();
            const int petsc_idx = node->getLocalPETScIndex();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                L_array[petsc_idx * NDIM + d] = std::cos(1.3 * lag_idx + 0.7 * d + 0.2);
            }
        }
        VecRestoreArray(L, &L_array);
    }
    A->apply(mv_x, mv_y);
    double L_dot_Ju;
    VecDot(L, Ju, &L_dot_Ju);
    double SL_dot_u = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        double cell_volume = 1.0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            cell_volume *= dx[d];
        }
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_save_idx);
        Pointer<SideData<NDIM, double>> A_u_data = patch->getPatchData(A_u_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                // Faces on the boundary of a patch that are not on a physical
                // boundary are shared with another patch.
                const bool on_lower_face = b()(axis) == patch_box.lower(axis);
                const bool on_upper_face = b()(axis) == patch_box.upper(axis) + 1;
                const bool shared = (on_lower_face && !pgeom->getTouchesRegularBoundary(axis, 0)) ||
                                    (on_upper_face && !pgeom->getTouchesRegularBoundary(axis, 1));
                const double weight = shared ? 0.5 * cell_volume : cell_volume;
                SL_dot_u -= weight * A_u_data->getArrayData(axis)(b(), 0) * u_data->getArrayData(axis)(b(), 0);
            }
        }
    }
    SL_dot_u = IBTK_MPI::sumReduction(SL_dot_u) / scale_spread;

    if (!IBTK_MPI::getRank())
    {
        output_file << "patches on the level = " << level->getNumberOfPatches() << '\n';
        output_file << std::setprecision(12) << "(L, J u) = " << L_dot_Ju << '\n' << "(S L, u) = " << SL_dot_u << '\n';
    }

    // Cleanup.
    A->deallocateOperatorState();
    VecDestroy(&Ju);
    VecDestroy(&mv_x);
    VecDestroy(&mv_y);
    PETScSAMRAIVectorReal::destroyPETScVector(u_p);
    PETScSAMRAIVectorReal::destroyPETScVector(g_f);
    VecDestroy(&L);
    VecDestroy(&U);
    VecDestroy(&V);
    VecDestroy(&F);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(u_save_idx);
    }
    x->deallocateVectorData();
    y->deallocateVectorData();
    return;
} // check_coupling_adjoint

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
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

#ifndef IBTK_HAVE_SILO
    // Suppress warnings caused by running without Silo.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    // Several parts of the code (such as LDataManager) expect mesh files,
    // specified in the input file, to exist in the current working
    // directory. Since tests are run in temporary directories we need to regenerate these input
    // to work.
    if (IBTK_MPI::getRank() == 0)
    {
        std::ifstream plate_vertex_stream(SOURCE_DIR "/plate2d.vertex");
        std::ofstream plate_vertex_cwd("plate2d.vertex");
        plate_vertex_cwd << plate_vertex_stream.rdbuf();
    }

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "CIB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Check the coupling adjoint identity in place of time stepping.
        const bool coupling_adjoint_check = input_db->getDatabase("Main")->keyExists("check_coupling_adjoint");

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
        const unsigned int num_structures = input_db->getIntegerWithDefault("num_structures", 1);
        Pointer<CIBMethod> ib_method_ops =
            new CIBMethod("CIBMethod", app_initializer->getComponentDatabase("CIBMethod"), num_structures);

        // Krylov solver for INS integrator that solves for [u,p,U,L]
        Pointer<CIBStaggeredStokesSolver> CIBSolver =
            new CIBStaggeredStokesSolver("CIBStaggeredStokesSolver",
                                         input_db->getDatabase("CIBStaggeredStokesSolver"),
                                         navier_stokes_integrator,
                                         ib_method_ops,
                                         "SP_");

        // Register the Krylov solver with INS integrator
        navier_stokes_integrator->setStokesSolver(CIBSolver);

        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        // Configure the IB solver.
        Pointer<IBStandardInitializer> ib_initializer = new IBStandardInitializer(
            "IBStandardInitializer", app_initializer->getComponentDatabase("IBStandardInitializer"));
        ib_method_ops->registerLInitStrategy(ib_initializer);

        // Specify structure kinematics
        FreeRigidDOFVector plate_free_dofs;
        plate_free_dofs << 0, 0, 0;
        ib_method_ops->setSolveRigidBodyVelocity(0, plate_free_dofs);

        ib_method_ops->registerExternalForceTorqueFunction(&NetExternalForceTorque, nullptr, 0);
        ib_method_ops->registerConstrainedVelocityFunction(nullptr, &ConstrainedCOMVel, nullptr, 0);

        // Create initial condition specification objects.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM>> visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
            ib_method_ops->registerVisItDataWriter(visit_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Create boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = nullptr;
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

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Set physical boundary operator used in spreading.
        ib_method_ops->setVelocityPhysBdryOp(time_integrator->getVelocityPhysBdryOp());

        // Register mobility matrices (if needed)
        std::string mobility_solver_type = input_db->getString("MOBILITY_SOLVER_TYPE");
        if (mobility_solver_type == "DIRECT")
        {
            std::string mat_name = "plate";
            std::vector<std::vector<unsigned>> struct_ids;
            std::vector<unsigned> prototype_structs;

            // Dense matrix type
            prototype_structs.push_back(0);

            // Dense matrix to operate upon
            struct_ids.push_back(prototype_structs);

            // Register the dense matrix with direct solver
            DirectMobilitySolver* direct_solvers = nullptr;
            CIBSolver->getSaddlePointSolver()->getCIBMobilitySolver()->getMobilitySolvers(
                nullptr, &direct_solvers, nullptr);

            direct_solvers->registerMobilityMat(
                mat_name, prototype_structs, EMPIRICAL, std::make_pair(LAPACK_LU, LAPACK_LU), 0);
            direct_solvers->registerStructIDsWithMobilityMat(mat_name, struct_ids);
        }
        navier_stokes_integrator->setStokesSolverNeedsInit();

        // Deallocate initialization objects.
        app_initializer.setNull();

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
        if (dump_postproc_data)
        {
            output_data(patch_hierarchy,
                        ib_method_ops->getLDataManager(),
                        iteration_num,
                        loop_time,
                        postproc_data_dump_dirname);
        }

        // File to write the output data
        ofstream output_file;
        if (!IBTK_MPI::getRank()) output_file.open("output");

        if (coupling_adjoint_check)
        {
            // As advanceHierarchy() does at the initial time, regrid first,
            // which also sets up the Lagrangian indices in the ghost cells that
            // spreading uses; the boundary conditions and the Lagrangian data
            // of the time step are set up at the start of the time step.
            time_integrator->regridHierarchy();
            time_integrator->preprocessIntegrateHierarchy(
                loop_time, loop_time + time_integrator->getMaximumTimeStepSize(), 1);
            check_coupling_adjoint(
                patch_hierarchy, navier_stokes_integrator, time_integrator, ib_method_ops, u_bc_coefs, output_file);
        }

        // Main time step loop (skipped by the coupling adjoint check).
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;

        while (!coupling_adjoint_check && !IBTK::rel_equal_eps(loop_time, loop_time_end) &&
               time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();

            pout << "Advancing hierarchy by timestep size dt = " << dt << "\n";
            if (loop_time > 0.11353)
            {
                IBAMR_DO_ONCE(plate_free_dofs[0] = 1; plate_free_dofs[1] = 1; plate_free_dofs[2] = 0;
                              ib_method_ops->setSolveRigidBodyVelocity(0, plate_free_dofs);
                              navier_stokes_integrator->setStokesSolverNeedsInit(););
            }

            if (time_integrator->atRegridPoint()) navier_stokes_integrator->setStokesSolverNeedsInit();
            if (ib_method_ops->flagRegrid())
            {
                time_integrator->regridHierarchy();
                navier_stokes_integrator->setStokesSolverNeedsInit();
            }
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n\nNet rigid force and torque on plate is : \n"
                 << ib_method_ops->getNetRigidGeneralizedForce(0) << "\n\n";

            if (!IBTK_MPI::getRank())
            {
                output_file << std::setprecision(13) << "\n\nNet rigid force and torque on plate is : \n"
                            << ib_method_ops->getNetRigidGeneralizedForce(0) << "\n\n";
            }

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
                silo_data_writer->writePlotData(iteration_num, loop_time);
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
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                output_data(patch_hierarchy,
                            ib_method_ops->getLDataManager(),
                            iteration_num,
                            loop_time,
                            postproc_data_dump_dirname);
            }
        }

        // Close file
        if (!IBTK_MPI::getRank()) output_file.close();

        // Cleanup boundary condition specification objects (when necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown
} // main

void
output_data(Pointer<PatchHierarchy<NDIM>> /*patch_hierarchy*/,
            LDataManager* /*l_data_manager*/,
            const int iteration_num,
            const double loop_time,
            const string& /*data_dump_dirname*/)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    return;
} // output_data
