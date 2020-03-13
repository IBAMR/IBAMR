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

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

#include <array>

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 LDataManager* l_data_manager,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

int finest_ln;
std::array<int, NDIM> N;
SAMRAI::tbox::Array<std::string> struct_list;
SAMRAI::tbox::Array<int> num_node;
SAMRAI::tbox::Array<double> ds;
int num_node_circum, num_node_radial;
double dr;
void
generate_structure(const unsigned int& strct_num,
                   const int& ln,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn)
{
    if (ln != finest_ln)
    {
        num_vertices = 0;
        vertex_posn.resize(num_vertices);
    }

    double R = 0.25;             // Radius of circle
    double perim = 2 * M_PI * R; // Perimenter of disc
    double dx = 1.0 / static_cast<double>(N[0]);
    if (struct_list[strct_num].compare("curve2d") == 0)
    {
        num_node[strct_num] = static_cast<int>(perim / dx);
        ds[strct_num] = 2.0 * M_PI * R / num_node[strct_num];
        num_vertices = num_node[strct_num];
        vertex_posn.resize(num_vertices);
        for (std::vector<IBTK::Point>::iterator it = vertex_posn.begin(); it != vertex_posn.end(); ++it)
        {
            Point& X = *it;
            int num = std::distance(vertex_posn.begin(), it);
            double theta = 2.0 * M_PI * num / num_vertices;
            X(0) = 0.5 + R * std::cos(theta);
            X(1) = 0.5 + R * std::sin(theta);
        }
    }
    return;
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

    // resize error vectors to hold data from u and p

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type =
            app_initializer->getComponentDatabase("Main")->getStringWithDefault("solver_type", "STAGGERED");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
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
        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        struct_list.resizeArray(input_db->getArraySize("STRUCTURE_LIST"));
        struct_list = input_db->getStringArray("STRUCTURE_LIST");
        std::vector<std::string> struct_list_vec(struct_list.getSize());
        for (int i = 0; i < struct_list.size(); ++i) struct_list_vec[i] = struct_list[i];
        ds.resizeArray(struct_list.size());
        num_node.resizeArray(struct_list.getSize());
        N[0] = N[1] = input_db->getInteger("N");
        finest_ln = input_db->getInteger("MAX_LEVELS") - 1;
        ib_initializer->setStructureNamesOnLevel(finest_ln, struct_list_vec);
        ib_initializer->registerInitStructureFunction(generate_structure);
        ib_method_ops->registerLInitStrategy(ib_initializer);

        // Create Eulerian initial condition specification objects.  These
        // objects also are used to specify exact solution values for error
        // analysis.
        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerVelocityInitialConditions(u_init);

        Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
            "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
        navier_stokes_integrator->registerPressureInitialConditions(p_init);

        Pointer<CartGridFunction> u_mask = new muParserCartGridFunction(
            "u_mask", app_initializer->getComponentDatabase("VelocityMask"), grid_geometry);
        Pointer<CartGridFunction> p_mask = new muParserCartGridFunction(
            "p_mask", app_initializer->getComponentDatabase("PressureMask"), grid_geometry);

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
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);

                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
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
        const int u_interp_idx = var_db->registerVariableAndContext(u_var, main_ctx, IntVector<NDIM>(n_ghosts));
        const int u_masked_idx = var_db->registerClonedPatchDataIndex(u_var, u_interp_idx);

        const Pointer<Variable<NDIM> > p_var = navier_stokes_integrator->getPressureVariable();
        const int p_interp_idx = var_db->registerVariableAndContext(p_var, main_ctx, IntVector<NDIM>(n_ghosts));
        const int p_masked_idx = var_db->registerClonedPatchDataIndex(p_var, p_interp_idx);
        visit_data_writer->registerPlotQuantity("P interp", "SCALAR", p_interp_idx);
        visit_data_writer->registerPlotQuantity("P masked", "SCALAR", p_masked_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_interp_idx);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p_interp_idx);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(u_masked_idx);
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(p_masked_idx);
        }

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();

        u_init->setDataOnPatchHierarchy(u_interp_idx, u_var, patch_hierarchy, loop_time);
        p_init->setDataOnPatchHierarchy(p_interp_idx, p_var, patch_hierarchy, loop_time);
        u_mask->setDataOnPatchHierarchy(u_masked_idx, u_var, patch_hierarchy, loop_time);
        p_mask->setDataOnPatchHierarchy(p_masked_idx, p_var, patch_hierarchy, loop_time);

        // Fill ghost data.
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> ghost_cell_components(4);
        ghost_cell_components[0] = ITC(u_interp_idx,
                                       "CONSERVATIVE_LINEAR_REFINE",
                                       true,
                                       "CONSERVATIVE_COARSEN",
                                       "LINEAR",
                                       false,
                                       {}, // u_bc_coefs
                                       nullptr);
        ghost_cell_components[1] = ITC(p_interp_idx,
                                       "LINEAR_REFINE",
                                       true,
                                       "CONSERVATIVE_COARSEN",
                                       "LINEAR",
                                       false,
                                       {}, // p_bc_coefs
                                       nullptr);
        ghost_cell_components[2] = ITC(u_masked_idx,
                                       "CONSERVATIVE_LINEAR_REFINE",
                                       true,
                                       "CONSERVATIVE_COARSEN",
                                       "LINEAR",
                                       false,
                                       {}, // u_bc_coefs
                                       nullptr);
        ghost_cell_components[3] = ITC(p_masked_idx,
                                       "LINEAR_REFINE",
                                       true,
                                       "CONSERVATIVE_COARSEN",
                                       "LINEAR",
                                       false,
                                       {}, // p_bc_coefs
                                       nullptr);
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
        ghost_fill_op.fillData(/*time*/ 0.0);

        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
        }

        LDataManager* l_data_manager = ib_method_ops->getLDataManager();
        Pointer<IBTK::LData> X_data = l_data_manager->getLData("X", finest_ln);
        Pointer<IBTK::LData> exact_U_data = l_data_manager->createLData("exact_U",
                                                                        finest_ln,
                                                                        NDIM,
                                                                        /*manage_data*/ true);
        Pointer<IBTK::LData> interp_U_data = l_data_manager->createLData("interp_U",
                                                                         finest_ln,
                                                                         NDIM,
                                                                         /*manage_data*/ true);
        Pointer<IBTK::LData> exact_P_data = l_data_manager->createLData("exact_P",
                                                                        finest_ln,
                                                                        1,
                                                                        /*manage_data*/ true);
        Pointer<IBTK::LData> interp_P_data = l_data_manager->createLData("interp_P",
                                                                         finest_ln,
                                                                         1,
                                                                         /*manage_data*/ true);

        // Set exact value for Lagrangian data
        const boost::multi_array_ref<double, 2>& X_array = *(X_data->getLocalFormVecArray());
        boost::multi_array_ref<double, 2>& U_array = *(exact_U_data->getLocalFormVecArray());
        boost::multi_array_ref<double, 2>& P_array = *(exact_P_data->getLocalFormVecArray());

        const Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();
        for (const auto& node_idx : local_nodes)
        {
            const int local_idx = node_idx->getLocalPETScIndex();
            const double* const X = &(X_array[local_idx][0]);
            double* const U = &(U_array[local_idx][0]);
            double* const P = &(P_array[local_idx][0]);

            for (int d = 0; d < NDIM; ++d)
            {
                U[d] = 1.0 + 3.5 * X[0] + 4.0 * X[1];
            }
            P[0] = 1.0 + 3.5 * X[0] + 4.0 * X[1];
        }

        // Interpolate Eulerian data
        l_data_manager->interp(u_interp_idx, u_masked_idx, interp_U_data, X_data, weighting_fcn, finest_ln);
        l_data_manager->interp(p_interp_idx, p_masked_idx, interp_P_data, X_data, weighting_fcn, finest_ln);

        // Find the error
        Vec vec_U_exact = exact_U_data->getVec();
        Vec vec_U_interp = interp_U_data->getVec();
        VecAYPX(vec_U_exact, -1.0, vec_U_interp);
        PetscScalar norm_1, norm_2, norm_max;
        VecNorm(vec_U_exact, NORM_1, &norm_1);
        VecNorm(vec_U_exact, NORM_2, &norm_2);
        VecNorm(vec_U_exact, NORM_INFINITY, &norm_max);

        pout << "Error in interpolated u:\n"
             << "  L1-norm:  " << std::setprecision(10) << norm_1 << "\n"
             << "  L2-norm:  " << norm_2 << "\n"
             << "  max-norm: " << norm_max << "\n";

        Vec vec_P_exact = exact_P_data->getVec();
        Vec vec_P_interp = interp_P_data->getVec();
        VecAYPX(vec_P_exact, -1.0, vec_P_interp);
        VecNorm(vec_P_exact, NORM_1, &norm_1);
        VecNorm(vec_P_exact, NORM_2, &norm_2);
        VecNorm(vec_P_exact, NORM_INFINITY, &norm_max);

        pout << "Error in interpolated p:\n"
             << "  L1-norm:  " << std::setprecision(10) << norm_1 << "\n"
             << "  L2-norm:  " << norm_2 << "\n"
             << "  max-norm: " << norm_max << "\n";

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
} // main

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
            LDataManager* l_data_manager,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();

    // Write Lagrangian data.
    const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
    Pointer<LData> X_data = l_data_manager->getLData("X", finest_hier_level);
    Vec X_petsc_vec = X_data->getVec();
    Vec X_lag_vec;
    VecDuplicate(X_petsc_vec, &X_lag_vec);
    l_data_manager->scatterPETScToLagrangian(X_petsc_vec, X_lag_vec, finest_hier_level);
    file_name = data_dump_dirname + "/" + "X.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    VecView(X_lag_vec, viewer);
    PetscViewerDestroy(&viewer);
    VecDestroy(&X_lag_vec);
    return;
} // output_data
