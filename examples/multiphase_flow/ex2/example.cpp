// ---------------------------------------------------------------------
//
// Copyright (c) 2018 - 2024 by the IBAMR developers
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
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/matlab_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AdvDiffPredictorCorrectorHierarchyIntegrator.h>
#include <ibamr/AdvDiffSemiImplicitHierarchyIntegrator.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEDirectForcingKinematics.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/IBHydrodynamicSurfaceForceEvaluator.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSVCStaggeredConservativeHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/LevelSetUtilities.h>
#include <ibamr/RelaxationLSMethod.h>
#include <ibamr/SurfaceTensionForceFunction.h>
#include <ibamr/vc_ins_utilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridFunctionSet.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibamr/app_namespaces.h>

// Application
#include "LSLocateCircularInterface.h"
#include "LSLocateGasInterface.h"

static double shift_x, shift_y;
void
coordinate_mapping_function(libMesh::Point& X, const libMesh::Point& s, void* /*ctx*/)
{
    X(0) = s(0) + shift_x;
    X(1) = s(1) + shift_y;
    ;
    return;
} // coordinate_mapping_function

void
cylinder_kinematics(double /*data_time*/, Eigen::Vector3d& U_com, Eigen::Vector3d& W_com, void* /*ctx*/)
{
    U_com.setZero();
    W_com.setZero();

    U_com[1] = -1.0;

    return;
} // cylinder_kinematics

// Function prototypes
void output_data(SAMRAIPointer<PatchHierarchyNd> patch_hierarchy,
                 SAMRAIPointer<INSVCStaggeredHierarchyIntegrator> navier_stokes_integrator,
                 LDataManager* l_data_manager,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

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
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();
#ifdef LIBMESH_HAVE_EXODUS_API
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
#else
        const bool uses_exodus = false;
        if (!app_initializer->getExodusIIFilename().empty())
        {
            plog << "WARNING: libMesh was compiled without Exodus support, so no "
                 << "Exodus output will be written in this program.\n";
        }
#endif
        const string exodus_filename = app_initializer->getExodusIIFilename();

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

        // Setup solid information
        CircularInterface circle;
        circle.R = input_db->getDouble("R");
        circle.X0[0] = input_db->getDouble("XCOM");
        circle.X0[1] = input_db->getDouble("YCOM");
#if (NDIM == 3)
        circle.X0[2] = input_db->getDouble("ZCOM");
#endif

        // Create a simple FE mesh.
        Mesh solid_mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        shift_x = circle.X0[0];
        shift_y = circle.X0[1];
        if (input_db->keyExists("XDA_FILENAME"))
        {
            TBOX_ASSERT(elem_type == "TRI3" || elem_type == "TRI6");

            std::string filename = input_db->getString("XDA_FILENAME");
            MatlabIO distmesh(solid_mesh);
            distmesh.read(filename);

            if (elem_type == "TRI6") solid_mesh.all_second_order();
        }
        else if (elem_type == "TRI3" || elem_type == "TRI6")
        {
#ifdef LIBMESH_HAVE_TRIANGLE
            const int num_circum_nodes = ceil(2.0 * M_PI * circle.R / ds);
            for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                solid_mesh.add_point(libMesh::Point(circle.R * cos(theta), circle.R * sin(theta)));
            }
            TriangleInterface triangle(solid_mesh);
            triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
            triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();

            if (elem_type == "TRI6") solid_mesh.all_second_order();
#else
            TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                       << "       but Triangle is required for TRI3 or TRI6 elements.\n");
#endif
        }
        else
        {
            // NOTE: number of segments along boundary is 4*2^r.
            const double num_circum_segments = 2.0 * M_PI * circle.R / ds;
            const int r = log2(0.25 * num_circum_segments);
            MeshTools::Generation::build_sphere(solid_mesh, circle.R, r, Utility::string_to_enum<ElemType>(elem_type));
        }

        solid_mesh.prepare_for_use();
        Mesh& mesh = solid_mesh;

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        SAMRAIPointer<INSVCStaggeredHierarchyIntegrator> navier_stokes_integrator;
        const string discretization_form =
            app_initializer->getComponentDatabase("Main")->getString("discretization_form");
        const bool conservative_form = (discretization_form == "CONSERVATIVE");
        if (conservative_form)
        {
            navier_stokes_integrator = new INSVCStaggeredConservativeHierarchyIntegrator(
                "INSVCStaggeredConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredConservativeHierarchyIntegrator"));
        }
        else if (!conservative_form)
        {
            navier_stokes_integrator = new INSVCStaggeredNonConservativeHierarchyIntegrator(
                "INSVCStaggeredNonConservativeHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSVCStaggeredNonConservativeHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << discretization_form << "\n"
                                                   << "Valid options are: CONSERVATIVE, NON_CONSERVATIVE");
        }

        // Set up the advection diffusion hierarchy integrator
        SAMRAIPointer<AdvDiffHierarchyIntegrator> adv_diff_integrator;
        const string adv_diff_solver_type = app_initializer->getComponentDatabase("Main")->getStringWithDefault(
            "adv_diff_solver_type", "PREDICTOR_CORRECTOR");
        if (adv_diff_solver_type == "PREDICTOR_CORRECTOR")
        {
            SAMRAIPointer<AdvectorExplicitPredictorPatchOps> predictor = new AdvectorExplicitPredictorPatchOps(
                "AdvectorExplicitPredictorPatchOps",
                app_initializer->getComponentDatabase("AdvectorExplicitPredictorPatchOps"));
            adv_diff_integrator = new AdvDiffPredictorCorrectorHierarchyIntegrator(
                "AdvDiffPredictorCorrectorHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffPredictorCorrectorHierarchyIntegrator"),
                predictor);
        }
        else if (adv_diff_solver_type == "SEMI_IMPLICIT")
        {
            adv_diff_integrator = new AdvDiffSemiImplicitHierarchyIntegrator(
                "AdvDiffSemiImplicitHierarchyIntegrator",
                app_initializer->getComponentDatabase("AdvDiffSemiImplicitHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << adv_diff_solver_type << "\n"
                                                   << "Valid options are: PREDICTOR_CORRECTOR, SEMI_IMPLICIT");
        }
        navier_stokes_integrator->registerAdvDiffHierarchyIntegrator(adv_diff_integrator);

        SAMRAIPointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        SAMRAIPointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);

        SAMRAIPointer<CartesianGridGeometryNd> grid_geometry = new CartesianGridGeometryNd(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<PatchHierarchyNd> patch_hierarchy = new PatchHierarchyNd("PatchHierarchy", grid_geometry);

        SAMRAIPointer<StandardTagAndInitializeNd> error_detector =
            new StandardTagAndInitializeNd("StandardTagAndInitialize",
                                           time_integrator,
                                           app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<BergerRigoutsosNd> box_generator = new BergerRigoutsosNd();
        SAMRAIPointer<LoadBalancerNd> load_balancer =
            new LoadBalancerNd("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<GriddingAlgorithmNd> gridding_algorithm =
            new GriddingAlgorithmNd("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

        // Create level sets for solid interface.
        const string& ls_name_solid = "level_set_solid";
        SAMRAIPointer<CellVariableNd<double> > phi_var_solid = new CellVariableNd<double>(ls_name_solid);
        SAMRAIPointer<RelaxationLSMethod> level_set_solid_ops =
            new RelaxationLSMethod(ls_name_solid, app_initializer->getComponentDatabase("LevelSet_Solid"));
        LSLocateCircularInterface setLSLocateCircularInterface(
            "LSLocateCircularInterface", adv_diff_integrator, phi_var_solid, &circle);
        level_set_solid_ops->registerInterfaceNeighborhoodLocatingFcn(
            &callLSLocateCircularInterfaceCallbackFunction, static_cast<void*>(&setLSLocateCircularInterface));

        // Create level sets for gas/liquid interface.
        const double fluid_height = input_db->getDouble("GAS_LS_INIT");
        const string& ls_name_gas = "level_set_gas";
        SAMRAIPointer<CellVariableNd<double> > phi_var_gas = new CellVariableNd<double>(ls_name_gas);
        SAMRAIPointer<RelaxationLSMethod> level_set_gas_ops =
            new RelaxationLSMethod(ls_name_gas, app_initializer->getComponentDatabase("LevelSet_Gas"));
        LSLocateGasInterface setLSLocateGasInterface(
            "LSLocateGasInterface", adv_diff_integrator, phi_var_gas, fluid_height);
        level_set_gas_ops->registerInterfaceNeighborhoodLocatingFcn(&callLSLocateGasInterfaceCallbackFunction,
                                                                    static_cast<void*>(&setLSLocateGasInterface));

        // Register the level sets with advection diffusion integrator.
        adv_diff_integrator->registerTransportedQuantity(phi_var_solid);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_solid, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_solid,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        adv_diff_integrator->registerTransportedQuantity(phi_var_gas);
        adv_diff_integrator->setDiffusionCoefficient(phi_var_gas, 0.0);
        adv_diff_integrator->setAdvectionVelocity(phi_var_gas,
                                                  navier_stokes_integrator->getAdvectionVelocityVariable());

        // Register the reinitialization functions for the level set variables
        IBAMR::LevelSetUtilities::SetLSProperties setSetLSPropertiesSolid("SetLSPropertiesSolid", level_set_solid_ops);
        adv_diff_integrator->registerResetFunction(phi_var_solid,
                                                   &IBAMR::LevelSetUtilities::setLSDataPatchHierarchy,
                                                   static_cast<void*>(&setSetLSPropertiesSolid));

        IBAMR::LevelSetUtilities::SetLSProperties setSetLSPropertiesGas("SetLSPropertiesGas", level_set_gas_ops);
        adv_diff_integrator->registerResetFunction(phi_var_gas,
                                                   &IBAMR::LevelSetUtilities::setLSDataPatchHierarchy,
                                                   static_cast<void*>(&setSetLSPropertiesGas));

        // LS initial conditions
        if (input_db->keyExists("LevelSetGasInitialConditions"))
        {
            SAMRAIPointer<CartGridFunction> phi_init_gas = new muParserCartGridFunction(
                "phi_init_gas", app_initializer->getComponentDatabase("LevelSetGasInitialConditions"), grid_geometry);
            adv_diff_integrator->setInitialConditions(phi_var_gas, phi_init_gas);
        }
        if (input_db->keyExists("LevelSetSolidInitialConditions"))
        {
            SAMRAIPointer<CartGridFunction> phi_init_solid =
                new muParserCartGridFunction("phi_init_solid",
                                             app_initializer->getComponentDatabase("LevelSetSolidInitialConditions"),
                                             grid_geometry);
            adv_diff_integrator->setInitialConditions(phi_var_solid, phi_init_solid);
        }

        // Setup the advected and diffused fluid quantities.
        SAMRAIPointer<CellVariableNd<double> > mu_var = new CellVariableNd<double>("mu");
        SAMRAIPointer<hier::VariableNd> rho_var;
        if (conservative_form)
        {
            rho_var = new SideVariableNd<double>("rho");
        }
        else
        {
            rho_var = new CellVariableNd<double>("rho");
        }

        navier_stokes_integrator->registerMassDensityVariable(rho_var);
        navier_stokes_integrator->registerViscosityVariable(mu_var);

        // Array for input into callback function
        const double rho_fluid = input_db->getDouble("RHO_F");
        const double rho_solid = input_db->getDouble("RHO_S");
        const double rho_gas = input_db->getDouble("RHO_G");
        const double mu_fluid = input_db->getDouble("MU_F");
        const double mu_gas = input_db->getDouble("MU_G");
        const double mu_solid = input_db->getDoubleWithDefault("MU_S", std::numeric_limits<double>::quiet_NaN());
        const bool set_mu_solid = input_db->getBool("SET_MU_S");
        const int num_solid_interface_cells = input_db->getDouble("NUM_SOLID_INTERFACE_CELLS");
        const int num_gas_interface_cells = input_db->getDouble("NUM_GAS_INTERFACE_CELLS");
        IBAMR::VCINSUtilities::SetFluidProperties setSetFluidProperties("SetFluidProperties",
                                                                        adv_diff_integrator,
                                                                        phi_var_gas,
                                                                        phi_var_solid,
                                                                        rho_fluid,
                                                                        rho_gas,
                                                                        rho_solid,
                                                                        mu_fluid,
                                                                        mu_gas,
                                                                        mu_solid,
                                                                        num_gas_interface_cells,
                                                                        num_solid_interface_cells,
                                                                        set_mu_solid);
        navier_stokes_integrator->registerResetFluidDensityFcn(&IBAMR::VCINSUtilities::callSetDensityCallbackFunction,
                                                               static_cast<void*>(&setSetFluidProperties));
        navier_stokes_integrator->registerResetFluidViscosityFcn(
            &IBAMR::VCINSUtilities::callSetViscosityCallbackFunction, static_cast<void*>(&setSetFluidProperties));

        // Register callback function for tagging refined cells for level set data
        const double tag_thresh = input_db->getDouble("LS_TAG_ABS_THRESH");
        const double tag_min_value = -tag_thresh;
        const double tag_max_value = tag_thresh;
        IBAMR::LevelSetUtilities::TagLSRefinementCells ls_tagger(
            adv_diff_integrator, phi_var_gas, tag_min_value, tag_max_value);
        time_integrator->registerApplyGradientDetectorCallback(&IBAMR::LevelSetUtilities::tagLSCells,
                                                               static_cast<void*>(&ls_tagger));

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            SAMRAIPointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            SAMRAIPointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVectorNd& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategyNd*> u_bc_coefs(NDIM);
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

        RobinBcCoefStrategyNd* rho_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("DensityBcCoefs"))
        {
            rho_bc_coef = new muParserRobinBcCoefs(
                "rho_bc_coef", app_initializer->getComponentDatabase("DensityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerMassDensityBoundaryConditions(rho_bc_coef);
        }

        RobinBcCoefStrategyNd* mu_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ViscosityBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("ViscosityBcCoefs"), grid_geometry);
            navier_stokes_integrator->registerViscosityBoundaryConditions(mu_bc_coef);
        }

        RobinBcCoefStrategyNd* phi_bc_coef = NULL;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("PhiBcCoefs"))
        {
            phi_bc_coef = new muParserRobinBcCoefs(
                "phi_bc_coef", app_initializer->getComponentDatabase("PhiBcCoefs"), grid_geometry);
        }
        adv_diff_integrator->setPhysicalBcCoef(phi_var_gas, phi_bc_coef);
        adv_diff_integrator->setPhysicalBcCoef(phi_var_solid, phi_bc_coef);

        // LS reinit boundary conditions, which is set to be the same as the BCs
        // for advection
        RobinBcCoefStrategyNd* ls_reinit_bcs = phi_bc_coef;
        level_set_solid_ops->registerPhysicalBoundaryCondition(ls_reinit_bcs);
        level_set_gas_ops->registerPhysicalBoundaryCondition(ls_reinit_bcs);

        // Initialize objects
        std::vector<double> grav_const(NDIM);
        input_db->getDoubleArray("GRAV_CONST", &grav_const[0], NDIM);
        SAMRAIPointer<CartGridFunction> grav_force =
            new IBAMR::VCINSUtilities::GravityForcing("GravityForcing", navier_stokes_integrator, grav_const);

        SAMRAIPointer<SurfaceTensionForceFunction> surface_tension_force =
            new SurfaceTensionForceFunction("SurfaceTensionForceFunction",
                                            app_initializer->getComponentDatabase("SurfaceTensionForceFunction"),
                                            adv_diff_integrator,
                                            phi_var_gas);

        SAMRAIPointer<CartGridFunctionSet> eul_forces = new CartGridFunctionSet("eulerian_forces");
        eul_forces->addFunction(grav_force);
        eul_forces->addFunction(surface_tension_force);
        time_integrator->registerBodyForceFunction(eul_forces);

        // Create IBFE direct forcing kinematics object.
        SAMRAIPointer<IBFEDirectForcingKinematics> df_kinematics_ops = new IBFEDirectForcingKinematics(
            "cylinder_dfk",
            app_initializer->getComponentDatabase("CylinderIBFEDirectForcingKinematics"),
            ib_method_ops,
            /*part*/ 0,
            /*register_for_restart*/ true);
        ib_method_ops->registerDirectForcingKinematics(df_kinematics_ops, /*part*/ 0);

        // Specify structure kinematics
        FreeRigidDOFVector solve_dofs;
        solve_dofs.setZero();
        input_db->getIntegerArray("FREE_DOFS", &solve_dofs[0], s_max_free_dofs);
        df_kinematics_ops->setSolveRigidBodyVelocity(solve_dofs);
        df_kinematics_ops->registerKinematicsFunction(&cylinder_kinematics, NULL);

        // Create the surface force evaluator object
        SAMRAIPointer<IBHydrodynamicSurfaceForceEvaluator> hydro_force_evaluator =
            new IBHydrodynamicSurfaceForceEvaluator("Cylinder_Hydroforce",
                                                    phi_var_solid,
                                                    adv_diff_integrator,
                                                    navier_stokes_integrator,
                                                    input_db->getDatabase("IBHydrodynamicSurfaceForceEvaluator"));
        IBTK::Vector3d pressure_force, viscous_force, pressure_torque, viscous_torque;
        hydro_force_evaluator->writeToFile();

        // Configure the IBFE solver.
        ib_method_ops->initializeFEEquationSystems();
        ib_method_ops->registerInitialCoordinateMappingFunction(coordinate_mapping_function);
        EquationSystems* equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

        // Set up visualization plot file writers.
        SAMRAIPointer<VisItDataWriterNd> visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                exodus_io->write_timestep(
                    exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }

        // File to write to for fluid mass data
        ofstream mass_file;
        if (!IBTK_MPI::getRank()) mass_file.open("mass_fluid.txt");
        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            pout << "Advancing hierarchy with timestep size dt = " << dt << "\n";
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            // Update the center of mass of the circle.
            const Eigen::Vector3d& structure_COM = df_kinematics_ops->getStructureCOM();
            circle.X0(0) = structure_COM[0];
            circle.X0(1) = structure_COM[1];
#if (NDIM == 3)
            circle.X0(2) = structure_COM[2];
#endif

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // Compute and print the hydrodynamic force and torque.
            hydro_force_evaluator->computeHydrodynamicForceTorque(
                pressure_force, viscous_force, pressure_torque, viscous_torque, circle.X0);

            // Compute the fluid mass in the domain from interpolated density
            const int rho_ins_idx = navier_stokes_integrator->getLinearOperatorRhoPatchDataIndex();
#if !defined(NDEBUG)
            TBOX_ASSERT(rho_ins_idx >= 0);
#endif
            const int coarsest_ln = 0;
            const int finest_ln = patch_hierarchy->getFinestLevelNumber();
            HierarchySideDataOpsRealNd<double> hier_rho_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
            HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
            hier_math_ops.setPatchHierarchy(patch_hierarchy);
            hier_math_ops.resetLevels(coarsest_ln, finest_ln);
            const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();
            const double mass_fluid = hier_rho_data_ops.integral(rho_ins_idx, wgt_sc_idx);

            // Write to file
            if (!IBTK_MPI::getRank())
            {
                mass_file << std::setprecision(13) << loop_time << "\t" << mass_fluid << std::endl;
            }

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    exodus_io->write_timestep(
                        exodus_filename, *equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                }
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

        // Close file
        if (!IBTK_MPI::getRank()) mass_file.close();

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        delete rho_bc_coef;
        delete mu_bc_coef;
        delete phi_bc_coef;

    } // cleanup dynamically allocated objects prior to shutdown
} // main
