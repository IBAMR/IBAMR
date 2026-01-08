// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
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

// Headers for basic PETSc objects
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/AcousticStreamingHierarchyIntegrator.h>
#include <ibamr/AcousticStreamingPETScMatUtilities.h>
#include <ibamr/AcousticStreamingPETScVecUtilities.h>
#include <ibamr/BrinkmanPenalizationMethod.h>
#include <ibamr/FOAcousticStreamingPETScLevelSolver.h>
#include <ibamr/SOAcousticStreamingBrinkmanPenalization.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <sstream>

#include <ibamr/namespaces.h>

// Routines to reset fluid properties
void
callSetFluidPropertyCallbackFunction(int dst_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > dst_var,
                                     SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                     const int /*cycle_num*/,
                                     const double time,
                                     const double /*current_time*/,
                                     const double /*new_time*/,
                                     void* ctx)
{
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    CartGridFunction* cart_fcn = static_cast<CartGridFunction*>(ctx);
    cart_fcn->setDataOnPatchHierarchy(dst_idx, dst_var, patch_hierarchy, time);
    return;
} // callSetFluidPropertyCallbackFunction

// Struct to reset solid level set
struct SolidLevelSetResetter
{
    Pointer<HierarchyIntegrator> time_integrator;
    int ls_idx;

    // "CYLINDER" geometry
    std::array<double, NDIM> center;
    double R;
};

// Struct to reset contour level set
struct ContourLevelSetResetter
{
    Pointer<HierarchyIntegrator> time_integrator;
    int contour_idx;

    // "Circular" geometry
    std::array<double, NDIM> center;
    double R;
};

void
reset_solid_level_set_callback_fcn(double /*current_time*/, double /*new_time*/, int /*cycle_num*/, void* ctx)
{
    SolidLevelSetResetter* resetter = static_cast<SolidLevelSetResetter*>(ctx);

    Pointer<PatchHierarchy<NDIM> > patch_hier = resetter->time_integrator->getPatchHierarchy();
    const int hier_finest_ln = patch_hier->getFinestLevelNumber();

    const int& patch_idx = resetter->ls_idx;

    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();

            Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(patch_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                const hier::Index<NDIM>& ci = it();
                IBTK::VectorNd coord = IBTK::VectorNd::Zero();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }

                const double distance = std::sqrt(std::pow(coord[0] - resetter->center[0], 2) +
                                                  std::pow(coord[1] - resetter->center[1], 2)) -
                                        resetter->R;
                (*ls_data)(ci) = distance;
            }
        }
    }

    return;
} // reset_solid_level_set_callback_fcn

void
reset_contour_level_set_callback_fcn(double /*current_time*/, double /*new_time*/, int /*num_cycles*/, void* ctx)
{
    ContourLevelSetResetter* resetter = static_cast<ContourLevelSetResetter*>(ctx);

    Pointer<PatchHierarchy<NDIM> > patch_hier = resetter->time_integrator->getPatchHierarchy();
    const int hier_finest_ln = patch_hier->getFinestLevelNumber();

    const int& patch_idx = resetter->contour_idx;

    for (int ln = 0; ln <= hier_finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > patch_level = patch_hier->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(patch_level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = patch_level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::Index<NDIM>& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();

            Pointer<CellData<NDIM, double> > contour_data = patch->getPatchData(patch_idx);
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                const hier::Index<NDIM>& ci = it();
                IBTK::VectorNd coord = IBTK::VectorNd::Zero();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }
                (*contour_data)(ci) = std::sqrt(std::pow(coord[0] - resetter->center[0], 2) +
                                                std::pow(coord[1] - resetter->center[1], 2)) -
                                      resetter->R;
            }
        }
    }
    return;
} // reset_contour_level_set_callback_fcn

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer =
            new AppInitializer(argc, argv, "acoustic_streaming_hier_integrator.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const std::string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const std::string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create acoustic hierarchy integrator
        Pointer<AcousticStreamingHierarchyIntegrator> time_integrator = new AcousticStreamingHierarchyIntegrator(
            "AcousticStreamingHierarchyIntegrator",
            app_initializer->getComponentDatabase("AcousticStreamingHierarchyIntegrator"),
            /*register_for_restart*/ true);

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
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

        // Create level sets for the solid interface.
        const std::string& ls_name_solid = "level_set_solid";
        Pointer<CellVariable<NDIM, double> > phi_var_solid = new CellVariable<NDIM, double>(ls_name_solid);

        // Create initial condition specification objects.
        Pointer<CartGridFunction> u1_init = new muParserCartGridFunction(
            "u1_init", app_initializer->getComponentDatabase("FOVelocityInitialConditions"), grid_geometry);
        time_integrator->registerFirstOrderVelocityInitialConditions(u1_init);
        Pointer<CartGridFunction> p1_init = new muParserCartGridFunction(
            "p1_init", app_initializer->getComponentDatabase("FOPressureInitialConditions"), grid_geometry);
        time_integrator->registerFirstOrderPressureInitialConditions(p1_init);

        Pointer<CartGridFunction> u2_init = new muParserCartGridFunction(
            "u2_init", app_initializer->getComponentDatabase("SOVelocityInitialConditions"), grid_geometry);
        time_integrator->registerSecondOrderVelocityInitialConditions(u2_init);
        Pointer<CartGridFunction> p2_init = new muParserCartGridFunction(
            "p2_init", app_initializer->getComponentDatabase("SOPressureInitialConditions"), grid_geometry);
        time_integrator->registerSecondOrderPressureInitialConditions(p2_init);

        Pointer<CartGridFunction> rho_init =
            new muParserCartGridFunction("rho_init", app_initializer->getComponentDatabase("rho"), grid_geometry);
        time_integrator->registerMassDensityInitialConditions(rho_init);
        Pointer<CartGridFunction> mu_init =
            new muParserCartGridFunction("mu_init", app_initializer->getComponentDatabase("mu"), grid_geometry);
        time_integrator->registerShearViscosityInitialConditions(mu_init);
        Pointer<CartGridFunction> lambda_init =
            new muParserCartGridFunction("lambda_init", app_initializer->getComponentDatabase("lambda"), grid_geometry);
        time_integrator->registerBulkViscosityInitialConditions(lambda_init);

        // Create boundary condition specification objects for the first-order system.
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        std::array<std::vector<RobinBcCoefStrategy<NDIM>*>, 2> u1_bc_coefs;
        u1_bc_coefs[0].resize(NDIM);
        u1_bc_coefs[1].resize(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (int comp = 0; comp < 2; ++comp)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    u1_bc_coefs[comp][d] = nullptr;
                }
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u1_real_bc_coefs_" + std::to_string(d);
                const std::string bc_coefs_db_name = "FOVelocityRealBcCoefs_" + std::to_string(d);

                u1_bc_coefs[0][d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }

            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u1_imag_bc_coefs_" + std::to_string(d);
                const std::string bc_coefs_db_name = "FOVelocityImagBcCoefs_" + std::to_string(d);

                u1_bc_coefs[1][d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
        }
        time_integrator->registerFirstOrderPhysicalBoundaryConditions(u1_bc_coefs);

        // Boundary conditions for density and viscosity
        std::vector<RobinBcCoefStrategy<NDIM>*> rho_bc_coefs(NDIM, nullptr);
        if (!(periodic_shift.min() > 0))
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                rho_bc_coefs[d] = new muParserRobinBcCoefs(
                    "rho_bc_coef",
                    app_initializer->getComponentDatabase("DensityBcCoefs_" + std::to_string(d)),
                    grid_geometry);
            }
        }

        RobinBcCoefStrategy<NDIM>* mu_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("ShearViscosityBcCoefs"))
        {
            mu_bc_coef = new muParserRobinBcCoefs(
                "mu_bc_coef", app_initializer->getComponentDatabase("ShearViscosityBcCoefs"), grid_geometry);
        }

        RobinBcCoefStrategy<NDIM>* lambda_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("BulkViscosityBcCoefs"))
        {
            lambda_bc_coef = new muParserRobinBcCoefs(
                "lambda_bc_coef", app_initializer->getComponentDatabase("BulkViscosityBcCoefs"), grid_geometry);
        }

        // Boundary conditions for the solid level set function
        RobinBcCoefStrategy<NDIM>* ls_bc_coef = nullptr;
        if (!(periodic_shift.min() > 0) && input_db->keyExists("LevelSetBcCoefs"))
        {
            ls_bc_coef = new muParserRobinBcCoefs(
                "ls_bc_coef", app_initializer->getComponentDatabase("LevelSetBcCoefs"), grid_geometry);
        }

        // Setup the integrator maintained material properties.
        Pointer<SideVariable<NDIM, double> > rho_var = new SideVariable<NDIM, double>("rho_var");
        time_integrator->registerMassDensityVariable(rho_var);
        time_integrator->registerMassDensityBoundaryConditions(rho_bc_coefs);

        Pointer<CellVariable<NDIM, double> > mu_var = new CellVariable<NDIM, double>("mu");
        time_integrator->registerShearViscosityVariable(mu_var);
        time_integrator->registerShearViscosityBoundaryConditions(mu_bc_coef);

        Pointer<CellVariable<NDIM, double> > lambda_var = new CellVariable<NDIM, double>("lambda");
        time_integrator->registerBulkViscosityVariable(lambda_var);
        time_integrator->registerBulkViscosityBoundaryConditions(lambda_bc_coef);

        // Callback fncs to reset fluid material properties
        time_integrator->registerResetFluidDensityFcn(&callSetFluidPropertyCallbackFunction,
                                                      static_cast<void*>(rho_init.getPointer()));
        time_integrator->registerResetFluidShearViscosityFcn(&callSetFluidPropertyCallbackFunction,
                                                             static_cast<void*>(mu_init.getPointer()));
        time_integrator->registerResetFluidBulkViscosityFcn(&callSetFluidPropertyCallbackFunction,
                                                            static_cast<void*>(lambda_init.getPointer()));

        // Reset solid geometry
        SolidLevelSetResetter solid_level_set_resetter;
        solid_level_set_resetter.time_integrator = time_integrator;
        solid_level_set_resetter.center[0] = input_db->getDouble("XCOM");
        solid_level_set_resetter.center[1] = input_db->getDouble("YCOM");
        solid_level_set_resetter.R = input_db->getDouble("Radius");

        auto fcn_ptr = &reset_solid_level_set_callback_fcn;
        time_integrator->registerIntegrateHierarchyCallback(fcn_ptr, static_cast<void*>(&solid_level_set_resetter));

        // Configure the Brinkman penalization objects for the first- and second-order solvers.
        Pointer<BrinkmanPenalizationMethod> fo_brinkman =
            new BrinkmanPenalizationMethod("First-order Brinkman Force",
                                           time_integrator,
                                           app_initializer->getComponentDatabase("FOBrinkmanPenalization"),
                                           /*register_for_restart*/ true);
        Pointer<SOAcousticStreamingBrinkmanPenalization> so_brinkman =
            new SOAcousticStreamingBrinkmanPenalization("Second-order Brinkman Force",
                                                        time_integrator,
                                                        app_initializer->getComponentDatabase("SOBrinkmanPenalization"),
                                                        /*register_for_restart*/ true);
        IBTK::FreeRigidDOFVector cylinder_dofs;
        input_db->getIntegerArray("FREE_DOFS", cylinder_dofs.data(), IBTK::s_max_free_dofs);
        const double cylinder_mass = input_db->getDouble("MASS");
        time_integrator->registerBrinkmanPenalizationStrategy(fo_brinkman,
                                                              so_brinkman,
                                                              phi_var_solid,
                                                              ls_bc_coef,
                                                              solid_level_set_resetter.center,
                                                              cylinder_dofs,
                                                              cylinder_mass);

        // Create level set for the contour integration and register it with acoustic integrator
        const std::string& ls_name_contour = "ARF";
        Pointer<CellVariable<NDIM, double> > phi_var_contour = new CellVariable<NDIM, double>(ls_name_contour);
        ContourLevelSetResetter contour_level_set_resetter;
        contour_level_set_resetter.time_integrator = time_integrator;
        contour_level_set_resetter.center[0] = input_db->getDouble("Contour_XCOM");
        contour_level_set_resetter.center[1] = input_db->getDouble("Contour_YCOM");
        contour_level_set_resetter.R = input_db->getDouble("Contour_Radius");
        time_integrator->registerContourVariable(phi_var_contour, ls_bc_coef, /*contour_value*/ 0.0);
        time_integrator->registerIntegrateHierarchyCallback(reset_contour_level_set_callback_fcn,
                                                            static_cast<void*>(&contour_level_set_resetter));

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Initialize the solid level set function
        solid_level_set_resetter.ls_idx = fo_brinkman->getLevelSetCurrentPatchIndex();
        fcn_ptr(0.0, 0.0, -1, static_cast<void*>(&solid_level_set_resetter));
        solid_level_set_resetter.ls_idx = fo_brinkman->getLevelSetNewPatchIndex();

        // Set the contour level set patch index
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        contour_level_set_resetter.contour_idx =
            var_db->mapVariableAndContextToIndex(phi_var_contour, time_integrator->getCurrentContext());
        reset_contour_level_set_callback_fcn(0.0, 0.0, -1, static_cast<void*>(&contour_level_set_resetter));

        // Remove the AppInitializer
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
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        const auto& arf = time_integrator->getAcousticRadiationForce();
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // Update the center of mass position of the immersed solid and integration contour
            auto& so_vel = time_integrator->getSORigidBodyVelocity();
            auto& X_com_old = time_integrator->getCenterOfMass();
            std::array<double, NDIM> X_com_new;
            int part = 0;
            for (int d = 0; d < NDIM; ++d)
            {
                X_com_new[d] = X_com_old[part][d] + dt * so_vel[part](d);
            }
            solid_level_set_resetter.center = X_com_new;
            contour_level_set_resetter.center = X_com_new;
            time_integrator->updateCenterOfMass(X_com_new, part);

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

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out_file("output");
            out_file << loop_time << '\t' << arf[0][0] << '\t' << arf[0][1] << '\t' << 0.0 << std::endl;
        }

        if (dump_viz_data && uses_visit)
        {
            pout << "Printing velocity and pressure differences as last visit output" << std::endl;
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }

        // Cleanup dumb pointers
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u1_bc_coefs[0][d];
            delete u1_bc_coefs[1][d];
            delete rho_bc_coefs[d];
        }

        delete mu_bc_coef;
        delete lambda_bc_coef;
        delete ls_bc_coef;
    } // cleanup dynamically allocated objects prior to shutdown
} // main
