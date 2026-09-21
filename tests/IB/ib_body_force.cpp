// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2026 by the IBAMR developers
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

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Set up application namespace declarations
#include <array>

#include <ibamr/app_namespaces.h>

int finest_ln;
void
generate_structure(const unsigned int& struct_num,
                   const int& ln,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn,
                   void* /*ctx*/)
{
    if (ln != finest_ln)
    {
        num_vertices = 0;
        vertex_posn.resize(num_vertices);
    }
    else
    {
        double shift = (struct_num == 0) ? -0.25 : +0.25;
        num_vertices = 1;
        vertex_posn.resize(num_vertices);
        vertex_posn[0] = Point(0.5 + shift, 0.5 + shift);
    }
    return;
}

int
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

#ifndef IBTK_HAVE_SILO
    // Suppress warnings caused by running without silo
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
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
        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        std::vector<std::string> struct_list = { "structure 0", "structure 1" };
        finest_ln = input_db->getInteger("MAX_LEVELS") - 1;
        ib_initializer->setStructureNamesOnLevel(finest_ln, struct_list);
        ib_initializer->registerInitStructureFunction(generate_structure);
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        ib_force_fcn->setUniformBodyForce({ +1.0, -0.25 }, 0, finest_ln);
        ib_force_fcn->setUniformBodyForce({ -1.0, +0.25 }, 1, finest_ln);
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Optionally check that the TimePoint overload of IBMethod::getForceData() returns the actual force data
        // computed by computeLagrangianForce(), rather than position data, and stop before time stepping. Print the
        // same per-vertex force components that the main time-stepping loop below prints, so the numerical values
        // are checked against the expected output.
        const bool check_time_point_data = input_db->getBoolWithDefault("CHECK_TIME_POINT_DATA", false);
        if (check_time_point_data)
        {
            const double current_time = time_integrator->getIntegratorTime();
            const double new_time = current_time + time_integrator->getMaximumTimeStepSize();
            ib_method_ops->preprocessIntegrateData(current_time, new_time, 1);
            const TimePoint time_points[] = { TimePoint::CURRENT_TIME, TimePoint::HALF_TIME };
            for (const auto time_pt : time_points)
            {
                const double time = ib_method_ops->convertTimeEnumToDouble(time_pt);
                ib_method_ops->computeLagrangianForce(time);
                std::vector<Pointer<LData>>* F_data = nullptr;
                bool* F_needs_ghost_fill = nullptr;
                ib_method_ops->getForceData(&F_data, &F_needs_ghost_fill, time_pt);
                auto* F_array = (*F_data)[finest_ln]->getVecArray();
                for (auto&& F : *F_array)
                {
                    pout << F[0] << " " << F[1] << "\n";
                }
                (*F_data)[finest_ln]->restoreArrays();
            }
            ib_method_ops->postprocessIntegrateData(current_time, new_time, 1);
        }

        // Main time step loop.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!check_time_point_data && !IBTK::rel_equal_eps(loop_time, loop_time_end) &&
               time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "At beginning of timestep # " << iteration_num << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;

            auto F_data = ib_method_ops->getLDataManager()->getLData("F", finest_ln);
            auto* F_array = F_data->getVecArray();
            for (auto&& F : *F_array)
            {
                pout << F[0] << " " << F[1] << "\n";
            }
            F_data->restoreArrays();
        }

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
} // main
