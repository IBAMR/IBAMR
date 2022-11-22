// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
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
#include <ibamr/GeneralizedIBMethod.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBKirchhoffRodForceGen.h>
#include <ibamr/IBRedundantInitializer.h>
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
        vertex_posn[0] = Point(0.5 + shift, 0.5 + shift, 0.5 + shift);
    }
    return;
}

void
generate_rods_and_directors(const unsigned int& /*struct_num*/,
                            const int& /*ln*/,
                            std::vector<std::vector<double> >& director_spec,
                            std::multimap<int, IBRedundantInitializer::Edge>& /*rod_edge_map*/,
                            std::map<IBRedundantInitializer::Edge,
                                     IBRedundantInitializer::RodSpec,
                                     IBRedundantInitializer::EdgeComp>& /*rod_spec*/,
                            void* /*ctx*/)
{
    director_spec.resize(1);
    director_spec[0].resize(9);
    director_spec[0][0] = director_spec[0][4] = director_spec[0][8] = 1.0;
    return;
}

void
generate_structure_file(const std::string& base_name, int struct_num)
{
    double shift = struct_num == 0 ? -0.25 : 0.25;
    ofstream file;
    file.open(base_name + ".vertex");
    file << "1\n";
    file << 0.5 + shift << " " << 0.5 + shift << " " << 0.5 + shift << "\n";
    file.close();

    file.open(base_name + ".director");
    file << "1\n";
    std::vector<std::vector<double> > director_spec(1);
    director_spec[0].resize(9);
    director_spec[0][0] = director_spec[0][4] = director_spec[0][8] = 1.0;
    for (const auto& directors : director_spec)
    {
        for (int i = 0; i < NDIM; ++i)
            file << directors[3 * i + 0] << " " << directors[3 * i + 1] << " " << directors[3 * i + 2] << "\n";
    }
    file.close();
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
        Pointer<GeneralizedIBMethod> ib_method_ops = new GeneralizedIBMethod(
            "GeneralizedIBMethod", app_initializer->getComponentDatabase("GeneralizedIBMethod"));
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

        const bool use_redundant = input_db->getBool("USE_REDUNDANT");
        // Configure the IB solver.
        Pointer<LInitStrategy> ib_initializer;
        if (use_redundant)
        {
            Pointer<IBRedundantInitializer> initializer = new IBRedundantInitializer(
                "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
            std::vector<std::string> struct_list = { "structure 0", "structure 1" };
            finest_ln = input_db->getInteger("MAX_LEVELS") - 1;
            initializer->setStructureNamesOnLevel(finest_ln, struct_list);
            initializer->registerInitStructureFunction(generate_structure);
            initializer->registerInitDirectorAndRodFunction(generate_rods_and_directors);
            ib_initializer = initializer;
        }
        else
        {
            generate_structure_file("structure_0", 0);
            generate_structure_file("structure_1", 1);
            ib_initializer = new IBStandardInitializer("IBStandardInitializer",
                                                       app_initializer->getComponentDatabase("IBStandardInitializer"));
        }
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBKirchhoffRodForceGen> ib_force_fcn = new IBKirchhoffRodForceGen();
        ib_force_fcn->setUniformBodyForce({ +1.0, -0.25, 0.0 }, 0, finest_ln);
        ib_force_fcn->setUniformBodyForce({ -1.0, 0.0, +0.25 }, 1, finest_ln);
        ib_method_ops->registerIBKirchhoffRodForceGen(ib_force_fcn);

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Main time step loop.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!IBTK::rel_equal_eps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
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
                pout << F[0] << " " << F[1] << " " << F[2] << "\n";
            }
            F_data->restoreArrays();
        }

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
} // main
