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
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/StableCentroidPartitioner.h>
#include <ibtk/libmesh_utilities.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// test stuff
#include "../tests.h"

// This test spreads a uniform fluid source from an FE structure that straddles
// the boundaries between patches (IBFEMethod::spreadFluidSource()) and checks
// that the integral of the spread source over the Cartesian grid equals the
// integral of the source over the structure, whatever the layout of the
// patches.

static double uniform_source_strength = 0.0;

void
uniform_source_function(double& Q,
                        const TensorValue<double>& /*FF*/,
                        const libMesh::Point& /*X*/,
                        const libMesh::Point& /*s*/,
                        Elem* /*elem*/,
                        const std::vector<const std::vector<double>*>& /*system_var_data*/,
                        const std::vector<const std::vector<VectorValue<double>>*>& /*system_grad_var_data*/,
                        double /*data_time*/,
                        void* /*ctx*/)
{
    Q = uniform_source_strength;
    return;
} // uniform_source_function

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create a square (or cube) of finite elements in the middle of the domain.
        std::vector<std::unique_ptr<ReplicatedMesh>> meshes;
        meshes.emplace_back(std::make_unique<ReplicatedMesh>(init.comm(), NDIM));
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        const auto elem_type = Utility::string_to_enum<ElemType>(input_db->getString("ELEM_TYPE"));
        const double lo = input_db->getDouble("STRUCTURE_LOWER");
        const double up = input_db->getDouble("STRUCTURE_UPPER");
        const int n_elem_per_side = static_cast<int>((up - lo) / ds);
        ReplicatedMesh& mesh = *meshes[0];
        if (NDIM == 2)
        {
            MeshTools::Generation::build_square(mesh, n_elem_per_side, n_elem_per_side, lo, up, lo, up, elem_type);
        }
        else
        {
            MeshTools::Generation::build_cube(
                mesh, n_elem_per_side, n_elem_per_side, n_elem_per_side, lo, up, lo, up, lo, up, elem_type);
        }
        mesh.prepare_for_use();
        StableCentroidPartitioner partitioner;
        partitioner.partition(mesh);

        // Create major algorithm and data objects that comprise the application.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"), false);
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy =
            new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry, false);
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();

        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"),
            false);

        std::vector<libMesh::MeshBase*> mesh_ptrs = { &mesh };
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           mesh_ptrs,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                           false);
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator,
                                              false);

        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer,
                                        false);

        ib_method_ops->initializeFEEquationSystems();
        uniform_source_strength = input_db->getDouble("SOURCE_STRENGTH");
        ib_method_ops->registerLagBodySourceFunction(uniform_source_function);
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Avoid inconsistent initial partitionings in parallel (see
        // IBFEMethod::d_skip_initial_workload_log).
        if (IBTK_MPI::getNodes() != 1) time_integrator->regridHierarchy();

        // Set up a cell-centered variable, with ghost cells, for the spread source.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<CellVariable<NDIM, double>> q_var = new CellVariable<NDIM, double>("q_test");
        const int n_ghosts = LEInteractor::getMinimumGhostWidth(input_db->getString("IB_DELTA_FUNCTION"));
        const int q_idx = var_db->registerVariableAndContext(q_var, var_db->getContext("q_test"), n_ghosts);
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(q_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<CellData<NDIM, double>> q_data = level->getPatch(p())->getPatchData(q_idx);
                q_data->fillAll(0.0);
            }
        }

        // The test: spread the source and integrate it over the grid.
        const double data_time = 0.0;
        ib_method_ops->computeLagrangianFluidSource(data_time);
        ib_method_ops->spreadFluidSource(q_idx, nullptr, {}, data_time);

        double eulerian_integral = 0.0;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> q_data = patch->getPatchData(q_idx);
                Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
                const double* const patch_dx = patch_geom->getDx();
                double cell_volume = 1.0;
                for (int d = 0; d < NDIM; ++d) cell_volume *= patch_dx[d];
                for (CellIterator<NDIM> ic(patch->getBox()); ic; ic++)
                {
                    eulerian_integral += (*q_data)(ic()) * cell_volume;
                }
            }
        }
        eulerian_integral = IBTK_MPI::sumReduction(eulerian_integral);

        double lagrangian_integral = 0.0;
        for (auto el = mesh.active_elements_begin(); el != mesh.active_elements_end(); ++el)
        {
            lagrangian_integral += uniform_source_strength * (*el)->volume();
        }

        Pointer<PatchLevel<NDIM>> finest_level = patch_hierarchy->getPatchLevel(finest_ln);
        const int num_patches = finest_level->getNumberOfPatches();
        plog << std::setprecision(10);
        plog << "Number of patches on the finest level = " << num_patches << "\n";
        plog << "Lagrangian source integral = " << lagrangian_integral << "\n";
        plog << "Eulerian source integral = " << eulerian_integral << "\n";

        if (std::abs(eulerian_integral - lagrangian_integral) > 1.0e-10 * std::abs(lagrangian_integral))
        {
            TBOX_ERROR("The integral of the spread source (" << eulerian_integral
                                                             << ") is not the integral of the source over the "
                                                                "structure ("
                                                             << lagrangian_integral << ")\n");
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
