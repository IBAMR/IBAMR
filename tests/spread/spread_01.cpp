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

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// other samrai stuff
#include <HierarchyDataOpsManager.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/linear_partitioner.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>
#include <libmesh/periodic_boundary.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// test stuff
#include "../tests.h"

// This file is the main driver for force spreading tests (i.e.,
// IBFEmethod::spreadForce). At the moment it simply prints out the force
// values.

// Coordinate mapping function.
void
coordinate_mapping_function(libMesh::Point& X, const libMesh::Point& s, void* /*ctx*/)
{
    // We have to be careful with how we pick these offsets since these tests
    // validate spreading by comparing pointwise values. In particular: with
    // an odd-order quadrature rule and a quad or hex mesh the cell midpoint
    // is a quadrature point. Due to roundoff this could work out to be, e.g.,
    // 0.5 +/- eps: this is problematic since we will get a different Eulerian
    // cell now based on roundoff accumulation.
    //
    // Get around these rounding issues by picking strange-looking shifts that
    // guarantee no quadrature point will be near 0.5, 0.75, or any other
    // possible Eulerian cell midpoint coordinate.
    X(0) = s(0) + 0.612345;
    X(1) = s(1) + 0.512345;
#if (NDIM == 3)
    X(2) = s(2) + 0.512345;
#endif
    return;
} // coordinate_mapping_function

static constexpr double R = 0.2;
static constexpr double w = 0.0625;

// alternative mapping function for periodic domains.
void
periodic_coordinate_mapping_function(libMesh::Point& X, const libMesh::Point& s, void* /*ctx*/)
{
    static constexpr double gamma = 0.15;
    X(0) = (R + s(1)) * cos(s(0) / R) + 0.505;
    X(1) = (R + gamma + s(1)) * sin(s(0) / R) + 0.505;
    return;
} // periodic_coordinate_mapping_function

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    PetscOptionsSetValue(nullptr, "-ksp_rtol", "1e-16");
    PetscOptionsSetValue(nullptr, "-ksp_atol", "1e-16");

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");

        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create a simple FE mesh.
        std::vector<std::unique_ptr<ReplicatedMesh> > meshes;
        meshes.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        const std::string elem_str = input_db->getString("ELEM_TYPE");
        const auto elem_type = Utility::string_to_enum<ElemType>(elem_str);

        const std::string geometry = input_db->getStringWithDefault("geometry", "sphere");

        if (geometry == "sphere")
        {
            ReplicatedMesh& mesh = *meshes[0];
            // libMesh circa version 1.5 fixed a bug with the way they read
            // Triangle input which results in vertices being numbered in a
            // different way after the patch. This actually makes a
            // substantial difference for this test since changing the vertex
            // numbering changes the quadrature points, which in turn changes
            // where we spread forces. Hence, unlike the examples, just avoid
            // Triangle altogether here.
            //
            // NOTE: number of segments along boundary is 4*2^r.
            const double num_circum_segments = 2.0 * M_PI * R / ds;
            const int r = log2(0.25 * num_circum_segments);
            MeshTools::Generation::build_sphere(mesh, R, r, elem_type);
        }
        else if (geometry == "cube")
        {
            ReplicatedMesh& mesh = *meshes[0];
            const double L = input_db->getDouble("L");
            if (NDIM == 2)
                MeshTools::Generation::build_square(mesh, 10, 12, 0.0, L, 0.0, L, elem_type);
            else
                MeshTools::Generation::build_cube(mesh, 10, 12, 14, 0.0, L, 0.0, L, 0.0, L, elem_type);
        }
        else if (geometry == "periodic")
        {
            const int n_x = ceil(2.0 * M_PI * R / ds);
            const int n_y = ceil(w / ds);
            MeshTools::Generation::build_square(*meshes[0], n_x, n_y, 0.0, 2.0 * M_PI * R, 0.0, w, elem_type);
        }
        else
        {
            TBOX_ASSERT(geometry == "composite_cube");
            const double L = input_db->getDouble("L");
            if (NDIM == 2)
            {
                ReplicatedMesh& mesh_0 = *meshes[0];
                meshes.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
                ReplicatedMesh& mesh_1 = *meshes[1];
                meshes.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
                ReplicatedMesh& mesh_2 = *meshes[2];
                meshes.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
                ReplicatedMesh& mesh_3 = *meshes[3];
                MeshTools::Generation::build_square(mesh_0, 10, 12, 0.0, L / 2, 0.0, L / 2, elem_type);
                MeshTools::Generation::build_square(mesh_1, 10, 12, L / 2, L, 0.0, L / 2, elem_type);
                MeshTools::Generation::build_square(mesh_2, 10, 12, 0.0, L / 2, L / 2, L, elem_type);
                MeshTools::Generation::build_square(mesh_3, 10, 12, L / 2, L, L / 2, L, elem_type);
            }
            else
            {
                ReplicatedMesh& mesh_0 = *meshes[0];
                meshes.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
                ReplicatedMesh& mesh_1 = *meshes[1];
                meshes.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
                ReplicatedMesh& mesh_2 = *meshes[2];
                meshes.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
                ReplicatedMesh& mesh_3 = *meshes[3];
                MeshTools::Generation::build_cube(mesh_0, 5, 6, 7, 0.0, L / 2, 0.0, L / 2, 0.0, L / 2, elem_type);
                MeshTools::Generation::build_cube(mesh_1, 5, 6, 7, L / 2, L, 0.0, L / 2, 0.0, L / 2, elem_type);
                MeshTools::Generation::build_cube(mesh_2, 5, 6, 7, 0.0, L / 2, L / 2, L, 0.0, L / 2, elem_type);
                MeshTools::Generation::build_cube(mesh_3, 5, 6, 7, L / 2, L, L / 2, L, 0.0, L / 2, elem_type);

                meshes.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
                ReplicatedMesh& mesh_4 = *meshes[4];
                meshes.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
                ReplicatedMesh& mesh_5 = *meshes[5];
                meshes.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
                ReplicatedMesh& mesh_6 = *meshes[6];
                meshes.emplace_back(new ReplicatedMesh(init.comm(), NDIM));
                ReplicatedMesh& mesh_7 = *meshes[7];
                MeshTools::Generation::build_cube(mesh_4, 5, 6, 7, 0.0, L / 2, 0.0, L / 2, L / 2, L, elem_type);
                MeshTools::Generation::build_cube(mesh_5, 5, 6, 7, L / 2, L, 0.0, L / 2, L / 2, L, elem_type);
                MeshTools::Generation::build_cube(mesh_6, 5, 6, 7, 0.0, L / 2, L / 2, L, L / 2, L, elem_type);
                MeshTools::Generation::build_cube(mesh_7, 5, 6, 7, L / 2, L, L / 2, L, L / 2, L, elem_type);
            }
        }

        // metis does a good job partitioning, but the partitioning relies on
        // random numbers: the seed changed in libMesh commit
        // 98cede90ca8837688ee13aac5e299a3765f083da (between 1.3.1 and
        // 1.4.0). Hence, to achieve consistent partitioning, use a simpler partitioning scheme:
        for (const auto& mesh : meshes)
        {
            mesh->prepare_for_use();
            LinearPartitioner partitioner;
            partitioner.partition(*mesh);
        }

        VectorValue<double> boundary_translation(2.0 * M_PI * R, 0.0, 0.0);
        PeriodicBoundary pbc(boundary_translation);
        pbc.myboundary = 3;
        pbc.pairedboundary = 1;

        std::size_t n_elem = 0;
        for (const auto& mesh : meshes) n_elem += mesh->n_active_elem();
        plog << "Number of elements: " << n_elem << std::endl;

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"), false);
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy =
            new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry, false);
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();

        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"),
                false);
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"),
                false);
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }

        std::vector<libMesh::MeshBase*> mesh_ptrs;
        for (auto& mesh : meshes)
        {
            mesh_ptrs.emplace_back(mesh.get());
        }
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

        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer,
                                        false);

        // Configure the IBFE solver.
        if (geometry == "sphere")
        {
            ib_method_ops->registerInitialCoordinateMappingFunction(coordinate_mapping_function);
        }
        else if (geometry == "periodic")
        {
            ib_method_ops->registerInitialCoordinateMappingFunction(periodic_coordinate_mapping_function);
        }
        else
        {
            // other meshes are already centered correctly
        }
        ib_method_ops->initializeFEEquationSystems();
        if (geometry == "periodic")
        {
            EquationSystems* equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();
            for (unsigned int k = 0; k < equation_systems->n_systems(); ++k)
            {
                System& system = equation_systems->get_system(k);
                system.get_dof_map().add_periodic_boundary(pbc);
            }
        }
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // only regrid if we are in parallel - do this to avoid weird problems
        // with inconsistent initial partitionings (see
        // IBFEMethod::d_skip_initial_workload_log)
        if (IBTK_MPI::getNodes() != 1) time_integrator->regridHierarchy();

        // Now for the actual test. Set up a new variable containing ghost data:
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        const Pointer<SAMRAI::hier::Variable<NDIM> > f_var = time_integrator->getBodyForceVariable();
        const Pointer<VariableContext> f_ghost_ctx = var_db->getContext("f_ghost");

        int n_ghosts = input_db->keyExists("IB_DELTA_FUNCTION") ?
                           LEInteractor::getMinimumGhostWidth(input_db->getString("IB_DELTA_FUNCTION")) :
                           3;
        if (app_initializer->getComponentDatabase("IBFEMethod")->keyExists("min_ghost_cell_width"))
        {
            n_ghosts = app_initializer->getComponentDatabase("IBFEMethod")->getInteger("min_ghost_cell_width");
        }
        const int f_ghost_idx = var_db->registerVariableAndContext(f_var, f_ghost_ctx, n_ghosts);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(f_ghost_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_ghost_idx);
                f_data->fillAll(0.0);
            }
        }

        // synch ghost data
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_cell_components(1);
        ghost_cell_components[0] = InterpolationTransactionComponent(f_ghost_idx,
                                                                     "CONSERVATIVE_LINEAR_REFINE",
                                                                     true,
                                                                     "CONSERVATIVE_COARSEN",
                                                                     "LINEAR",
                                                                     false,
                                                                     {}, // f_bc_coefs
                                                                     nullptr);
        HierarchyGhostCellInterpolation ghost_fill_op;
        ghost_fill_op.initializeOperatorState(ghost_cell_components, patch_hierarchy);
        ghost_fill_op.fillData(/*time*/ 0.0);

        const double dt = time_integrator->getMaximumTimeStepSize();
        time_integrator->preprocessIntegrateHierarchy(
            time_integrator->getIntegratorTime(), time_integrator->getIntegratorTime() + dt, 1 /*???*/);
        for (unsigned int part_n = 0; part_n < meshes.size(); ++part_n)
        {
            auto& fe_data_manager = *ib_method_ops->getFEDataManager(part_n);
            auto& equation_systems = *fe_data_manager.getEquationSystems();
            auto& force_system = equation_systems.get_system(ib_method_ops->getForceSystemName());
            auto& half_f_vector = dynamic_cast<libMesh::PetscVector<double>&>(force_system.get_vector("half"));
            for (unsigned int i = half_f_vector.first_local_index(); i < half_f_vector.last_local_index(); ++i)
            {
                half_f_vector.set(i, i % 10);
            }
            half_f_vector.close();
        }

        // the partitioning isn't relevant when we only have one processor
        if (IBTK_MPI::getNodes() != 1)
        {
            print_partitioning_on_plog_0(
                patch_hierarchy, patch_hierarchy->getFinestLevelNumber(), patch_hierarchy->getFinestLevelNumber());
        }

        // Test the accumulation code when we have meshes that are against the
        // boundary (i.e., the cube and composite cube geometries)
        const double data_time = time_integrator->getIntegratorTime() + dt / 2;
        RobinPhysBdryPatchStrategy* bdry_op = geometry == "sphere" ? nullptr : time_integrator->getVelocityPhysBdryOp();
        if (input_db->getBoolWithDefault("spread_with_ibfemethod", true))
            ib_method_ops->spreadForce(f_ghost_idx, bdry_op, {}, data_time);
        else
        {
            TBOX_ASSERT(meshes.size() == 1); // not implemented for multiple parts
            auto& fe_data_manager = *ib_method_ops->getFEDataManager(0);
            auto& equation_systems = *fe_data_manager.getEquationSystems();
            auto& force_system = equation_systems.get_system(ib_method_ops->getForceSystemName());
            auto& F_vec = dynamic_cast<libMesh::PetscVector<double>&>(force_system.get_vector("half"));
            auto& position_system = equation_systems.get_system(ib_method_ops->getCurrentCoordinatesSystemName());
            auto& X_vec = dynamic_cast<libMesh::PetscVector<double>&>(*position_system.current_local_solution);
            fe_data_manager.spread(
                f_ghost_idx, F_vec, X_vec, ib_method_ops->getForceSystemName(), nullptr, data_time, false, false);

            // here's the real test for the bug in this case: FEDataManager::spread()
            // previously did not copy values spread into ghost regions, so the values
            // computed by these two combined calls would be wrong.
            bdry_op->setPatchDataIndex(f_ghost_idx);
            const int ln = patch_hierarchy->getFinestLevelNumber();
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                const Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<PatchData<NDIM> > f_data = patch->getPatchData(f_ghost_idx);
                bdry_op->accumulateFromPhysicalBoundaryData(*patch, data_time, f_data->getGhostCellWidth());
            }
        }
        const double cutoff = input_db->getDoubleWithDefault("output_cutoff_value", 0.0);
        std::ostringstream out;
        {
            const int ln = patch_hierarchy->getFinestLevelNumber();
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);

            // We don't need to print this if we are running in serial
            if (IBTK_MPI::getNodes() != 1)
            {
                out << "\nrank: " << IBTK_MPI::getRank() << '\n';
            }
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                bool printed_value = false;
                std::ostringstream patch_out;
                patch_out << "patch number " << p() << '\n';
                patch_out.precision(16);
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_ghost_idx);
                const Box<NDIM> patch_box = patch->getBox();

                // same as SideData::print, but elides zero values. We don't
                // print any information about the patch when no values are
                // above the cutoff.
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    patch_out << "Array side normal = " << axis << std::endl;
                    for (int d = 0; d < f_data->getDepth(); ++d)
                    {
                        patch_out << "Array depth = " << d << std::endl;
                        const ArrayData<NDIM, double>& data = f_data->getArrayData(axis);
                        for (SideIterator<NDIM> i(patch_box, axis); i; i++)
                        {
                            const double value = data(i(), d);
                            if (std::abs(value) > cutoff)
                            {
                                patch_out << "array" << i() << " = " << value << '\n';
                                printed_value = true;
                            }
                        }
                    }
                }
                if (printed_value) out << patch_out.str();
            }
        }
        IBTK_MPI::barrier();

        print_strings_on_plog_0(out.str());
    } // cleanup dynamically allocated objects prior to shutdown
} // main
