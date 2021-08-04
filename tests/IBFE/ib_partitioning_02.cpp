// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/IBFESurfaceMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/StableCentroidPartitioner.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/centroid_partitioner.h>
#include <libmesh/equation_systems.h>
#include <libmesh/gmv_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

#include <petscsys.h>
#include <petscvec.h>

#include <boost/multi_array.hpp>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <HierarchyDataOpsManager.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>
#include <mpi.h>

#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

// This file is the main driver for parallel IB element partitioning for surface meshes.

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");

        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create a simple FE mesh.
        Mesh solid_mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        string elem_type = input_db->getString("ELEM_TYPE");
        const double R = 0.3;
        // we want R/2^n_refinements = dx so that we have roughly equal spacing for both elements and cells
        const int n_refinements = int(std::log2(R / dx));
        const std::string shape = input_db->getStringWithDefault("SHAPE", "SPHERE");
        // no support for tetrahedral spheres in 3D yet
        if (shape == "CUBE")
        {
            MeshTools::Generation::build_cube(solid_mesh,
                                              n_refinements,
                                              n_refinements,
                                              n_refinements,
                                              -R,
                                              R,
                                              -R,
                                              R,
                                              -R,
                                              R,
                                              Utility::string_to_enum<ElemType>(elem_type),
                                              false);
        }
        else if (shape == "SQUARE")
        {
            MeshTools::Generation::build_square(solid_mesh,
                                                n_refinements,
                                                n_refinements,
                                                -R,
                                                R,
                                                -R,
                                                R,
                                                Utility::string_to_enum<ElemType>(elem_type),
                                                false);
        }
        else
        {
            TBOX_ASSERT(shape == "SPHERE");
            MeshTools::Generation::build_sphere(
                solid_mesh, R, n_refinements, Utility::string_to_enum<ElemType>(elem_type), 10);
        }
        solid_mesh.prepare_for_use();
        // metis does a good job partitioning, but the partitioning relies on
        // random numbers: the seed changed in libMesh commit
        // 98cede90ca8837688ee13aac5e299a3765f083da (between 1.3.1 and
        // 1.4.0). Hence, to achieve consistent partitioning, use a simpler partitioning scheme:
        IBTK::StableCentroidPartitioner partitioner;
        partitioner.partition(solid_mesh);

        for (auto node = solid_mesh.nodes_begin(); node != solid_mesh.nodes_end(); ++node)
        {
            (**node)(0) += 0.6;
            (**node)(1) += 0.5;
#if NDIM == 3
            (**node)(2) += 0.5;
#endif
        }

        BoundaryMesh boundary_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        BoundaryInfo& boundary_info = solid_mesh.get_boundary_info();
        boundary_info.sync(boundary_mesh);
        boundary_mesh.prepare_for_use();

        bool use_boundary_mesh = input_db->getBoolWithDefault("USE_BOUNDARY_MESH", false);
        Mesh& mesh = use_boundary_mesh ? boundary_mesh : solid_mesh;

        GMVIO gmv_out(mesh);
        gmv_out.write("mesh.gmv");

        plog << "Number of elements: " << mesh.n_active_elem() << std::endl;

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

        Pointer<IBStrategy> ib_ops;
        if (use_boundary_mesh)
            ib_ops = new IBFESurfaceMethod(
                "IBFEMethod",
                app_initializer->getComponentDatabase("IBFEMethod"),
                &mesh,
                app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                /*register_for_restart*/ false);
        else
            ib_ops =
                new IBFEMethod("IBFEMethod",
                               app_initializer->getComponentDatabase("IBFEMethod"),
                               &mesh,
                               app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                               false);
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_ops,
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
        if (use_boundary_mesh)
        {
            Pointer<IBFESurfaceMethod> ibfe_ops = ib_ops;
            ibfe_ops->initializeFEEquationSystems();
            ibfe_ops->initializeFEData();
        }
        else
        {
            Pointer<IBFEMethod> ibfe_ops = ib_ops;
            ibfe_ops->initializeFEEquationSystems();
            ibfe_ops->initializeFEData();
        }
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        time_integrator->regridHierarchy();

        FEDataManager* fe_data_manager = nullptr;
        if (use_boundary_mesh)
        {
            Pointer<IBFESurfaceMethod> ibfe_ops = ib_ops;
            fe_data_manager = ibfe_ops->getFEDataManager();
        }
        else
        {
            Pointer<IBFEMethod> ibfe_ops = ib_ops;
            fe_data_manager = ibfe_ops->getFEDataManager();
        }
        const std::string& displacement_name = IBFEMethod::COORDS_SYSTEM_NAME;
        std::unique_ptr<libMesh::PetscVector<double> > ib_vector =
            fe_data_manager->buildIBGhostedVector(displacement_name);
        Vec petsc_vec = ib_vector->vec();

        const int range_tag = 0;
        const int length_tag = 1;
        const int index_tag = 2;

        const int rank = IBTK_MPI::getRank();
        const int n_nodes = IBTK_MPI::getNodes();

        if (rank == 0)
        {
            for (int r = 0; r < n_nodes; ++r)
            {
                PetscInt bounds[2];
                std::vector<PetscInt> indices;

                if (r == 0)
                {
                    int ierr = VecGetOwnershipRange(petsc_vec, &bounds[0], &bounds[1]);
                    TBOX_ASSERT(ierr == 0);

                    ISLocalToGlobalMapping mapping;
                    ierr = VecGetLocalToGlobalMapping(petsc_vec, &mapping);
                    TBOX_ASSERT(ierr == 0);
                    TBOX_ASSERT(mapping);
                    PetscInt n;
                    ierr = ISLocalToGlobalMappingGetSize(mapping, &n);
                    TBOX_ASSERT(ierr == 0);
                    const PetscInt* petsc_indices;
                    ierr = ISLocalToGlobalMappingGetIndices(mapping, &petsc_indices);
                    TBOX_ASSERT(ierr == 0);

                    indices.insert(indices.begin(), petsc_indices, petsc_indices + n);
                }
                else
                {
                    int ierr =
                        MPI_Recv(bounds, 2, MPIU_INT, r, range_tag, IBTK_MPI::getCommunicator(), MPI_STATUS_IGNORE);
                    TBOX_ASSERT(ierr == 0);

                    PetscInt n = -1;
                    ierr = MPI_Recv(&n, 1, MPIU_INT, r, length_tag, IBTK_MPI::getCommunicator(), MPI_STATUS_IGNORE);
                    TBOX_ASSERT(ierr == 0);
                    indices.resize(n);
                    ierr = MPI_Recv(
                        indices.data(), n, MPIU_INT, r, index_tag, IBTK_MPI::getCommunicator(), MPI_STATUS_IGNORE);
                    TBOX_ASSERT(ierr == 0);
                }

                plog << "locally owned range on " << r << ": [" << bounds[0] << ", " << bounds[1] << ")\n";
                plog << "total length on " << r << ": " << indices.size() << '\n';
                plog << "global index ranges on " << r << ":\n";
                // Print out intervals instead of every single number
                for (PetscInt i = 0; i < long(indices.size()); ++i)
                {
                    if (i == 0)
                    {
                        plog << '[' << indices[i] << ", ";
                    }
                    else
                    {
                        // if we started a new range of consecutive integers
                        // then start a new interval on a new line
                        if (indices[i] != indices[i - 1] + 1)
                        {
                            plog << indices[i - 1] + 1 << ')' << '\n';
                            plog << '[' << indices[i] << ", ";
                        }
                    }
                }
                plog << indices.back() + 1 << ")\n";
            }
        }
        else
        {
            // send local range
            PetscInt bounds[2];
            int ierr = VecGetOwnershipRange(petsc_vec, &bounds[0], &bounds[1]);
            TBOX_ASSERT(ierr == 0);
            ierr = MPI_Send(bounds, 2, MPIU_INT, 0, range_tag, IBTK_MPI::getCommunicator());
            TBOX_ASSERT(ierr == 0);

            // send locally stored global indices
            ISLocalToGlobalMapping mapping;
            ierr = VecGetLocalToGlobalMapping(petsc_vec, &mapping);
            TBOX_ASSERT(ierr == 0);
            TBOX_ASSERT(mapping);
            PetscInt n;
            ierr = ISLocalToGlobalMappingGetSize(mapping, &n);
            TBOX_ASSERT(ierr == 0);
            const PetscInt* indices;
            ierr = ISLocalToGlobalMappingGetIndices(mapping, &indices);
            TBOX_ASSERT(ierr == 0);

            ierr = MPI_Send(&n, 1, MPIU_INT, 0, length_tag, IBTK_MPI::getCommunicator());
            TBOX_ASSERT(ierr == 0);

            ierr = MPI_Send(indices, n, MPIU_INT, 0, index_tag, IBTK_MPI::getCommunicator());
            TBOX_ASSERT(ierr == 0);
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
