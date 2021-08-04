// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2021 by the IBAMR developers
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
#include <libmesh/exodusII_io.h>
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
#include <ibtk/StableCentroidPartitioner.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <boost/multi_array.hpp>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

// test stuff
#include "../tests.h"

// Like spread_02, but explicitly check the values spread from a constant field
// when the structure lives on multiple levels. This shows that the correct value is spread away from
IBTK::Point
exact_forcing(const IBTK::Point&)
{
    IBTK::Point result;
    result[0] = 1.0;
    result[1] = 1.0;
#if NDIM == 3
    result[2] = 1.0;
#endif
    return result;
}

int
main(int argc, char** argv)
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    // suppress warnings caused by using a refinement ratio of 4 and not
    // setting up coarsening correctly
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);

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

        // This test only supports one type of mesh.
        {
            ReplicatedMesh& mesh = *meshes[0];
            const double L = input_db->getDouble("L");
            if (NDIM == 2)
                MeshTools::Generation::build_square(
                    mesh, int(L / ds), int(L / (4.0 * ds)), 0.0, L, 0.0, L / 4.0, elem_type);
            else
                MeshTools::Generation::build_cube(mesh,
                                                  int(L / ds),
                                                  int(L / (4.0 * ds)),
                                                  int(L / (4.0 * ds)),
                                                  0.0,
                                                  L,
                                                  0.0,
                                                  L / 4.0,
                                                  0.0,
                                                  L / 4.0,
                                                  elem_type);

            // Assign boundary elements to a finer patch level
            MeshBase::element_iterator el_end = mesh.active_elements_end();
            for (MeshBase::element_iterator el = mesh.active_elements_begin(); el != el_end; ++el)
            {
                const auto centroid = (*el)->centroid();
                if (centroid(0) < 0.5)
                {
                    (*el)->subdomain_id() = 2;
                }
                else
                {
                    (*el)->subdomain_id() = 1;
                }
            }
        }

        // metis does a good job partitioning, but the partitioning relies on
        // random numbers: the seed changed in libMesh commit
        // 98cede90ca8837688ee13aac5e299a3765f083da (between 1.3.1 and
        // 1.4.0). Hence, to achieve consistent partitioning, use a simpler partitioning scheme:
        for (const auto& mesh : meshes)
        {
            mesh->prepare_for_use();
            StableCentroidPartitioner partitioner;
            partitioner.partition(*mesh);
        }

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
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        time_integrator->registerVisItDataWriter(visit_data_writer);

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

        ib_method_ops->initializeFEEquationSystems();
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
            const std::string& F_name = IBFEMethod::FORCE_SYSTEM_NAME;
            const std::string& X_name = IBFEMethod::COORDS_SYSTEM_NAME;
            EquationSystems* equation_systems = ib_method_ops->getFEDataManager(part_n)->getEquationSystems();
            System& F_system = equation_systems->get_system(F_name);
            System& X_system = equation_systems->get_system(X_name);

            const unsigned int n_vars = F_system.n_vars();
            TBOX_ASSERT(n_vars == NDIM);
            const DofMap& F_dof_map = F_system.get_dof_map();
            const DofMap& X_dof_map = X_system.get_dof_map();

            PetscVector<double>& current_F = dynamic_cast<libMesh::PetscVector<double>&>(F_system.get_vector("half"));
            PetscVector<double>& current_X = *dynamic_cast<PetscVector<double>*>(X_system.current_local_solution.get());

            std::vector<dof_id_type> X_idxs;
            std::vector<dof_id_type> F_idxs;
            // Unlike elements, nodes are uniquely allocated to different processors
            const ReplicatedMesh& mesh = *meshes[part_n];
            for (auto node_iter = mesh.local_nodes_begin(); node_iter != mesh.local_nodes_end(); ++node_iter)
            {
                IBTK::Point X;
                for (unsigned int var_n = 0; var_n < n_vars; ++var_n)
                {
                    IBTK::get_nodal_dof_indices(X_dof_map, *node_iter, var_n, X_idxs);
                    X[var_n] = current_X(X_idxs[0]);
                }

                const auto F = exact_forcing(X);
                for (unsigned int var_n = 0; var_n < n_vars; ++var_n)
                {
                    IBTK::get_nodal_dof_indices(F_dof_map, *node_iter, var_n, F_idxs);
                    current_F.set(F_idxs[0], F[var_n]);
                }
            }
            current_F.close();
            *F_system.solution = current_F;
        }

        std::ostringstream out;
        if (IBTK_MPI::getNodes() != 1)
        {
            // partitioning is only relevant when there are multiple processors
            Pointer<PatchLevel<NDIM> > patch_level =
                patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
            const BoxArray<NDIM> boxes = patch_level->getBoxes();
            plog << "hierarchy boxes:\n";
            for (int i = 0; i < boxes.size(); ++i) plog << boxes[i] << '\n';
            // rank is only relevant when there are multiple processors
            out << "\nrank: " << IBTK_MPI::getRank() << '\n';
        }

        // Here is the actual test:

        const double data_time = time_integrator->getIntegratorTime() + dt / 2;
        RobinPhysBdryPatchStrategy* bdry_op = time_integrator->getVelocityPhysBdryOp();
        ib_method_ops->spreadForce(f_ghost_idx, bdry_op, {}, data_time);

        // velocity interpolation requires that ghost data be present
        ghost_fill_op.fillData(/*time*/ 0.0);

        // Now put it back on the structure:
        ib_method_ops->interpolateVelocity(f_ghost_idx, {}, {}, data_time);

        // The rest is just bookkeeping and some visualization.

        // convert SC to CC for plotting purposes:
        Pointer<CellVariable<NDIM, double> > f_cc_var = new CellVariable<NDIM, double>("f_cc", NDIM);
        const int f_ghost_cc_idx = var_db->registerVariableAndContext(f_cc_var, f_ghost_ctx, IntVector<NDIM>(0));
        visit_data_writer->registerPlotQuantity("f_ghost", "VECTOR", f_ghost_cc_idx);

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity("f_ghost_" + std::to_string(d), "SCALAR", f_ghost_cc_idx, d);
        }

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        // This is what we get for SAMRAI inventing its own type system and then
        // only half-implementing it
        Pointer<SideVariable<NDIM, double> > temp = f_var;
        TBOX_ASSERT(temp);
        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(f_ghost_cc_idx, 0.0);
        }

        hier_math_ops.interp(f_ghost_cc_idx, f_cc_var, f_ghost_idx, temp, nullptr, 0.0, /*synch_cf_interface=*/true);

        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

#ifdef LIBMESH_HAVE_EXODUS_API
        {
            std::unique_ptr<ExodusII_IO> exodus_io(new ExodusII_IO(*meshes[0]));
            EquationSystems* equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();
            exodus_io->write_timestep("out.ex2", *equation_systems, 1, 0.0);
        }
#endif

        {
            const int ln = patch_hierarchy->getFinestLevelNumber();
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                bool printed_value = false;
                std::ostringstream patch_out;
                patch_out << "patch number " << p() << '\n';
                patch_out.precision(16);
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                Pointer<SideData<NDIM, double> > f_data = patch->getPatchData(f_ghost_idx);
                const Box<NDIM> patch_box = patch->getBox();

                // same as SideData::print, but elides values close to 1.
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
                            if (std::abs(value - 1.0) > 1e-2)
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
