// ---------------------------------------------------------------------
//
// Copyright (c) 2021 - 2021 by the IBAMR developers
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
#include <ibtk/config.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CCLaplaceOperator.h>
#include <ibtk/CCPoissonSolverManager.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/PETScAugmentedKrylovLinearSolver.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>

#include "libmesh/boundary_info.h"
#include "libmesh/boundary_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/explicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/node.h"
#include "libmesh/petsc_vector.h"
#include "libmesh/utility.h"

#include <petscsys.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/app_namespaces.h>

// Local includes
#include "CartLaplaceOperator.h"

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
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_poisson.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create variables and register them with the variable database.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");

        Pointer<CellVariable<NDIM, double> > u_cc_var = new CellVariable<NDIM, double>("u_cc");
        Pointer<CellVariable<NDIM, double> > f_cc_var = new CellVariable<NDIM, double>("f_cc");
        Pointer<CellVariable<NDIM, double> > e_cc_var = new CellVariable<NDIM, double>("e_cc");
        Pointer<CellVariable<NDIM, double> > r_cc_var = new CellVariable<NDIM, double>("r_cc");
        Pointer<CellVariable<NDIM, double> > ls_var = new CellVariable<NDIM, double>("ls");

        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVector<NDIM>(2));
        const int f_cc_idx = var_db->registerVariableAndContext(f_cc_var, ctx, IntVector<NDIM>(2));
        const int e_cc_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(2));
        const int r_cc_idx = var_db->registerVariableAndContext(r_cc_var, ctx, IntVector<NDIM>(2));
        const int ls_idx = var_db->registerVariableAndContext(ls_var, ctx, IntVector<NDIM>(2));

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cc_var->getName(), "SCALAR", u_cc_idx);
        visit_data_writer->registerPlotQuantity(f_cc_var->getName(), "SCALAR", f_cc_idx);
        visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "SCALAR", e_cc_idx);
        visit_data_writer->registerPlotQuantity(r_cc_var->getName(), "SCALAR", r_cc_idx);
        visit_data_writer->registerPlotQuantity(ls_var->getName(), "SCALAR", ls_idx);

        // Initialize the AMR patch hierarchy.
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int tag_buffer = 1;
        int level_number = 0;
        bool done = false;
        while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        // Allocate data on each level of the patch hierarchy.
        plog << "Allocating data\n";
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(f_cc_idx, 0.0);
            level->allocatePatchData(e_cc_idx, 0.0);
            level->allocatePatchData(r_cc_idx, 0.0);
            level->allocatePatchData(ls_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int h_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> f_vec("f", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> r_vec("r", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        u_vec.addComponent(u_cc_var, u_cc_idx, h_cc_idx);
        f_vec.addComponent(f_cc_var, f_cc_idx, h_cc_idx);
        e_vec.addComponent(e_cc_var, e_cc_idx, h_cc_idx);
        r_vec.addComponent(r_cc_var, r_cc_idx, h_cc_idx);

        u_vec.setToScalar(0.0);
        f_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);
        r_vec.setToScalar(1.0);

        // Setup exact solutions.
        plog << "Setting up solution data\n";
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);
        muParserCartGridFunction ls_fcn("ls", app_initializer->getComponentDatabase("ls"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(e_cc_idx, e_cc_var, patch_hierarchy, 0.0);
        f_fcn.setDataOnPatchHierarchy(f_cc_idx, f_cc_var, patch_hierarchy, 0.0);
        ls_fcn.setDataOnPatchHierarchy(ls_idx, ls_var, patch_hierarchy, 0.0);

        // Set up the finite element mesh
        // Note we use this to create "augmented" dofs.
        plog << "Creating finite element mesh\n";
        const double R = input_db->getDouble("RADIUS");
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        const std::string elem_type = input_db->getString("ELEM_TYPE");
        Mesh mesh(init.comm(), NDIM);
        const double num_circum_segments = 2.0 * M_PI * R / ds;
        const int r = std::log2(0.25 * num_circum_segments);
        MeshTools::Generation::build_sphere(mesh, R, r, Utility::string_to_enum<ElemType>(elem_type));
        // Ensure the nodes on the surface are on the analytic boundary.
        MeshBase::element_iterator el_end = mesh.elements_end();
        for (MeshBase::element_iterator el = mesh.elements_begin(); el != el_end; ++el)
        {
            Elem* const elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor_ptr(side);
                if (!at_mesh_bdry) continue;
                for (unsigned int k = 0; k < elem->n_nodes(); ++k)
                {
                    if (!elem->is_node_on_side(k, side)) continue;
                    Node& n = elem->node_ref(k);
                    n = R * n.unit();
                }
            }
        }
        mesh.prepare_for_use();
        // Extract boundary mesh
        plog << "Extracting boundary mesh\n";
        BoundaryMesh bdry_mesh(mesh.comm(), mesh.mesh_dimension() - 1);
        BoundaryInfo bdry_info = mesh.get_boundary_info();
        bdry_info.sync(bdry_mesh);
        bdry_mesh.prepare_for_use();

        // Set up dummy systems for drawing
        auto vol_io = libmesh_make_unique<ExodusII_IO>(mesh);
        auto bdry_io = libmesh_make_unique<ExodusII_IO>(bdry_mesh);
        EquationSystems vol_eq_sys(mesh), bdry_eq_sys(bdry_mesh);

        vol_eq_sys.add_system<ExplicitSystem>("X");
        bdry_eq_sys.add_system<ExplicitSystem>("X");
        vol_eq_sys.init();
        bdry_eq_sys.init();

        // Set up the augmented Vec
        // Use libMesh's PetscVector. Create a vector with n_nodes elements.
        // Note we will have to be careful about parallel implementation!
        plog << "Creating PETSc Vector\n";
        PetscVector<double> aug_vec(bdry_mesh.comm(), bdry_mesh.n_nodes());

        plog << "Creating Operator\n";
        CartLaplaceOperator lap_op("LapOp", &bdry_mesh, 2.0 * dx);
        lap_op.setPatchHierarchy(patch_hierarchy);
        lap_op.setLS(ls_idx);
        plog << "Sorting DOFs\n";
        lap_op.sortLagDOFsToCells();
        plog << "Outputting data\n";
        lap_op.printPtMap(pout);

        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
        vol_io->write_timestep("vol_mesh.ex", vol_eq_sys, 1, 0.0);
        bdry_io->write_timestep("bdry_mesh.ex", bdry_eq_sys, 1, 0.0);

    } // cleanup dynamically allocated objects prior to shutdown
} // main
