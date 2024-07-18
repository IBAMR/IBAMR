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

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/matlab_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_triangle_interface.h>

// Application specific includes
#include <ibamr/FESurfaceDistanceEvaluator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

static double shift_x, shift_y;
#if (NDIM == 3)
static double shift_z;
#endif

double dx, ds;
int coarsest_ln = 0;

// Struct to maintain the properties of the circular interface
struct CircularInterface
{
    Eigen::Vector3d X0;
    double R;
};
CircularInterface circle;

void
calculate_distance_analytically(SAMRAIPointer<PatchHierarchyNd> patch_hierarchy, int E_idx)
{
    int hier_finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= hier_finest_ln; ++ln)
    {
        SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            SAMRAIPointer<PatchNd> patch = level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();
            SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
            const double* patch_X_lower = patch_geom->getXLower();
            const hier::IndexNd& patch_lower_idx = patch_box.lower();
            const double* const patch_dx = patch_geom->getDx();

            SAMRAIPointer<CellDataNd<double> > E_data = patch->getPatchData(E_idx);
            for (BoxNd::Iterator it(patch_box); it; it++)
            {
                CellIndexNd ci(it());

                // Get physical coordinates
                IBTK::Vector coord = IBTK::Vector::Zero();
                for (int d = 0; d < NDIM; ++d)
                {
                    coord[d] = patch_X_lower[d] + patch_dx[d] * (static_cast<double>(ci(d) - patch_lower_idx(d)) + 0.5);
                }
                const double distance =
                    std::sqrt(std::pow((coord[0] - circle.X0(0)), 2.0) + std::pow((coord[1] - circle.X0(1)), 2.0)
#if (NDIM == 3)
                              + std::pow((coord[2] - circle.X0(2)), 2.0)
#endif
                    );

                (*E_data)(ci) = distance - circle.R;
            }
        }
    }

    return;
} // calculate_distance_analytically

void
calculate_error_near_band(SAMRAIPointer<PatchHierarchyNd> patch_hierarchy,
                          int E_idx,
                          int d_idx,
                          int W_idx,
                          double& E_interface,
                          int& num_interface_pts,
                          double& volume_near_interface)

{
    E_interface = 0.0;
    num_interface_pts = 0;
    volume_near_interface = 0.0;
    int hier_finest_ln = patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= hier_finest_ln; ++ln)
    {
        SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevelNd::Iterator p(level); p; p++)
        {
            SAMRAIPointer<PatchNd> patch = level->getPatch(p());
            const BoxNd& patch_box = patch->getBox();
            SAMRAIPointer<CellDataNd<double> > D_data = patch->getPatchData(d_idx);
            SAMRAIPointer<CellDataNd<double> > E_data = patch->getPatchData(E_idx);
            SAMRAIPointer<CellDataNd<double> > W_data = patch->getPatchData(W_idx);
            for (BoxNd::Iterator it(patch_box); it; it++)
            {
                CellIndexNd ci(it());
                const double phi = (*D_data)(ci);
                const double err = (*E_data)(ci);
                const double dV = (*W_data)(ci);
                SAMRAIPointer<CartesianPatchGeometryNd> patch_geom = patch->getPatchGeometry();
                const double* const patch_dx = patch_geom->getDx();

                if (std::abs(err) <= 4.0 * patch_dx[0])
                {
                    if (std::abs(phi) <= 1.2 * patch_dx[0])
                    {
                        E_interface += std::abs(err - phi) * dV;
                        volume_near_interface += dV;
                        num_interface_pts++;
                    }
                    (*E_data)(ci) = std::abs(err - phi);
                }
                else
                {
                    (*E_data)(ci) = 0.0;
                }
            }
        }
    }
    num_interface_pts = IBTK_MPI::sumReduction(num_interface_pts);
    E_interface = IBTK_MPI::sumReduction(E_interface);
    volume_near_interface = IBTK_MPI::sumReduction(volume_near_interface);

    return;
} // calculate_error_near_band

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

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        auto app_initializer = make_samrai_shared<AppInitializer>(argc, argv, "FESurfaceDistance.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
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

        // Setup solid information
        circle.R = input_db->getDouble("R");
        circle.X0[0] = input_db->getDouble("XCOM");
        circle.X0[1] = input_db->getDouble("YCOM");
#if (NDIM == 3)
        circle.X0[2] = input_db->getDouble("ZCOM");
#endif

        // Create a simple FE mesh.
        Mesh solid_mesh(init.comm(), NDIM);

        // Create mesh based upon input file
        dx = input_db->getDouble("DX");
        ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        const bool use_vol_extracted_bdry_mesh = input_db->getBool("USE_VOLUME_EXTRACTED_BOUNDARY_MESH");
        shift_x = circle.X0[0];
        shift_y = circle.X0[1];
#if (NDIM == 3)
        shift_z = circle.X0[2];
#endif
        if (input_db->keyExists("XDA_FILENAME"))
        {
            TBOX_ASSERT(elem_type == "TRI3" || elem_type == "TRI6");

            std::string filename = input_db->getString("XDA_FILENAME");
            MatlabIO distmesh(solid_mesh);
            distmesh.read(filename);

            if (elem_type == "TRI6") solid_mesh.all_second_order();
        }
        else if (input_db->keyExists("GMSH_FILENAME"))
        {
            TBOX_ASSERT(elem_type == "TRI3" || elem_type == "TRI6");

            std::string filename = input_db->getString("GMSH_FILENAME");
            solid_mesh.read(filename);

            if (elem_type == "TRI6") solid_mesh.all_second_order();
        }
        else if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
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
            triangle.elem_type() = Utility::string_to_enum<ElemType>(elem_type);
            triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();
#else
            TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                       << "       but Triangle is required for TRI3 or TRI6 elements.\n");
#endif
        }
        else if (NDIM == 3)
        {
            // NOTE: number of segments along boundary is 4*2^r.
            const double num_circum_segments = 2.0 * M_PI * circle.R / ds;
            const int r = log2(0.25 * num_circum_segments);
            MeshTools::Generation::build_sphere(solid_mesh, circle.R, r, HEX8);
        }
        else
        {
            TBOX_ERROR("ERROR: Invalid element type\n"
                       << "currently this routine supports only TRI3 or EDGE2 at the boundaries");
        }

        // Translate the mesh to a new location.
        // This can also be achieved by libMesh's canned routine MeshTools::Modification::translate().
        for (MeshBase::node_iterator it = solid_mesh.nodes_begin(); it != solid_mesh.nodes_end(); ++it)
        {
            Node* n = *it;
            libMesh::Point& x = *n;
            x(0) += shift_x;
            x(1) += shift_y;
#if (NDIM == 3)
            x(2) += shift_z;
#endif
        }
        solid_mesh.prepare_for_use();

        // Create boundary mesh
        BoundaryMesh boundary_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        BoundaryInfo& boundary_info = solid_mesh.get_boundary_info();
        boundary_info.sync(boundary_mesh);
        boundary_mesh.prepare_for_use();
#if (NDIM == 3)
        // Convert the QUAD4 at the boundaries to TRI3
        MeshTools::Modification::all_tri(boundary_mesh);
#endif

        Mesh& mesh = solid_mesh;

        auto grid_geometry = make_samrai_shared<CartesianGridGeometryNd>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        auto patch_hierarchy = make_samrai_shared<PatchHierarchyNd>("PatchHierarchy", grid_geometry);
        auto error_detector = make_samrai_shared<StandardTagAndInitializeNd>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        auto box_generator = make_samrai_shared<BergerRigoutsosNd>();
        auto load_balancer =
            make_samrai_shared<LoadBalancerNd>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        auto gridding_algorithm =
            make_samrai_shared<GriddingAlgorithmNd>("GriddingAlgorithm",
                                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                                    error_detector,
                                                    box_generator,
                                                    load_balancer);

        // Compute distance function from FE mesh.
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        SAMRAIPointer<CellVariableNd<double> > n_var = make_samrai_shared<CellVariableNd<double> >("num_elements", 1);
        SAMRAIPointer<CellVariableNd<double> > d_var = make_samrai_shared<CellVariableNd<double> >("distance", 1);
        const IntVectorNd no_width = 0;
        SAMRAIPointer<VariableContext> main_ctx = var_db->getContext("Main");
        const int n_idx = var_db->registerVariableAndContext(n_var, main_ctx, no_width);
        const int d_idx = var_db->registerVariableAndContext(d_var, main_ctx, no_width);

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

        int hier_finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= hier_finest_ln; ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(n_idx, 0.0);
            level->allocatePatchData(d_idx, 0.0);
        }
        HierarchyCellDataOpsRealNd<double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, hier_finest_ln);
        hier_cc_data_ops.setToScalar(n_idx, 0.0);
        hier_cc_data_ops.setToScalar(d_idx, 5 * dx);

        // Set up visualization plot file writers.
        SAMRAIPointer<VisItDataWriterNd> visit_data_writer = app_initializer->getVisItDataWriter();
        std::unique_ptr<ExodusII_IO> exodus_io(uses_exodus ? new ExodusII_IO(mesh) : nullptr);

        visit_data_writer->registerPlotQuantity(d_var->getName(), "SCALAR", d_idx);
        const int gcw = input_db->getIntegerWithDefault("GCW", 1);
        FESurfaceDistanceEvaluator surface_distance_eval("FESurfaceDistanceEvaluator",
                                                         patch_hierarchy,
                                                         mesh,
                                                         boundary_mesh,
                                                         /*gcw*/ gcw,
                                                         use_vol_extracted_bdry_mesh);
        pout << "Started mapping intersections" << std::endl;
        surface_distance_eval.mapIntersections();
        pout << "Finished mapping intersections" << std::endl;
        pout << "Computing face normal" << std::endl;
        surface_distance_eval.calculateSurfaceNormals();
        pout << "Finished calculation of face normal" << std::endl;
        pout << "Computing distances" << std::endl;
        surface_distance_eval.computeSignedDistance(n_idx, d_idx);
        pout << "Finished computing distances" << std::endl;
        pout << "Updating sign" << std::endl;
        surface_distance_eval.updateSignAwayFromInterface(d_idx, patch_hierarchy, 5 * dx);
        pout << "Finished updating sign" << std::endl;

        // Compute the error.
        SAMRAIPointer<CellVariableNd<double> > E_var = make_samrai_shared<CellVariableNd<double> >("E");
        const int E_idx = var_db->registerVariableAndContext(E_var, main_ctx);
        for (int ln = 0; ln <= hier_finest_ln; ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            if (!level->checkAllocated(E_idx)) level->allocatePatchData(E_idx, 0.0);
        }

        auto hier_math_ops =
            make_samrai_shared<HierarchyMathOps>("HierarchyMathOps", patch_hierarchy, coarsest_ln, hier_finest_ln);
        const int wgt_cc_idx = hier_math_ops->getCellWeightPatchDescriptorIndex();

        double E_interface = 0.0;
        int num_interface_pts = 0;
        double volume_near_interface = 0.0;
        calculate_distance_analytically(patch_hierarchy, E_idx);
        calculate_error_near_band(
            patch_hierarchy, E_idx, d_idx, wgt_cc_idx, E_interface, num_interface_pts, volume_near_interface);
        pout << "Error in D near interface after level set initialization:" << std::endl
             << "L1-norm:  " << std::setprecision(10) << E_interface / volume_near_interface << std::endl;
        pout << "Number of points within the interface (used to compute interface error):" << std::endl
             << num_interface_pts << std::endl;

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "Number of boundary elements = " << boundary_mesh.n_elem() << std::endl;
            out << "Error in distance near interface: L1-norm = " << std::setprecision(10)
                << E_interface / volume_near_interface << std::endl;
            out << "Number of points used to compute the interface error = " << num_interface_pts << std::endl;
        }
        // Disable plotting by default in the test suite:
        if (false)
        {
            const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();
            if (dump_viz_data)
            {
                pout << "\n\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
                }
                if (uses_exodus)
                {
                    exodus_io->write(exodus_filename);
                }
            }
        }
    }
}
