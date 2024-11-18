// ---------------------------------------------------------------------
//
// Copyright (c) 2022 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Check the 'interpolate cell/side/face to nodal' routines.

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAIVectorReal.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

#include "../tests.h"

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

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "nodalinterp.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        const std::string var_type = input_db->getStringWithDefault("var_type", "SIDE");
        const auto N = input_db->getInteger("N");

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
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

        Pointer<hier::Variable<NDIM> > u_var, e_var;
        Pointer<CellVariable<NDIM, double> > u_cc_var = new CellVariable<NDIM, double>("u_cc", NDIM);
        Pointer<SideVariable<NDIM, double> > u_sc_var = new SideVariable<NDIM, double>("u_sc");
        Pointer<FaceVariable<NDIM, double> > u_fc_var = new FaceVariable<NDIM, double>("u_fc");

        if (var_type == "CELL")
        {
            u_var = u_cc_var;
        }
        else if (var_type == "SIDE")
        {
            u_var = u_sc_var;
        }
        else if (var_type == "FACE")
        {
            u_var = u_fc_var;
        }
        else
        {
            TBOX_ERROR("not implemented");
        }

        const int u_idx = var_db->registerVariableAndContext(u_var, ctx, IntVector<NDIM>(1));

        const bool fine_boundary_represents_var = input_db->getBoolWithDefault("fine_boundary_represents_var", false);
        Pointer<NodeVariable<NDIM, double> > u_nc_var =
            new NodeVariable<NDIM, double>("u_nc", NDIM, fine_boundary_represents_var);
        Pointer<NodeVariable<NDIM, double> > e_nc_var =
            new NodeVariable<NDIM, double>("e_nc", NDIM, fine_boundary_represents_var);
        Pointer<CellVariable<NDIM, double> > e_cc_var = new CellVariable<NDIM, double>("e_cc", NDIM);

        // Don't add more ghosts unless we need them to interpolate_back
        const bool interp_back = input_db->getBoolWithDefault("interp_back", false);
        const int u_nc_idx = var_db->registerVariableAndContext(u_nc_var, ctx, IntVector<NDIM>(interp_back ? 1 : 0));

        int e_idx = IBTK::invalid_index;
        if (interp_back)
        {
            // only one implemented right now in HierarchyMathOps
            TBOX_ASSERT(var_type == "CELL");
            e_idx = var_db->registerVariableAndContext(e_cc_var, ctx, IntVector<NDIM>(0));
            e_var = e_cc_var;
        }
        else
        {
            e_idx = var_db->registerVariableAndContext(e_nc_var, ctx, IntVector<NDIM>(0));
            e_var = e_nc_var;
        }

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_nc_var->getName(), "VECTOR", u_nc_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(u_nc_var->getName() + std::to_string(d), "SCALAR", u_nc_idx, d);
        }

        visit_data_writer->registerPlotQuantity(e_cc_var->getName(), "VECTOR", e_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(e_var->getName() + std::to_string(d), "SCALAR", e_idx, d);
        }

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
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_idx, 0.0);
            level->allocatePatchData(u_nc_idx, 0.0);
            level->allocatePatchData(e_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);

        SAMRAIVectorReal<NDIM, double> u_vec("u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        // if we are interpolating back then we want to use the variable we
        // started with - otherwise compare to MMS with node data.
        if (interp_back)
        {
            u_vec.addComponent(u_var, u_idx);
        }
        else
        {
            u_vec.addComponent(u_nc_var, u_nc_idx);
        }

        e_vec.addComponent(e_var, e_idx);
        u_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        u_fcn.setDataOnPatchHierarchy(u_idx, u_var, patch_hierarchy, 0.0);
        u_fcn.setDataOnPatchHierarchy(e_idx, e_var, patch_hierarchy, 0.0);

        // We need updated ghost values to interpolate from side or face data:
        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> output_transaction_comps;
        output_transaction_comps.emplace_back(
            u_idx, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false);
        HierarchyGhostCellInterpolation hier_bdry_fill;
        hier_bdry_fill.initializeOperatorState(output_transaction_comps, patch_hierarchy);
        hier_bdry_fill.fillData(0.0);

        // interpolate from side/face-centered to nodal:
        const bool synch_dst_cf_interface = input_db->getBoolWithDefault("synch_dst_cf_interface", false);
        const bool synch_src_cf_interface = true;
        if (var_type == "CELL")
        {
            hier_math_ops.interp(u_nc_idx, u_nc_var, synch_dst_cf_interface, u_idx, u_cc_var, nullptr, 0.0);
        }
        else if (var_type == "SIDE")
        {
            hier_math_ops.interp(
                u_nc_idx, u_nc_var, synch_dst_cf_interface, u_idx, u_sc_var, nullptr, 0.0, synch_src_cf_interface);
        }
        else if (var_type == "FACE")
        {
            hier_math_ops.interp(
                u_nc_idx, u_nc_var, synch_dst_cf_interface, u_idx, u_fc_var, nullptr, 0.0, synch_src_cf_interface);
        }

        if (interp_back)
        {
            // only node-to-cell is implemented and it requires updated ghost values
            std::vector<ITC> comps_2;
            comps_2.emplace_back(u_nc_idx, "LINEAR_REFINE", false, "CONSTANT_COARSEN", "LINEAR", false);
            HierarchyGhostCellInterpolation hier_bdry_fill_2;
            hier_bdry_fill_2.initializeOperatorState(comps_2, patch_hierarchy);
            hier_bdry_fill_2.fillData(0.0);

            hier_math_ops.interp(e_idx, e_cc_var, u_nc_idx, u_nc_var, nullptr, 0.0, true);
        }
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&e_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double> >(&u_vec, false));

        // Compute error and print error norms.
        const double max_norm = e_vec.maxNorm();

        std::ostringstream out;
        out << std::setprecision(16);
        if (IBTK_MPI::getRank() == 0)
        {
            out << "|e|_oo = " << max_norm << "\n";
        }

        // Set invalid values on coarse levels (i.e., coarse-grid values that
        // are covered by finer grid patches) to equal zero.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber() - 1; ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            BoxArray<NDIM> refined_region_boxes;
            Pointer<PatchLevel<NDIM> > next_finer_level = patch_hierarchy->getPatchLevel(ln + 1);
            refined_region_boxes = next_finer_level->getBoxes();
            refined_region_boxes.coarsen(next_finer_level->getRatioToCoarserLevel());
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<NodeData<NDIM, double> > e_nc_data = patch->getPatchData(e_idx);
                if (e_nc_data)
                {
                    for (int i = 0; i < refined_region_boxes.getNumberOfBoxes(); ++i)
                    {
                        const Box<NDIM> refined_box = refined_region_boxes[i];
                        const Box<NDIM> intersection = Box<NDIM>::grow(patch_box, 1) * refined_box;
                        if (!intersection.empty())
                        {
                            e_nc_data->fillAll(0.0, intersection);
                        }
                    }
                }
            }
        }

        if (input_db->getBoolWithDefault("print_cf_plane", false))
        {
            const int d = input_db->getIntegerWithDefault("cf_plane_axis", 1);
            TBOX_ASSERT(0 <= d && d <= NDIM);
            const int depth = input_db->getIntegerWithDefault("component_to_print", 0);
            for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
            {
                out << "level = " << ln << '\n';
                Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM> > patch = level->getPatch(p());
                    Pointer<NodeData<NDIM, double> > u_data = patch->getPatchData(u_nc_idx);
                    Pointer<PatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();

                    for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
                    {
                        const NodeIndex<NDIM>& idx = ni();
                        IBTK::VectorNd point;
                        point[0] = double(idx[0]) / (N * pgeom->getRatio()[0]);
                        point[1] = double(idx[1]) / (N * pgeom->getRatio()[1]);
#if NDIM == 3
                        point[2] = double(idx[2]) / (N * pgeom->getRatio()[2]);
#endif
                        const double value = (*u_data)(idx, depth);

#if NDIM == 2
                        if (std::abs(point[d] - 0.5) < 1e-6)
                        {
                            out << point[0] << ", " << point[1] << ": " << value << '\n';
                        }
#else
                        if (std::abs(point[d] - 0.5) < 1e-6)
                        {
                            out << point[0] << ", " << point[1] << ", " << point[2] << ": " << value << '\n';
                        }
#endif
                    }
                }
            }
        }

        print_strings_on_plog_0(out.str());

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

    } // cleanup dynamically allocated objects prior to shutdown
} // run_example
