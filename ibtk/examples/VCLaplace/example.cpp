// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2022 by the IBAMR developers
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
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for major SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <HierarchyDataOpsManager.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/ibtk_enums.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

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
        auto app_initializer = make_samrai_shared<AppInitializer>(argc, argv, "vc_laplace.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
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

        // Create variables and register them with the variable database.
        VariableDatabaseNd* var_db = VariableDatabaseNd::getDatabase();
        SAMRAIPointer<VariableContext> ctx = var_db->getContext("context");

        SAMRAIPointer<SideVariableNd<double> > u_side_var = make_samrai_shared<SideVariableNd<double> >("u_side");
        SAMRAIPointer<SideVariableNd<double> > f_side_var = make_samrai_shared<SideVariableNd<double> >("f_side");
        SAMRAIPointer<SideVariableNd<double> > e_side_var = make_samrai_shared<SideVariableNd<double> >("e_side");

        const int u_side_idx = var_db->registerVariableAndContext(u_side_var, ctx, IntVectorNd(1));
        const int f_side_idx = var_db->registerVariableAndContext(f_side_var, ctx, IntVectorNd(1));
        const int e_side_idx = var_db->registerVariableAndContext(e_side_var, ctx, IntVectorNd(1));

        SAMRAIPointer<CellVariableNd<double> > u_cell_var = make_samrai_shared<CellVariableNd<double> >("u_cell", NDIM);
        SAMRAIPointer<CellVariableNd<double> > f_cell_var = make_samrai_shared<CellVariableNd<double> >("f_cell", NDIM);
        SAMRAIPointer<CellVariableNd<double> > e_cell_var = make_samrai_shared<CellVariableNd<double> >("e_cell", NDIM);

        const int u_cell_idx = var_db->registerVariableAndContext(u_cell_var, ctx, IntVectorNd(0));
        const int f_cell_idx = var_db->registerVariableAndContext(f_cell_var, ctx, IntVectorNd(0));
        const int e_cell_idx = var_db->registerVariableAndContext(e_cell_var, ctx, IntVectorNd(0));

#if (NDIM == 2)
        SAMRAIPointer<NodeVariableNd<double> > mu_node_var = make_samrai_shared<NodeVariableNd<double> >("mu_node");
        const int mu_node_idx = var_db->registerVariableAndContext(mu_node_var, ctx, IntVectorNd(1));
#endif
#if (NDIM == 3)
        SAMRAIPointer<EdgeVariableNd<double> > mu_edge_var = make_samrai_shared<EdgeVariableNd<double> >("mu_edge");
        const int mu_edge_idx = var_db->registerVariableAndContext(mu_edge_var, ctx, IntVectorNd(1));
        SAMRAIPointer<CellVariableNd<double> > mu_cell_var = make_samrai_shared<CellVariableNd<double> >("mu_cell");
        const int mu_cell_idx = var_db->registerVariableAndContext(mu_cell_var, ctx, IntVectorNd(0));
#endif
        // Register variables for plotting.
        SAMRAIPointer<VisItDataWriterNd> visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);

        visit_data_writer->registerPlotQuantity(u_cell_var->getName(), "VECTOR", u_cell_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(u_cell_var->getName() + std::to_string(d), "SCALAR", u_cell_idx, d);
        }

        visit_data_writer->registerPlotQuantity(f_cell_var->getName(), "VECTOR", f_cell_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(f_cell_var->getName() + std::to_string(d), "SCALAR", f_cell_idx, d);
        }

        visit_data_writer->registerPlotQuantity(e_cell_var->getName(), "VECTOR", e_cell_idx);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            visit_data_writer->registerPlotQuantity(e_cell_var->getName() + std::to_string(d), "SCALAR", e_cell_idx, d);
        }

#if (NDIM == 2)
        visit_data_writer->registerPlotQuantity(mu_node_var->getName(), "SCALAR", mu_node_idx);
#endif
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

        // Set the simulation time to be zero.
        const double data_time = 0.0;

        // Allocate data on each level of the patch hierarchy.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_side_idx, data_time);
            level->allocatePatchData(f_side_idx, data_time);
            level->allocatePatchData(e_side_idx, data_time);
            level->allocatePatchData(u_cell_idx, data_time);
            level->allocatePatchData(f_cell_idx, data_time);
            level->allocatePatchData(e_cell_idx, data_time);
#if (NDIM == 2)
            level->allocatePatchData(mu_node_idx, data_time);
#endif
#if (NDIM == 3)
            level->allocatePatchData(mu_edge_idx, data_time);
            level->allocatePatchData(mu_cell_idx, data_time);
#endif
        }

        // Setup exact solution data.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        muParserCartGridFunction f_fcn("f", app_initializer->getComponentDatabase("f"), grid_geometry);
        muParserCartGridFunction mu_fcn("mu", app_initializer->getComponentDatabase("mu"), grid_geometry);

        u_fcn.setDataOnPatchHierarchy(u_side_idx, u_side_var, patch_hierarchy, data_time);
        f_fcn.setDataOnPatchHierarchy(e_side_idx, e_side_var, patch_hierarchy, data_time);
#if (NDIM == 2)
        mu_fcn.setDataOnPatchHierarchy(mu_node_idx, mu_node_var, patch_hierarchy, data_time);
#endif
#if (NDIM == 3)
        mu_fcn.setDataOnPatchHierarchy(mu_edge_idx, mu_edge_var, patch_hierarchy, data_time);
#endif
        // Create an object to communicate ghost cell data.
        typedef HierarchyGhostCellInterpolation::InterpolationTransactionComponent InterpolationTransactionComponent;
        InterpolationTransactionComponent u_transaction(
            u_side_idx, "CONSERVATIVE_LINEAR_REFINE", true, "CONSERVATIVE_COARSEN", "LINEAR", false);
#if (NDIM == 2)
        InterpolationTransactionComponent mu_transaction(
            mu_node_idx, "LINEAR_REFINE", false, "CONSTANT_COARSEN", "QUADRATIC", false);
        vector<InterpolationTransactionComponent> transactions(2);
        transactions[0] = u_transaction;
        transactions[1] = mu_transaction;
        auto bdry_fill_op = make_samrai_shared<HierarchyGhostCellInterpolation>();
        bdry_fill_op->initializeOperatorState(transactions, patch_hierarchy);
#endif
#if (NDIM == 3)
        InterpolationTransactionComponent mu_transaction(
            mu_edge_idx, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false);
        vector<InterpolationTransactionComponent> transactions(2);
        transactions[0] = u_transaction;
        transactions[1] = mu_transaction;
        auto bdry_fill_op = make_samrai_shared<HierarchyGhostCellInterpolation>();
        bdry_fill_op->initializeOperatorState(transactions, patch_hierarchy);
#endif
        // Fill ghost cells
        bdry_fill_op->fillData(0.0);

        // Create the math operations object and get the patch data index for
        // the side-centered weighting factor.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int dx_side_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

#if (NDIM == 2)
        // Compute f := div mu (grad(u) + grad(u)^T).
        hier_math_ops.vc_laplace(f_side_idx,
                                 f_side_var,
                                 1.0,
                                 0.0,
                                 mu_node_idx,
                                 mu_node_var,
                                 u_side_idx,
                                 u_side_var,
                                 bdry_fill_op,
                                 data_time,
                                 VC_AVERAGE_INTERP);
#endif

#if (NDIM == 3)
        // Compute f := div mu (grad(u) + grad(u)^T).
        hier_math_ops.vc_laplace(f_side_idx,
                                 f_side_var,
                                 1.0,
                                 0.0,
                                 mu_edge_idx,
                                 mu_edge_var,
                                 u_side_idx,
                                 u_side_var,
                                 bdry_fill_op,
                                 data_time,
                                 VC_AVERAGE_INTERP);
#endif
        // Compute error and print error norms.
        SAMRAIPointer<HierarchyDataOpsRealNd<double> > hier_side_data_ops =
            HierarchyDataOpsManagerNd::getManager()->getOperationsDouble(u_side_var, patch_hierarchy, true);
        hier_side_data_ops->subtract(e_side_idx, e_side_idx, f_side_idx); // computes e := e - f
        pout << "|e|_oo = " << hier_side_data_ops->maxNorm(e_side_idx, dx_side_idx) << "\n";
        pout << "|e|_2  = " << hier_side_data_ops->L2Norm(e_side_idx, dx_side_idx) << "\n";
        pout << "|e|_1  = " << hier_side_data_ops->L1Norm(e_side_idx, dx_side_idx) << "\n";

        // Interpolate the side-centered data to cell centers for output.
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(u_cell_idx, u_cell_var, u_side_idx, u_side_var, nullptr, data_time, synch_cf_interface);
        hier_math_ops.interp(f_cell_idx, f_cell_var, f_side_idx, f_side_var, nullptr, data_time, synch_cf_interface);
        hier_math_ops.interp(e_cell_idx, e_cell_var, e_side_idx, e_side_var, nullptr, data_time, synch_cf_interface);

#if (NDIM == 3)
        // Interpolate the edge-centered data to cell centers for output.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<PatchLevelNd> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevelNd::Iterator p(level); p; p++)
            {
                SAMRAIPointer<PatchNd> patch = level->getPatch(p());
                const BoxNd& box = patch->getBox();
                SAMRAIPointer<EdgeDataNd<double> > mu_ec_data = patch->getPatchData(mu_edge_idx);
                SAMRAIPointer<CellDataNd<double> > mu_cc_data = patch->getPatchData(mu_cell_idx);
                for (BoxNd::Iterator it(CellGeometryNd::toCellBox(box)); it; it++)
                {
                    CellIndexNd ci(it());
                    BoxNd edge_box(ci, ci);
                    double avg_mu = 0.0;
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (EdgeIteratorNd e(edge_box, axis); e; e++)
                        {
                            EdgeIndexNd ei(e());
                            avg_mu += (*mu_ec_data)(ei);
                        }
                    }
                    (*mu_cc_data)(ci) = avg_mu / 12.0;
                }
            }
        }
#endif

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, data_time);

    } // cleanup dynamically allocated objects prior to shutdown
} // run_example`
