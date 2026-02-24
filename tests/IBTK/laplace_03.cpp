// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2024 by the IBAMR developers
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

#include "ibtk/samrai_compatibility_names.h"

#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// SAMRAI INCLUDES
#include "SAMRAIBergerRigoutsos.h"
#include "SAMRAIBox.h"
#include "SAMRAICartesianGridGeometry.h"
#include "SAMRAICellData.h"
#include "SAMRAICellGeometry.h"
#include "SAMRAICellIndex.h"
#include "SAMRAICellVariable.h"
#include "SAMRAIEdgeData.h"
#include "SAMRAIEdgeIndex.h"
#include "SAMRAIEdgeIterator.h"
#include "SAMRAIEdgeVariable.h"
#include "SAMRAIGriddingAlgorithm.h"
#include "SAMRAIHierarchyDataOpsManager.h"
#include "SAMRAIHierarchyDataOpsReal.h"
#include "SAMRAIIntVector.h"
#include "SAMRAILoadBalancer.h"
#include "SAMRAINodeVariable.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAISideVariable.h"
#include "SAMRAIStandardTagAndInitialize.h"
#include "SAMRAIVariableDatabase.h"
#include "SAMRAIVisItDataWriter.h"

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_enums.h>
#include <ibtk/muParserCartGridFunction.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

// A test program to check that the side-centered variable-coefficient
// generalized Laplace operator discretization yields the expected order of
// accuracy.

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
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "vc_laplace.log");
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database.
        SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = new SAMRAIPatchHierarchy("PatchHierarchy", grid_geometry);
        SAMRAIPointer<SAMRAIStandardTagAndInitialize> error_detector = new SAMRAIStandardTagAndInitialize(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<SAMRAIBergerRigoutsos> box_generator = new SAMRAIBergerRigoutsos();
        SAMRAIPointer<SAMRAILoadBalancer> load_balancer =
            new SAMRAILoadBalancer("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<SAMRAIGriddingAlgorithm> gridding_algorithm =
            new SAMRAIGriddingAlgorithm("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Create variables and register them with the variable database.
        SAMRAIVariableDatabase* var_db = SAMRAIVariableDatabase::getDatabase();
        SAMRAIPointer<VariableContext> ctx = var_db->getContext("context");

        SAMRAIPointer<SAMRAISideVariable<double>> u_side_var = new SAMRAISideVariable<double>("u_side");
        SAMRAIPointer<SAMRAISideVariable<double>> f_side_var = new SAMRAISideVariable<double>("f_side");
        SAMRAIPointer<SAMRAISideVariable<double>> e_side_var = new SAMRAISideVariable<double>("e_side");

        const int u_side_idx = var_db->registerVariableAndContext(u_side_var, ctx, SAMRAIIntVector(1));
        const int f_side_idx = var_db->registerVariableAndContext(f_side_var, ctx, SAMRAIIntVector(1));
        const int e_side_idx = var_db->registerVariableAndContext(e_side_var, ctx, SAMRAIIntVector(1));

        SAMRAIPointer<SAMRAICellVariable<double>> u_cell_var = new SAMRAICellVariable<double>("u_cell", NDIM);
        SAMRAIPointer<SAMRAICellVariable<double>> f_cell_var = new SAMRAICellVariable<double>("f_cell", NDIM);
        SAMRAIPointer<SAMRAICellVariable<double>> e_cell_var = new SAMRAICellVariable<double>("e_cell", NDIM);

        const int u_cell_idx = var_db->registerVariableAndContext(u_cell_var, ctx, SAMRAIIntVector(0));
        const int f_cell_idx = var_db->registerVariableAndContext(f_cell_var, ctx, SAMRAIIntVector(0));
        const int e_cell_idx = var_db->registerVariableAndContext(e_cell_var, ctx, SAMRAIIntVector(0));

#if (NDIM == 2)
        SAMRAIPointer<SAMRAINodeVariable<double>> mu_node_var = new SAMRAINodeVariable<double>("mu_node");
        const int mu_node_idx = var_db->registerVariableAndContext(mu_node_var, ctx, SAMRAIIntVector(1));
#endif
#if (NDIM == 3)
        SAMRAIPointer<SAMRAIEdgeVariable<double>> mu_edge_var = new SAMRAIEdgeVariable<double>("mu_edge");
        const int mu_edge_idx = var_db->registerVariableAndContext(mu_edge_var, ctx, SAMRAIIntVector(1));
        SAMRAIPointer<SAMRAICellVariable<double>> mu_cell_var = new SAMRAICellVariable<double>("mu_cell");
        const int mu_cell_idx = var_db->registerVariableAndContext(mu_cell_var, ctx, SAMRAIIntVector(0));
#endif
        // Register variables for plotting.
        SAMRAIPointer<SAMRAIVisItDataWriter> visit_data_writer = app_initializer->getVisItDataWriter();
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
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
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
        SAMRAIPointer<HierarchyGhostCellInterpolation> bdry_fill_op = new HierarchyGhostCellInterpolation();
        bdry_fill_op->initializeOperatorState(transactions, patch_hierarchy);
#endif
#if (NDIM == 3)
        InterpolationTransactionComponent mu_transaction(
            mu_edge_idx, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false);
        vector<InterpolationTransactionComponent> transactions(2);
        transactions[0] = u_transaction;
        transactions[1] = mu_transaction;
        SAMRAIPointer<HierarchyGhostCellInterpolation> bdry_fill_op = new HierarchyGhostCellInterpolation();
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
        SAMRAIPointer<SAMRAIHierarchyDataOpsReal<double>> hier_side_data_ops =
            SAMRAIHierarchyDataOpsManager::getManager()->getOperationsDouble(u_side_var, patch_hierarchy, true);
        hier_side_data_ops->subtract(e_side_idx, e_side_idx, f_side_idx); // computes e := e - f
        const double max_norm = hier_side_data_ops->maxNorm(e_side_idx, dx_side_idx);
        const double l2_norm = hier_side_data_ops->L2Norm(e_side_idx, dx_side_idx);
        const double l1_norm = hier_side_data_ops->L1Norm(e_side_idx, dx_side_idx);

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "|e|_oo = " << max_norm << "\n";
            out << "|e|_2  = " << l2_norm << "\n";
            out << "|e|_1  = " << l1_norm << "\n";
        }

        // Interpolate the side-centered data to cell centers for output.
        static const bool synch_cf_interface = true;
        hier_math_ops.interp(u_cell_idx, u_cell_var, u_side_idx, u_side_var, nullptr, data_time, synch_cf_interface);
        hier_math_ops.interp(f_cell_idx, f_cell_var, f_side_idx, f_side_var, nullptr, data_time, synch_cf_interface);
        hier_math_ops.interp(e_cell_idx, e_cell_var, e_side_idx, e_side_var, nullptr, data_time, synch_cf_interface);

#if (NDIM == 3)
        // Interpolate the edge-centered data to cell centers for output.
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            SAMRAIPointer<SAMRAIPatchLevel> level = patch_hierarchy->getPatchLevel(ln);
            for (SAMRAIPatchLevel::Iterator p(level); p; p++)
            {
                SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
                const SAMRAIBox& box = patch->getBox();
                SAMRAIPointer<SAMRAIEdgeData<double>> mu_ec_data = patch->getPatchData(mu_edge_idx);
                SAMRAIPointer<SAMRAICellData<double>> mu_cc_data = patch->getPatchData(mu_cell_idx);
                for (SAMRAIBox::Iterator it(SAMRAICellGeometry::toCellBox(box)); it; it++)
                {
                    SAMRAICellIndex ci(it());
                    SAMRAIBox edge_box(ci, ci);
                    double avg_mu = 0.0;
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (SAMRAIEdgeIterator e(edge_box, axis); e; e++)
                        {
                            SAMRAIEdgeIndex ei(e());
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
