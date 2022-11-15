// ---------------------------------------------------------------------
//
// Copyright (c) 2022 - 2022 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Check the 'interpolate side/face to nodal' routines.

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
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "curl_01.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        const std::string src_var_type = input_db->getStringWithDefault("src_var_type", "SIDE");
        const std::string dst_var_type = input_db->getStringWithDefault("dst_var_type", "NODE");

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

        Pointer<Variable<NDIM> > u_var;
        Pointer<SideVariable<NDIM, double> > u_sc_var = new SideVariable<NDIM, double>("u_sc");

        if (src_var_type == "SIDE")
        {
            u_var = u_sc_var;
        }
        else
        {
            TBOX_ERROR("not implemented");
        }

        const bool fine_boundary_represents_var = input_db->getBoolWithDefault("fine_boundary_represents_var", false);
        const unsigned int curl_dim = (NDIM == 2 ? 1 : NDIM);
        Pointer<Variable<NDIM> > curl_u_var, e_var;
        Pointer<NodeVariable<NDIM, double> > curl_u_nc_var =
            new NodeVariable<NDIM, double>("curl_u_nc", curl_dim, fine_boundary_represents_var);
        Pointer<NodeVariable<NDIM, double> > e_nc_var =
            new NodeVariable<NDIM, double>("e_nc", curl_dim, fine_boundary_represents_var);

        if (dst_var_type == "NODE")
        {
            curl_u_var = curl_u_nc_var;
            e_var = e_nc_var;
        }
        else
        {
            TBOX_ERROR("not implemented");
        }

        const int u_idx = var_db->registerVariableAndContext(u_var, ctx, IntVector<NDIM>(1));
        const int curl_u_idx = var_db->registerVariableAndContext(curl_u_var, ctx, IntVector<NDIM>(0));
        const int e_idx = var_db->registerVariableAndContext(e_var, ctx, IntVector<NDIM>(0));

        // Register variables for plotting.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        TBOX_ASSERT(visit_data_writer);
        if (dst_var_type == "NODE")
        {
            if (curl_dim > 1)
            {
                visit_data_writer->registerPlotQuantity(curl_u_var->getName(), "VECTOR", curl_u_idx);
                visit_data_writer->registerPlotQuantity(e_var->getName(), "VECTOR", e_idx);
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    visit_data_writer->registerPlotQuantity(
                        curl_u_var->getName() + std::to_string(d), "SCALAR", curl_u_idx, d);
                    visit_data_writer->registerPlotQuantity(e_var->getName() + std::to_string(d), "SCALAR", e_idx, d);
                }
            }
            else
            {
                visit_data_writer->registerPlotQuantity(curl_u_var->getName(), "SCALAR", curl_u_idx);
                visit_data_writer->registerPlotQuantity(e_var->getName(), "SCALAR", e_idx);
            }
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
            level->allocatePatchData(curl_u_idx, 0.0);
            level->allocatePatchData(e_idx, 0.0);
        }

        // Setup vector objects.
        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);

        SAMRAIVectorReal<NDIM, double> curl_u_vec(
            "curl u", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        SAMRAIVectorReal<NDIM, double> e_vec("e", patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        curl_u_vec.addComponent(curl_u_var, curl_u_idx);
        e_vec.addComponent(e_var, e_idx);

        curl_u_vec.setToScalar(0.0);
        e_vec.setToScalar(0.0);

        // Setup exact solutions.
        muParserCartGridFunction u_fcn("u", app_initializer->getComponentDatabase("u"), grid_geometry);
        u_fcn.setDataOnPatchHierarchy(u_idx, u_var, patch_hierarchy, 0.0);
        muParserCartGridFunction curl_u_fcn("curl u", app_initializer->getComponentDatabase("curl_u"), grid_geometry);
        curl_u_fcn.setDataOnPatchHierarchy(e_idx, e_var, patch_hierarchy, 0.0);

        // We need updated ghost values to evaluate the curl:
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> output_transaction_comps;
        output_transaction_comps.emplace_back(
            u_idx, "CONSERVATIVE_LINEAR_REFINE", false, "CONSERVATIVE_COARSEN", "LINEAR", false);
        Pointer<HierarchyGhostCellInterpolation> hier_bdry_fill = new HierarchyGhostCellInterpolation();
        hier_bdry_fill->initializeOperatorState(output_transaction_comps, patch_hierarchy);
        hier_bdry_fill->fillData(0.0);

        // evaluate the curl:
        const bool synch_src_cf_interface = true;
        const bool synch_dst_cf_interface = input_db->getBoolWithDefault("synch_dst_cf_interface", false);
        if (dst_var_type == "NODE" && src_var_type == "SIDE")
        {
            hier_math_ops.curl(
                curl_u_idx, curl_u_nc_var, synch_dst_cf_interface, u_idx, u_sc_var, NULL, synch_src_cf_interface, 0.0);
        }
        else
        {
            TBOX_ERROR("not implemented");
        }

        // Compute error and print error norms.
        e_vec.subtract(Pointer<SAMRAIVectorReal<NDIM, double> >(&e_vec, false),
                       Pointer<SAMRAIVectorReal<NDIM, double> >(&curl_u_vec, false));
        const double max_norm = e_vec.maxNorm();

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            out << "|e|_oo = " << max_norm << "\n";
        }

        // Output data for plotting.
        visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);

    } // cleanup dynamically allocated objects prior to shutdown
} // run_example
