// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/CartGridPointwiseFunction.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/muParserCartGridFunction.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SAMRAI_config.h>
#include <StandardTagAndInitialize.h>

#include <ibtk/app_namespaces.h>

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
    // Initialize PETSc, MPI, and SAMRAI.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Suppress a warning
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "adv_diff.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        Pointer<Database> main_db = app_initializer->getComponentDatabase("Main");
        const string reaction_exodus_filename = app_initializer->getExodusIIFilename("reaction");

        // Get various standard options set in the input file.
        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const std::string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        Pointer<CellVariable<NDIM, double> > cc_1_var = new CellVariable<NDIM, double>("cc_1");
        Pointer<CellVariable<NDIM, double> > cc_d_var = new CellVariable<NDIM, double>("cc_d", NDIM);
        Pointer<CellVariable<NDIM, double> > cc_d2_var = new CellVariable<NDIM, double>("cc_d2", NDIM * 2);
        Pointer<NodeVariable<NDIM, double> > nc_var = new NodeVariable<NDIM, double>("nc");
        Pointer<SideVariable<NDIM, double> > sc_var = new SideVariable<NDIM, double>("sc");
        Pointer<FaceVariable<NDIM, double> > fc_var = new FaceVariable<NDIM, double>("fc");
        Pointer<EdgeVariable<NDIM, double> > ec_var = new EdgeVariable<NDIM, double>("ec");

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("CTX");
        IntVector<NDIM> no_ghosts(0);
        const int cc_1_idx = var_db->registerVariableAndContext(cc_1_var, ctx, no_ghosts);
        const int cc_1_cloned_idx = var_db->registerClonedPatchDataIndex(cc_1_var, cc_1_idx);
        const int cc_d_idx = var_db->registerVariableAndContext(cc_d_var, ctx, no_ghosts);
        const int cc_d_cloned_idx = var_db->registerClonedPatchDataIndex(cc_d_var, cc_d_idx);
        const int cc_d2_idx = var_db->registerVariableAndContext(cc_d2_var, ctx, no_ghosts);
        const int cc_d2_cloned_idx = var_db->registerClonedPatchDataIndex(cc_d2_var, cc_d2_idx);
        const int nc_idx = var_db->registerVariableAndContext(nc_var, ctx, no_ghosts);
        const int nc_cloned_idx = var_db->registerClonedPatchDataIndex(nc_var, nc_idx);
        const int sc_idx = var_db->registerVariableAndContext(sc_var, ctx, no_ghosts);
        const int sc_cloned_idx = var_db->registerClonedPatchDataIndex(sc_var, sc_idx);
        const int fc_idx = var_db->registerVariableAndContext(fc_var, ctx, no_ghosts);
        const int fc_cloned_idx = var_db->registerClonedPatchDataIndex(fc_var, fc_idx);
        const int ec_idx = var_db->registerVariableAndContext(ec_var, ctx, no_ghosts);
        const int ec_cloned_idx = var_db->registerClonedPatchDataIndex(ec_var, ec_idx);

        ComponentSelector comps;
        comps.setFlag(cc_1_idx);
        comps.setFlag(cc_1_cloned_idx);
        comps.setFlag(cc_d_idx);
        comps.setFlag(cc_d_cloned_idx);
        comps.setFlag(cc_d2_idx);
        comps.setFlag(cc_d2_cloned_idx);
        comps.setFlag(nc_idx);
        comps.setFlag(nc_cloned_idx);
        comps.setFlag(sc_idx);
        comps.setFlag(sc_cloned_idx);
        comps.setFlag(fc_idx);
        comps.setFlag(fc_cloned_idx);
        comps.setFlag(ec_idx);
        comps.setFlag(ec_cloned_idx);

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        int tag_buffer = 1;
        int ln = 0;
        bool done = false;
        while (!done && (gridding_algorithm->levelCanBeRefined(ln)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(ln);
            ++ln;
        }

        // Read functions from input file
        muParserCartGridFunction mu_scalar_fcn(
            "scalar", app_initializer->getComponentDatabase("ScalarFcn"), grid_geometry);
        muParserCartGridFunction mu_vector_fcn(
            "vector", app_initializer->getComponentDatabase("VectorFcn"), grid_geometry);
        muParserCartGridFunction mu_other_fcn(
            "other", app_initializer->getComponentDatabase("OtherFcn"), grid_geometry);

        PointwiseFunctions::ScalarFcn scalar_fcn = [](const double /*Q*/, const VectorNd& x, const double t) -> double
        { return x(0) + x(1) + t; };

        PointwiseFunctions::VectorFcn vector_fcn =
            [](const VectorNd& /*Q*/, const VectorNd& x, const double t) -> VectorNd
        {
            VectorNd ret;
            for (int d = 0; d < NDIM; ++d) ret[d] = x(0) + x(1) + t;
            return ret;
        };

        PointwiseFunctions::OtherFcn other_fcn = [](const VectorXd& Q, const VectorNd& x, const double t) -> VectorXd
        {
            VectorXd ret = Q;
            for (int d = 0; d < ret.rows(); ++d) ret[d] = x(0) + x(1) + t;
            return ret;
        };

        PointwiseFunctions::StaggeredFcn vector_fcn_side =
            [](const double /*Q*/, const VectorNd& x, const double t, const int /*axis*/) -> double
        { return x(0) + x(1) + t; };

        CartGridPointwiseFunction pt_scalar_fcn("scalar", scalar_fcn);
        CartGridPointwiseFunction pt_vector_fcn("vector", vector_fcn);
        CartGridPointwiseFunction pt_other_fcn("other", other_fcn);
        CartGridPointwiseFunction pt_vector_side_fcn("vector_side", vector_fcn_side);

        // Allocate data
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(comps);
        }

        HierarchyCellDataOpsReal<NDIM, double> hier_cc_data_ops(patch_hierarchy);
        HierarchyNodeDataOpsReal<NDIM, double> hier_nc_data_ops(patch_hierarchy);
        HierarchyFaceDataOpsReal<NDIM, double> hier_fc_data_ops(patch_hierarchy);
        HierarchySideDataOpsReal<NDIM, double> hier_sc_data_ops(patch_hierarchy);
        HierarchyEdgeDataOpsReal<NDIM, double> hier_ec_data_ops(patch_hierarchy);

        const double time = 1.0;

        plog << "Cell centered scalar\n";
        mu_scalar_fcn.setDataOnPatchHierarchy(cc_1_idx, cc_1_var, patch_hierarchy, time, false);
        pt_scalar_fcn.setDataOnPatchHierarchy(cc_1_cloned_idx, cc_1_var, patch_hierarchy, time, false);
        hier_cc_data_ops.subtract(cc_1_idx, cc_1_idx, cc_1_cloned_idx);
        plog << "Max difference: " << hier_cc_data_ops.maxNorm(cc_1_idx, IBTK::invalid_index) << "\n";

        plog << "Cell centered vector\n";
        mu_vector_fcn.setDataOnPatchHierarchy(cc_d_idx, cc_d_var, patch_hierarchy, time, false);
        pt_vector_fcn.setDataOnPatchHierarchy(cc_d_cloned_idx, cc_d_var, patch_hierarchy, time, false);
        hier_cc_data_ops.subtract(cc_d_idx, cc_d_idx, cc_d_cloned_idx);
        plog << "Max difference: " << hier_cc_data_ops.maxNorm(cc_d_idx, IBTK::invalid_index) << "\n";

        plog << "Cell centered other\n";
        mu_other_fcn.setDataOnPatchHierarchy(cc_d2_idx, cc_d2_var, patch_hierarchy, time, false);
        pt_other_fcn.setDataOnPatchHierarchy(cc_d2_cloned_idx, cc_d2_var, patch_hierarchy, time, false);
        hier_cc_data_ops.subtract(cc_d2_idx, cc_d2_idx, cc_d2_cloned_idx);
        plog << "Max difference: " << hier_cc_data_ops.maxNorm(cc_d2_idx, IBTK::invalid_index) << "\n";

        plog << "Node centered\n";
        mu_scalar_fcn.setDataOnPatchHierarchy(nc_idx, nc_var, patch_hierarchy, time, false);
        pt_scalar_fcn.setDataOnPatchHierarchy(nc_cloned_idx, nc_var, patch_hierarchy, time, false);
        hier_nc_data_ops.subtract(nc_idx, nc_idx, nc_cloned_idx);
        plog << "Max difference: " << hier_nc_data_ops.maxNorm(nc_idx, IBTK::invalid_index) << "\n";

        plog << "Side centered\n";
        mu_vector_fcn.setDataOnPatchHierarchy(sc_idx, sc_var, patch_hierarchy, time, false);
        pt_vector_side_fcn.setDataOnPatchHierarchy(sc_cloned_idx, sc_var, patch_hierarchy, time, false);
        hier_sc_data_ops.subtract(sc_idx, sc_idx, sc_cloned_idx);
        plog << "Max difference: " << hier_sc_data_ops.maxNorm(sc_idx, IBTK::invalid_index) << "\n";

        plog << "Face centered\n";
        mu_vector_fcn.setDataOnPatchHierarchy(fc_idx, fc_var, patch_hierarchy, time, false);
        pt_vector_side_fcn.setDataOnPatchHierarchy(fc_cloned_idx, fc_var, patch_hierarchy, time, false);
        hier_fc_data_ops.subtract(fc_idx, fc_idx, fc_cloned_idx);
        plog << "Max difference: " << hier_fc_data_ops.maxNorm(fc_idx, IBTK::invalid_index) << "\n";

        plog << "Edge centered\n";
        mu_vector_fcn.setDataOnPatchHierarchy(ec_idx, ec_var, patch_hierarchy, time, false);
        pt_vector_side_fcn.setDataOnPatchHierarchy(ec_cloned_idx, ec_var, patch_hierarchy, time, false);
        hier_ec_data_ops.subtract(ec_idx, ec_idx, ec_cloned_idx);
        plog << "Max difference: " << hier_ec_data_ops.maxNorm(ec_idx, IBTK::invalid_index) << "\n";

        // Deallocate data
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(comps);
        }

    } // cleanup dynamically allocated objects prior to shutdown
    return 0;
} // main
