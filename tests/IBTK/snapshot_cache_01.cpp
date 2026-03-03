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

// Test verifying that we store, retrieve, and interpolate snapshots.

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/SecondaryHierarchy.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/samrai_compatibility_names.h>
#include <ibtk/snapshot_utilities.h>

#include <SAMRAIBergerRigoutsos.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICellVariable.h>
#include <SAMRAIEdgeVariable.h>
#include <SAMRAIFaceVariable.h>
#include <SAMRAIGridGeometry.h>
#include <SAMRAIGriddingAlgorithm.h>
#include <SAMRAIHierarchyCellDataOpsReal.h>
#include <SAMRAIHierarchyDataOpsManager.h>
#include <SAMRAIHierarchyDataOpsReal.h>
#include <SAMRAILoadBalancer.h>
#include <SAMRAINodeVariable.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAISideVariable.h>
#include <SAMRAIStandardTagAndInitialize.h>
#include <SAMRAIVariable.h>
#include <SAMRAIVariableDatabase.h>

#include <array>
#include <utility>

#include <ibtk/app_namespaces.h>

std::unique_ptr<SnapshotCache> fill_data(const std::string& test_name,
                                         Pointer<SAMRAIVariable> var,
                                         Pointer<SAMRAIPatchHierarchy> hierarchy,
                                         Pointer<muParserCartGridFunction> fcn,
                                         const std::set<double>& time_pts,
                                         bool register_for_restart);

void test_data(Pointer<SAMRAIVariable> var,
               const std::unique_ptr<SnapshotCache>& snapshot_cache,
               Pointer<SAMRAIPatchHierarchy> new_hierarchy,
               Pointer<muParserCartGridFunction> fcn,
               const std::set<double>& interp_pts,
               const std::string& refine_type,
               const int wgt_idx);

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        // Parse command line options, set some standard options from the input
        // file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        std::string restart_file_name = app_initializer->getRestartDumpDirectory();

        // Create major algorithm and data objects that comprise the
        // application. These objects are configured from the input
        // database. Nearly all SAMRAI applications (at least those in IBAMR)
        // start by setting up the same half-dozen objects.
        // Note we generate two patch hierarchies. We want to test our code that interpolates snapshots from one
        // hierarchy onto a different hierarchy.
        Pointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<SAMRAIPatchHierarchy> old_patch_hierarchy =
            new SAMRAIPatchHierarchy("OldPatchHierarchy", grid_geometry);
        Pointer<SAMRAIPatchHierarchy> new_patch_hierarchy =
            new SAMRAIPatchHierarchy("NewPatchHierarchy", grid_geometry);
        Pointer<SAMRAIStandardTagAndInitialize> old_error_detector =
            new SAMRAIStandardTagAndInitialize("OldStandardTagAndInitialize",
                                               nullptr,
                                               app_initializer->getComponentDatabase("OldStandardTagAndInitialize"));
        Pointer<SAMRAIStandardTagAndInitialize> new_error_detector =
            new SAMRAIStandardTagAndInitialize("NewStandardTagAndInitialize",
                                               nullptr,
                                               app_initializer->getComponentDatabase("NewStandardTagAndInitialize"));
        Pointer<SAMRAIBergerRigoutsos> box_generator = new SAMRAIBergerRigoutsos();
        Pointer<SAMRAILoadBalancer> load_balancer =
            new SAMRAILoadBalancer("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<SAMRAIGriddingAlgorithm> old_gridding_algorithm =
            new SAMRAIGriddingAlgorithm("OldGriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        old_error_detector,
                                        box_generator,
                                        load_balancer);
        Pointer<SAMRAIGriddingAlgorithm> new_gridding_alg =
            new SAMRAIGriddingAlgorithm("NewGriddingAlg",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        new_error_detector,
                                        box_generator,
                                        load_balancer);

        // Generate the two grids.
        std::array<std::pair<Pointer<SAMRAIGriddingAlgorithm>, Pointer<SAMRAIPatchHierarchy>>, 2> grid_hier_pairs = {
            std::make_pair(old_gridding_algorithm, old_patch_hierarchy),
            std::make_pair(new_gridding_alg, new_patch_hierarchy)
        };
        for (auto& grid_hier_pair : grid_hier_pairs)
        {
            grid_hier_pair.first->makeCoarsestLevel(grid_hier_pair.second, 0.0);
            const int tag_buffer = std::numeric_limits<int>::max();
            int level_number = 0;
            while (grid_hier_pair.first->levelCanBeRefined(level_number))
            {
                grid_hier_pair.first->makeFinerLevel(grid_hier_pair.second, 0.0, 0.0, tag_buffer);
                ++level_number;
            }
        }

        Pointer<muParserCartGridFunction> fcn =
            new muParserCartGridFunction("Fcn", app_initializer->getComponentDatabase("fcn"), grid_geometry);
        const double t_start = input_db->getDouble("t_start");
        const double t_end = input_db->getDouble("t_end");
        const unsigned int num_snaps = input_db->getInteger("num_snaps");
        const double dt = (t_end - t_start) / static_cast<double>(num_snaps);
        std::set<double> time_pts, interp_pts;
        for (unsigned int i = 0; i < num_snaps; ++i) time_pts.insert(t_start + static_cast<double>(i) * dt);
        for (unsigned int i = 0; i < num_snaps - 1; ++i)
            interp_pts.insert(t_start + (static_cast<double>(i) + 0.5) * dt);

        bool register_for_restart = input_db->getBool("register_for_restart");

        // First fill in the data
        Pointer<SAMRAICellVariable<double>> c_var = new SAMRAICellVariable<double>("c_var");
        std::unique_ptr<SnapshotCache> c_cache = nullptr;
        Pointer<SAMRAINodeVariable<double>> n_var = new SAMRAINodeVariable<double>("n_var");
        std::unique_ptr<SnapshotCache> n_cache = nullptr;
        Pointer<SAMRAISideVariable<double>> s_var = new SAMRAISideVariable<double>("s_var");
        std::unique_ptr<SnapshotCache> s_cache = nullptr;
        Pointer<SAMRAIEdgeVariable<double>> e_var = new SAMRAIEdgeVariable<double>("e_var");
        std::unique_ptr<SnapshotCache> e_cache = nullptr;
        Pointer<SAMRAIFaceVariable<double>> f_var = new SAMRAIFaceVariable<double>("face");
        std::unique_ptr<SnapshotCache> f_cache = nullptr;

        auto var_db = SAMRAIVariableDatabase::getDatabase();
        var_db->registerVariableAndContext(c_var, var_db->getContext("ctx"), 1);
        var_db->registerVariableAndContext(c_var, var_db->getContext("err"));
        var_db->registerVariableAndContext(n_var, var_db->getContext("ctx"), 1);
        var_db->registerVariableAndContext(n_var, var_db->getContext("err"));
        var_db->registerVariableAndContext(s_var, var_db->getContext("ctx"), 1);
        var_db->registerVariableAndContext(s_var, var_db->getContext("err"));
        var_db->registerVariableAndContext(e_var, var_db->getContext("ctx"), 1);
        var_db->registerVariableAndContext(e_var, var_db->getContext("err"));
        var_db->registerVariableAndContext(f_var, var_db->getContext("ctx"), 1);
        var_db->registerVariableAndContext(f_var, var_db->getContext("err"));

        bool from_restart = RestartManager::getManager()->isFromRestart();
        if (RestartManager::getManager()->isFromRestart())
        {
            // Read data from restart.
            // Generate the cache directly
            c_cache = std::make_unique<SnapshotCache>(
                "cell::SnapshotCache", c_var, nullptr, grid_geometry, register_for_restart);
            n_cache = std::make_unique<SnapshotCache>(
                "node::SnapshotCache", n_var, nullptr, grid_geometry, register_for_restart);
            s_cache = std::make_unique<SnapshotCache>(
                "side::SnapshotCache", s_var, nullptr, grid_geometry, register_for_restart);
            e_cache = std::make_unique<SnapshotCache>(
                "edge::SnapshotCache", e_var, nullptr, grid_geometry, register_for_restart);
            f_cache = std::make_unique<SnapshotCache>(
                "face::SnapshotCache", f_var, nullptr, grid_geometry, register_for_restart);
        }
        else
        {
            // Allocate data as normal.
            c_cache = fill_data("cell", c_var, old_patch_hierarchy, fcn, time_pts, register_for_restart);
            n_cache = fill_data("node", n_var, old_patch_hierarchy, fcn, time_pts, register_for_restart);
            s_cache = fill_data("side", s_var, old_patch_hierarchy, fcn, time_pts, register_for_restart);
            e_cache = fill_data("edge", e_var, old_patch_hierarchy, fcn, time_pts, register_for_restart);
            f_cache = fill_data("face", f_var, old_patch_hierarchy, fcn, time_pts, register_for_restart);
        }

        // Write restart files if we need to
        if (register_for_restart && !from_restart) RestartManager::getManager()->writeRestartFile(restart_file_name, 1);

        // Now test data
        pout << "Testing with cell variable\n";
        test_data(c_var, c_cache, new_patch_hierarchy, fcn, interp_pts, "CONSERVATIVE_LINEAR_REFINE", -1);

        pout << "Testing with node variable\n";
        test_data(n_var, n_cache, new_patch_hierarchy, fcn, interp_pts, "LINEAR_REFINE", -1);

        pout << "Testing with side variable\n";
        test_data(s_var, s_cache, new_patch_hierarchy, fcn, interp_pts, "CONSERVATIVE_LINEAR_REFINE", -1);

        pout << "Testing with edge variable\n";
        test_data(e_var, e_cache, new_patch_hierarchy, fcn, interp_pts, "CONSERVATIVE_LINEAR_REFINE", -1);

        pout << "Testing with face variable\n";
        test_data(f_var, f_cache, new_patch_hierarchy, fcn, interp_pts, "CONSERVATIVE_LINEAR_REFINE", -1);
    }
} // main()

std::unique_ptr<SnapshotCache>
fill_data(const std::string& test_name,
          Pointer<SAMRAIVariable> var,
          Pointer<SAMRAIPatchHierarchy> hierarchy,
          Pointer<muParserCartGridFunction> fcn,
          const std::set<double>& time_pts,
          bool register_for_restart)
{
    // Actually do the test.
    // Create the index that we fill with data
    auto var_db = SAMRAIVariableDatabase::getDatabase();
    int var_idx = var_db->registerVariableAndContext(var, var_db->getContext("ctx"), 1 /*ghosts*/);

    // Create a SnapshotCache to store snapshots on the "old" hierarchy.
    Pointer<SAMRAIGridGeometry> grid_geom = hierarchy->getGridGeometry();
    std::unique_ptr<SnapshotCache> snapshot_cache =
        std::make_unique<SnapshotCache>(test_name + "::SnapshotCache", var, nullptr, grid_geom, register_for_restart);

    // Fill in snapshot cache with several values.
    for (const auto& t : time_pts)
    {
        // Allocate patch data
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<SAMRAIPatchLevel> level = hierarchy->getPatchLevel(ln);
            level->allocatePatchData(var_idx, t);
        }
        // Fill in patch data
        fcn->setDataOnPatchHierarchy(var_idx, var, hierarchy, t);

        snapshot_cache->storeSnapshot(var_idx, t, hierarchy);
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<SAMRAIPatchLevel> level = hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(var_idx);
        }
    }

    return snapshot_cache;
}

void
test_data(Pointer<SAMRAIVariable> var,
          const std::unique_ptr<SnapshotCache>& snapshot_cache,
          Pointer<SAMRAIPatchHierarchy> new_hierarchy,
          Pointer<muParserCartGridFunction> fcn,
          const std::set<double>& interp_pts,
          const std::string& refine_type,
          const int wgt_idx)
{
    // Actually do the test.
    // Create the index that we fill with data
    auto var_db = SAMRAIVariableDatabase::getDatabase();
    int var_idx = var_db->registerVariableAndContext(var, var_db->getContext("ctx"), 1 /*ghosts*/);
    int err_idx = var_db->registerVariableAndContext(var, var_db->getContext("err"));

    // Create a SnapshotCache to store snapshots on the "old" hierarchy.
    auto hier_math_ops = SAMRAIHierarchyDataOpsManager::getManager();
    Pointer<SAMRAIHierarchyDataOpsReal<double>> hier_data_ops = hier_math_ops->getOperationsDouble(var, new_hierarchy);

    // Interpolate those values to other points on the "new" hierarchy
    for (const auto& t : interp_pts)
    {
        for (int ln = 0; ln <= new_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<SAMRAIPatchLevel> level = new_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(err_idx, t);
            level->allocatePatchData(var_idx, t);
        }
        fill_snapshot_at_time(*snapshot_cache, var_idx, t, err_idx, new_hierarchy, refine_type, hier_data_ops);

        // Compute and output error
        fcn->setDataOnPatchHierarchy(err_idx, var, new_hierarchy, t);
        hier_data_ops->subtract(err_idx, var_idx, err_idx);
        pout << "  L1-norm:  " << hier_data_ops->L1Norm(err_idx, wgt_idx) << "\n";
        pout << "  L2-norm:  " << hier_data_ops->L2Norm(err_idx, wgt_idx) << "\n";
        pout << "  max-norm: " << hier_data_ops->maxNorm(err_idx, wgt_idx) << "\n";

        for (int ln = 0; ln <= new_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<SAMRAIPatchLevel> level = new_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(err_idx);
            level->deallocatePatchData(var_idx);
        }
    }
}
