// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2022 by the IBAMR developers
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
#include <ibtk/snapshot_utilities.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchyDataOpsManager.h>
#include <HierarchyDataOpsReal.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <array>
#include <utility>

#include <ibtk/app_namespaces.h>

std::unique_ptr<SnapshotCache> fill_data(const std::string& test_name,
                                         Pointer<Variable<NDIM> > var,
                                         Pointer<PatchHierarchy<NDIM> > hierarchy,
                                         Pointer<muParserCartGridFunction> fcn,
                                         const std::set<double>& time_pts,
                                         bool register_for_restart);

void test_data(Pointer<Variable<NDIM> > var,
               const std::unique_ptr<SnapshotCache>& snapshot_cache,
               Pointer<PatchHierarchy<NDIM> > new_hierarchy,
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
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > old_patch_hierarchy =
            new PatchHierarchy<NDIM>("OldPatchHierarchy", grid_geometry);
        Pointer<PatchHierarchy<NDIM> > new_patch_hierarchy =
            new PatchHierarchy<NDIM>("NewPatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > old_error_detector = new StandardTagAndInitialize<NDIM>(
            "OldStandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("OldStandardTagAndInitialize"));
        Pointer<StandardTagAndInitialize<NDIM> > new_error_detector = new StandardTagAndInitialize<NDIM>(
            "NewStandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("NewStandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > old_gridding_algorithm =
            new GriddingAlgorithm<NDIM>("OldGriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        old_error_detector,
                                        box_generator,
                                        load_balancer);
        Pointer<GriddingAlgorithm<NDIM> > new_gridding_alg =
            new GriddingAlgorithm<NDIM>("NewGriddingAlg",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        new_error_detector,
                                        box_generator,
                                        load_balancer);

        // Generate the two grids.
        std::array<std::pair<Pointer<GriddingAlgorithm<NDIM> >, Pointer<PatchHierarchy<NDIM> > >, 2> grid_hier_pairs = {
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
        Pointer<CellVariable<NDIM, double> > c_var = new CellVariable<NDIM, double>("c_var");
        std::unique_ptr<SnapshotCache> c_cache = nullptr;
        Pointer<NodeVariable<NDIM, double> > n_var = new NodeVariable<NDIM, double>("n_var");
        std::unique_ptr<SnapshotCache> n_cache = nullptr;
        Pointer<SideVariable<NDIM, double> > s_var = new SideVariable<NDIM, double>("s_var");
        std::unique_ptr<SnapshotCache> s_cache = nullptr;
        Pointer<EdgeVariable<NDIM, double> > e_var = new EdgeVariable<NDIM, double>("e_var");
        std::unique_ptr<SnapshotCache> e_cache = nullptr;
        Pointer<FaceVariable<NDIM, double> > f_var = new FaceVariable<NDIM, double>("face");
        std::unique_ptr<SnapshotCache> f_cache = nullptr;

        auto var_db = VariableDatabase<NDIM>::getDatabase();
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
            c_cache.reset(
                new SnapshotCache("cell::SnapshotCache", c_var, nullptr, grid_geometry, register_for_restart));
            n_cache.reset(
                new SnapshotCache("node::SnapshotCache", n_var, nullptr, grid_geometry, register_for_restart));
            s_cache.reset(
                new SnapshotCache("side::SnapshotCache", s_var, nullptr, grid_geometry, register_for_restart));
            e_cache.reset(
                new SnapshotCache("edge::SnapshotCache", e_var, nullptr, grid_geometry, register_for_restart));
            f_cache.reset(
                new SnapshotCache("face::SnapshotCache", f_var, nullptr, grid_geometry, register_for_restart));
        }
        else
        {
            // Allocate data as normal.
            c_cache = std::move(fill_data("cell", c_var, old_patch_hierarchy, fcn, time_pts, register_for_restart));
            n_cache = std::move(fill_data("node", n_var, old_patch_hierarchy, fcn, time_pts, register_for_restart));
            s_cache = std::move(fill_data("side", s_var, old_patch_hierarchy, fcn, time_pts, register_for_restart));
            e_cache = std::move(fill_data("edge", e_var, old_patch_hierarchy, fcn, time_pts, register_for_restart));
            f_cache = std::move(fill_data("face", f_var, old_patch_hierarchy, fcn, time_pts, register_for_restart));
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
          Pointer<Variable<NDIM> > var,
          Pointer<PatchHierarchy<NDIM> > hierarchy,
          Pointer<muParserCartGridFunction> fcn,
          const std::set<double>& time_pts,
          bool register_for_restart)
{
    // Actually do the test.
    // Create the index that we fill with data
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    int var_idx = var_db->registerVariableAndContext(var, var_db->getContext("ctx"), 1 /*ghosts*/);

    // Create a SnapshotCache to store snapshots on the "old" hierarchy.
    Pointer<GridGeometry<NDIM> > grid_geom = hierarchy->getGridGeometry();
    std::unique_ptr<SnapshotCache> snapshot_cache(
        new SnapshotCache(test_name + "::SnapshotCache", var, nullptr, grid_geom, register_for_restart));

    // Fill in snapshot cache with several values.
    for (const auto& t : time_pts)
    {
        // Allocate patch data
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            level->allocatePatchData(var_idx, t);
        }
        // Fill in patch data
        fcn->setDataOnPatchHierarchy(var_idx, var, hierarchy, t);

        snapshot_cache->storeSnapshot(var_idx, t, hierarchy);
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(var_idx);
        }
    }

    return snapshot_cache;
}

void
test_data(Pointer<Variable<NDIM> > var,
          const std::unique_ptr<SnapshotCache>& snapshot_cache,
          Pointer<PatchHierarchy<NDIM> > new_hierarchy,
          Pointer<muParserCartGridFunction> fcn,
          const std::set<double>& interp_pts,
          const std::string& refine_type,
          const int wgt_idx)
{
    // Actually do the test.
    // Create the index that we fill with data
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    int var_idx = var_db->registerVariableAndContext(var, var_db->getContext("ctx"), 1 /*ghosts*/);
    int err_idx = var_db->registerVariableAndContext(var, var_db->getContext("err"));

    // Create a SnapshotCache to store snapshots on the "old" hierarchy.
    auto hier_math_ops = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double> > hier_data_ops = hier_math_ops->getOperationsDouble(var, new_hierarchy);

    // Interpolate those values to other points on the "new" hierarchy
    for (const auto& t : interp_pts)
    {
        for (int ln = 0; ln <= new_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = new_hierarchy->getPatchLevel(ln);
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
            Pointer<PatchLevel<NDIM> > level = new_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(err_idx);
            level->deallocatePatchData(var_idx);
        }
    }
}
