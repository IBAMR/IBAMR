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

// Test verifying that we store, retrieve, and interpolate snapshots.

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyTimeInterpolator.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/SecondaryHierarchy.h>
#include <ibtk/muParserCartGridFunction.h>

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

template <class VariableType>
void test(const std::string& test_name,
          Pointer<PatchHierarchy<NDIM> > cur_hierarchy,
          Pointer<PatchHierarchy<NDIM> > new_hierarchy,
          Pointer<muParserCartGridFunction> fcn,
          const std::set<double>& time_pts,
          const std::set<double>& interp_pts,
          const double t_start,
          const double t_end,
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
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cc_laplace.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Create major algorithm and data objects that comprise the
        // application. These objects are configured from the input
        // database. Nearly all SAMRAI applications (at least those in IBAMR)
        // start by setting up the same half-dozen objects.
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
        for (unsigned int i = 0; i < num_snaps; ++i) interp_pts.insert(t_start + (static_cast<double>(i) + 0.5) * dt);

        pout << "Testing with cell variable\n";
        test<CellVariable<NDIM, double> >("cell",
                                          old_patch_hierarchy,
                                          new_patch_hierarchy,
                                          fcn,
                                          time_pts,
                                          interp_pts,
                                          t_start,
                                          t_end,
                                          "CONSERVATIVE_LINEAR_REFINE",
                                          -1);
        pout << "Testing with node variable\n";
        test<NodeVariable<NDIM, double> >("node",
                                          old_patch_hierarchy,
                                          new_patch_hierarchy,
                                          fcn,
                                          time_pts,
                                          interp_pts,
                                          t_start,
                                          t_end,
                                          "LINEAR_REFINE",
                                          -1);
        pout << "Testing with side variable\n";
        test<SideVariable<NDIM, double> >("side",
                                          old_patch_hierarchy,
                                          new_patch_hierarchy,
                                          fcn,
                                          time_pts,
                                          interp_pts,
                                          t_start,
                                          t_end,
                                          "CONSERVATIVE_LINEAR_REFINE",
                                          -1);
        pout << "Testing with edge variable\n";
        test<EdgeVariable<NDIM, double> >("edge",
                                          old_patch_hierarchy,
                                          new_patch_hierarchy,
                                          fcn,
                                          time_pts,
                                          interp_pts,
                                          t_start,
                                          t_end,
                                          "CONSERVATIVE_LINEAR_REFINE",
                                          -1);
        pout << "Testing with face variable\n";
        test<FaceVariable<NDIM, double> >("face",
                                          old_patch_hierarchy,
                                          new_patch_hierarchy,
                                          fcn,
                                          time_pts,
                                          interp_pts,
                                          t_start,
                                          t_end,
                                          "CONSERVATIVE_LINEAR_REFINE",
                                          -1);
    }
} // main()

template <class VariableType>
void
test(const std::string& test_name,
     Pointer<PatchHierarchy<NDIM> > cur_hierarchy,
     Pointer<PatchHierarchy<NDIM> > new_hierarchy,
     Pointer<muParserCartGridFunction> fcn,
     const std::set<double>& time_pts,
     const std::set<double>& interp_pts,
     const double t_start,
     const double t_end,
     const std::string& refine_type,
     const int wgt_idx)
{
    // Actually do the test.
    // Create the index that we fill with data
    Pointer<VariableType> var = new VariableType(test_name);
    auto var_db = VariableDatabase<NDIM>::getDatabase();
    int var_idx = var_db->registerVariableAndContext(var, var_db->getContext("ctx"), 1 /*ghosts*/);
    int err_idx = var_db->registerVariableAndContext(var, var_db->getContext("err"));

    // Create a hierarchy interpolator.
    std::unique_ptr<HierarchyTimeInterpolator<VariableType> > interpolator(new HierarchyTimeInterpolator<VariableType>(
        test_name + "::HierarchyTimeInterpolator", new_hierarchy, time_pts, t_start, t_end));

    // Fill in snapshot cache with several values.
    for (const auto& t : time_pts)
    {
        // Allocate patch data
        for (int ln = 0; ln <= cur_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = cur_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(var_idx, t);
        }
        // Fill in patch data
        fcn->setDataOnPatchHierarchy(var_idx, var, cur_hierarchy, t);

        interpolator->updateTimeAveragedSnapshot(var_idx, t, cur_hierarchy, refine_type, -1, 1.0e-8);
        for (int ln = 0; ln <= cur_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = cur_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(var_idx);
        }
    }

    // Interpolate those values to other points
    auto hier_math_ops = HierarchyDataOpsManager<NDIM>::getManager();
    Pointer<HierarchyDataOpsReal<NDIM, double> > hier_data_ops = hier_math_ops->getOperationsDouble(var, new_hierarchy);
    // Now interpolate these values to var_idx
    for (const auto& t : interp_pts)
    {
        for (int ln = 0; ln <= new_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM> > level = new_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(err_idx, t);
            level->allocatePatchData(var_idx, t);
        }
        interpolator->fillSnapshotAtTime(var_idx, t, new_hierarchy, refine_type);

        // Compute and output error
        hier_data_ops->setToScalar(err_idx, 0.0);
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
