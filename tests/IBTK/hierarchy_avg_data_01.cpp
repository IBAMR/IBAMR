// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------
#include <ibtk/samrai_compatibility_names.h>
// SAMRAI INCLUDES
#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyAveragedDataManager.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/SecondaryHierarchy.h>
#include <ibtk/muParserCartGridFunction.h>

#include <SAMRAIBergerRigoutsos.h>
#include <SAMRAICartesianGridGeometry.h>
#include <SAMRAICartesianPatchGeometry.h>
#include <SAMRAICellData.h>
#include <SAMRAICellIndex.h>
#include <SAMRAICellIterator.h>
#include <SAMRAICellVariable.h>
#include <SAMRAIEdgeData.h>
#include <SAMRAIEdgeIndex.h>
#include <SAMRAIEdgeIterator.h>
#include <SAMRAIEdgeVariable.h>
#include <SAMRAIFaceData.h>
#include <SAMRAIFaceIndex.h>
#include <SAMRAIFaceIterator.h>
#include <SAMRAIFaceVariable.h>
#include <SAMRAIGriddingAlgorithm.h>
#include <SAMRAIHierarchyCellDataOpsReal.h>
#include <SAMRAIHierarchyDataOpsManager.h>
#include <SAMRAIHierarchyDataOpsReal.h>
#include <SAMRAIIndex.h>
#include <SAMRAILoadBalancer.h>
#include <SAMRAINodeData.h>
#include <SAMRAINodeIndex.h>
#include <SAMRAINodeIterator.h>
#include <SAMRAINodeVariable.h>
#include <SAMRAIPatch.h>
#include <SAMRAIPatchHierarchy.h>
#include <SAMRAIPatchLevel.h>
#include <SAMRAIPointer.h>
#include <SAMRAISideData.h>
#include <SAMRAISideIndex.h>
#include <SAMRAISideIterator.h>
#include <SAMRAISideVariable.h>
#include <SAMRAIStandardTagAndInitialize.h>
#include <SAMRAIVariableDatabase.h>

#include <array>
#include <random>
#include <utility>

#include <ibtk/app_namespaces.h>

// Need to roll our own distribution to make it machine independent
class uniform_double_distribution
{
public:
    uniform_double_distribution(const double a, const double b) : d_a(a), d_b(b)
    {
        // intentionally blank
    }

    template <class Generator>
    double operator()(Generator& g) const
    {
        return rnd(g, d_a, d_b);
    }

private:
    template <class Generator>
    double rnd(Generator& g, const double a, const double b) const
    {
        double val = static_cast<double>(g() - g.min()) / static_cast<double>(g.max() - g.min()) * (b - a) + a;
        return val;
    }

    double d_a, d_b;
};

void fill_data(int idx,
               SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
               double time,
               std::mt19937& gen,
               uniform_double_distribution& dis);

void test(const std::string& test_name,
          const int idx,
          const std::unique_ptr<HierarchyAveragedDataManager>& avg_manager,
          SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
          const std::set<double>& times,
          const int max_periods,
          const std::string& refine_type,
          const double noise_max,
          const bool do_restart,
          const std::string restart_dir_name);

double
fcn(const VectorNd& x, const double t, const double noise)
{
    double val = 0.0;
    for (int d = 0; d < NDIM; ++d) val += std::sin(2.0 * M_PI * (x[d] - t)) + noise;
    return val;
}

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
        SAMRAIPointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
        SAMRAIPointer<Database> input_db = app_initializer->getInputDatabase();
        std::string restart_file_name = app_initializer->getRestartDumpDirectory();

        // Create major algorithm and data objects that comprise the
        // application. These objects are configured from the input
        // database. Nearly all SAMRAI applications (at least those in IBAMR)
        // start by setting up the same half-dozen objects.
        // Note we generate two patch hierarchies. We want to test our code that interpolates snapshots from one
        // hierarchy onto a different hierarchy.
        SAMRAIPointer<SAMRAICartesianGridGeometry> grid_geometry = new SAMRAICartesianGridGeometry(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        SAMRAIPointer<SAMRAIPatchHierarchy> patch_hierarchy = new SAMRAIPatchHierarchy("PatchHierarchy", grid_geometry);
        SAMRAIPointer<SAMRAIStandardTagAndInitialize> error_detector = new SAMRAIStandardTagAndInitialize(
            "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        SAMRAIPointer<SAMRAIBergerRigoutsos> box_generator = new SAMRAIBergerRigoutsos();
        SAMRAIPointer<SAMRAILoadBalancer> load_balancer =
            new SAMRAILoadBalancer("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        SAMRAIPointer<SAMRAIGriddingAlgorithm> gridding_algorithm =
            new SAMRAIGriddingAlgorithm("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Generate the grid.
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = std::numeric_limits<int>::max();
        int level_number = 0;
        while (gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }

        const double t_start = input_db->getDouble("t_start");
        const double t_end = input_db->getDouble("t_end");
        const unsigned int num_snaps = input_db->getInteger("num_snaps");
        const double dt = (t_end - t_start) / static_cast<double>(num_snaps);
        std::set<double> time_pts, interp_pts;
        for (unsigned int i = 0; i < num_snaps; ++i) time_pts.insert(t_start + (static_cast<double>(i)) * dt);
        for (unsigned int i = 0; i < num_snaps - 1; ++i)
            interp_pts.insert(t_start + (static_cast<double>(i) + 0.5) * dt);

        const unsigned int max_periods = input_db->getInteger("max_periods");
        const double noise_max = input_db->getDouble("noise_max");

        // First fill in the data
        SAMRAIPointer<Database> avg_manager_db = new InputDatabase("AvgManager");
        avg_manager_db->putBool("output_data", false);
        avg_manager_db->putBool("enable_logging", false);
        SAMRAIPointer<SAMRAICellVariable<double>> c_var = new SAMRAICellVariable<double>("c_var");
        std::unique_ptr<HierarchyAveragedDataManager> c_avg_manager(new HierarchyAveragedDataManager(
            "CellAvgManager", c_var, avg_manager_db, time_pts, t_start, t_end, 0.2, grid_geometry, true));
        SAMRAIPointer<SAMRAINodeVariable<double>> n_var = new SAMRAINodeVariable<double>("n_var");
        std::unique_ptr<HierarchyAveragedDataManager> n_avg_manager(new HierarchyAveragedDataManager(
            "NodeAvgManager", n_var, avg_manager_db, time_pts, t_start, t_end, 0.2, grid_geometry, true));
        SAMRAIPointer<SAMRAISideVariable<double>> s_var = new SAMRAISideVariable<double>("s_var");
        std::unique_ptr<HierarchyAveragedDataManager> s_avg_manager(new HierarchyAveragedDataManager(
            "SideAvgManager", s_var, avg_manager_db, time_pts, t_start, t_end, 0.2, grid_geometry, true));
        SAMRAIPointer<SAMRAIEdgeVariable<double>> e_var = new SAMRAIEdgeVariable<double>("e_var");
        std::unique_ptr<HierarchyAveragedDataManager> e_avg_manager(new HierarchyAveragedDataManager(
            "EdgeAvgManager", e_var, avg_manager_db, time_pts, t_start, t_end, 0.2, grid_geometry, true));
        SAMRAIPointer<SAMRAIFaceVariable<double>> f_var = new SAMRAIFaceVariable<double>("f_var");
        std::unique_ptr<HierarchyAveragedDataManager> f_avg_manager(new HierarchyAveragedDataManager(
            "FaceAvgManager", f_var, avg_manager_db, time_pts, t_start, t_end, 0.2, grid_geometry, true));

        auto var_db = SAMRAIVariableDatabase::getDatabase();
        int c_idx = var_db->registerVariableAndContext(c_var, var_db->getContext("ctx"), 1);
        int n_idx = var_db->registerVariableAndContext(n_var, var_db->getContext("ctx"), 1);
        int s_idx = var_db->registerVariableAndContext(s_var, var_db->getContext("ctx"), 1);
        int e_idx = var_db->registerVariableAndContext(e_var, var_db->getContext("ctx"), 1);
        int f_idx = var_db->registerVariableAndContext(f_var, var_db->getContext("ctx"), 1);

        const bool restart = input_db->getBool("RESTART");
        const std::string restart_dir_name = app_initializer->getRestartDumpDirectory();
        if (restart) Utilities::recursiveMkdir(restart_dir_name);

        // Do tests
        test("Cell",
             c_idx,
             c_avg_manager,
             patch_hierarchy,
             time_pts,
             max_periods,
             "CONSERVATIVE_LINEAR_REFINE",
             noise_max,
             restart,
             restart_dir_name);
        test("Node",
             n_idx,
             n_avg_manager,
             patch_hierarchy,
             time_pts,
             max_periods,
             "LINEAR_REFINE",
             noise_max,
             restart,
             restart_dir_name);
        test("Side",
             s_idx,
             s_avg_manager,
             patch_hierarchy,
             time_pts,
             max_periods,
             "CONSERVATIVE_LINEAR_REFINE",
             noise_max,
             restart,
             restart_dir_name);
        test("Edge",
             e_idx,
             e_avg_manager,
             patch_hierarchy,
             time_pts,
             max_periods,
             "CONSERVATIVE_LINEAR_REFINE",
             noise_max,
             restart,
             restart_dir_name);
        test("Face",
             f_idx,
             f_avg_manager,
             patch_hierarchy,
             time_pts,
             max_periods,
             "CONSERVATIVE_LINEAR_REFINE",
             noise_max,
             restart,
             restart_dir_name);
    }
} // main()

void
test(const std::string& test_name,
     const int idx,
     const std::unique_ptr<HierarchyAveragedDataManager>& avg_manager,
     SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
     const std::set<double>& times,
     const int max_periods,
     const std::string& refine_type,
     const double noise_max,
     const bool do_restart,
     const std::string restart_dir_name)
{
    // Fill in data.
    std::mt19937 gen(1);
    uniform_double_distribution dis(-noise_max, noise_max);
    std::vector<std::pair<bool, int>> steady_state_vec(times.size(), std::make_pair(false, 0));
    pout << "Testing " << test_name << "\n";
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    int period = from_restart ? 15 : 0;
    if (from_restart) plog << "Starting after " << period << " periods.\n";
    for (; period < max_periods; ++period)
    {
        int snap_num = 0;
        for (const auto& t : times)
        {
            if (!steady_state_vec[snap_num].first)
            {
                fill_data(idx, hierarchy, t, gen, dis);
                bool at_steady_state =
                    avg_manager->updateTimeAveragedSnapshot(idx, t, hierarchy, refine_type, -1, 1.0e-8);
                if (at_steady_state) steady_state_vec[snap_num].first = true;
                steady_state_vec[snap_num].second++;
            }
            ++snap_num;
        }

        if (do_restart && !from_restart && period == 15)
            RestartManager::getManager()->writeRestartFile(restart_dir_name, period);
        // Determine if we've hit all the steady states
        bool found_all_steady_states = true;
        for (const auto& at_steady_state : steady_state_vec)
            found_all_steady_states = found_all_steady_states && at_steady_state.first;
        if (found_all_steady_states) break;
    }
    for (const auto& steady_state : steady_state_vec)
    {
        if (steady_state.first)
            pout << "Found steady state with " << steady_state.second << " periods\n";
        else
            pout << "Did not find steady state after " << steady_state.second << " periods\n";
    }
}

void
fill_data(const int idx,
          SAMRAIPointer<SAMRAIPatchHierarchy> hierarchy,
          const double time,
          std::mt19937& gen,
          uniform_double_distribution& dis)
{
    // Fill idx with the function value
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        SAMRAIPointer<SAMRAIPatchLevel> level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(idx)) level->allocatePatchData(idx, time);
        for (SAMRAIPatchLevel::Iterator p(level); p; p++)
        {
            SAMRAIPointer<SAMRAIPatch> patch = level->getPatch(p());
            SAMRAIPointer<SAMRAICartesianPatchGeometry> pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double* const xlow = pgeom->getXLower();
            const SAMRAIIndex& idx_low = patch->getBox().lower();
            SAMRAIPointer<SAMRAICellData<double>> c_data = patch->getPatchData(idx);
            SAMRAIPointer<SAMRAINodeData<double>> n_data = patch->getPatchData(idx);
            SAMRAIPointer<SAMRAIEdgeData<double>> e_data = patch->getPatchData(idx);
            SAMRAIPointer<SAMRAISideData<double>> s_data = patch->getPatchData(idx);
            SAMRAIPointer<SAMRAIFaceData<double>> f_data = patch->getPatchData(idx);

            if (c_data)
            {
                for (SAMRAICellIterator ci(patch->getBox()); ci; ci++)
                {
                    const SAMRAICellIndex& idx = ci();
                    VectorNd x;
                    for (unsigned int d = 0; d < NDIM; ++d)
                        x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                    (*c_data)(idx) = fcn(x, time, dis(gen));
                }
            }
            else if (n_data)
            {
                for (SAMRAINodeIterator ni(patch->getBox()); ni; ni++)
                {
                    const SAMRAINodeIndex& idx = ni();
                    VectorNd x;
                    for (unsigned int d = 0; d < NDIM; ++d)
                        x[d] = xlow[d] + dx[d] * static_cast<double>(idx(d) - idx_low(d));
                    (*n_data)(idx) = fcn(x, time, dis(gen));
                }
            }
            else if (e_data)
            {
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SAMRAIEdgeIterator ei(patch->getBox(), axis); ei; ei++)
                    {
                        const SAMRAIEdgeIndex& idx = ei();
                        VectorNd x;
                        for (unsigned int d = 0; d < NDIM; ++d)
                            x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + d == axis ? 0.5 : 0.0);
                        (*e_data)(idx) = fcn(x, time, dis(gen));
                    }
                }
            }
            else if (s_data)
            {
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SAMRAISideIterator si(patch->getBox(), axis); si; si++)
                    {
                        const SAMRAISideIndex& idx = si();
                        VectorNd x;
                        for (unsigned int d = 0; d < NDIM; ++d)
                            x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + d == axis ? 0.0 : 0.5);
                        (*s_data)(idx) = fcn(x, time, dis(gen));
                    }
                }
            }
            else if (f_data)
            {
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (SAMRAIFaceIterator fi(patch->getBox(), axis); fi; fi++)
                    {
                        const SAMRAIFaceIndex& idx = fi();
                        VectorNd x;
                        for (unsigned int d = 0; d < NDIM; ++d)
                            x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + d == axis ? 0.0 : 0.5);
                        (*f_data)(idx) = fcn(x, time, dis(gen));
                    }
                }
            }
        }
    }
}
