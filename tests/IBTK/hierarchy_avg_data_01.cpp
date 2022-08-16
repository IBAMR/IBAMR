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
#include <ibtk/HierarchyAveragedDataManager.h>
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
               Pointer<PatchHierarchy<NDIM> > hierarchy,
               double time,
               std::mt19937& gen,
               uniform_double_distribution& dis);

void test(const std::string& test_name,
          const int idx,
          const std::unique_ptr<HierarchyAveragedDataManager>& avg_manager,
          Pointer<PatchHierarchy<NDIM> > hierarchy,
          const std::set<double>& times,
          const int max_periods,
          const std::string& refine_type,
          const double noise_max);

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
        for (unsigned int i = 0; i < num_snaps; ++i) time_pts.insert(t_start + static_cast<double>(i) * dt);
        for (unsigned int i = 0; i < num_snaps - 1; ++i)
            interp_pts.insert(t_start + (static_cast<double>(i) + 0.5) * dt);

        const unsigned int max_periods = input_db->getInteger("max_periods");
        const double noise_max = input_db->getDouble("noise_max");

        // First fill in the data
        Pointer<Database> avg_manager_db = new InputDatabase("AvgManager");
        avg_manager_db->putBool("output_data", false);
        avg_manager_db->putBool("enable_logging", true);
        Pointer<CellVariable<NDIM, double> > c_var = new CellVariable<NDIM, double>("c_var");
        std::unique_ptr<HierarchyAveragedDataManager> c_avg_manager(new HierarchyAveragedDataManager(
            "CellAvgManager", c_var, avg_manager_db, patch_hierarchy, time_pts, t_start, t_end, grid_geometry, false));
        c_avg_manager->setSteadyStateThreshold(0.2);
        Pointer<NodeVariable<NDIM, double> > n_var = new NodeVariable<NDIM, double>("n_var");
        std::unique_ptr<HierarchyAveragedDataManager> n_avg_manager(new HierarchyAveragedDataManager(
            "NodeAvgManager", n_var, avg_manager_db, patch_hierarchy, time_pts, t_start, t_end, grid_geometry, false));
        n_avg_manager->setSteadyStateThreshold(0.2);
        Pointer<SideVariable<NDIM, double> > s_var = new SideVariable<NDIM, double>("s_var");
        std::unique_ptr<HierarchyAveragedDataManager> s_avg_manager(new HierarchyAveragedDataManager(
            "SideAvgManager", s_var, avg_manager_db, patch_hierarchy, time_pts, t_start, t_end, grid_geometry, false));
        s_avg_manager->setSteadyStateThreshold(0.2);
        Pointer<EdgeVariable<NDIM, double> > e_var = new EdgeVariable<NDIM, double>("e_var");
        std::unique_ptr<HierarchyAveragedDataManager> e_avg_manager(new HierarchyAveragedDataManager(
            "EdgeAvgManager", e_var, avg_manager_db, patch_hierarchy, time_pts, t_start, t_end, grid_geometry, false));
        e_avg_manager->setSteadyStateThreshold(0.2);
        Pointer<FaceVariable<NDIM, double> > f_var = new FaceVariable<NDIM, double>("f_var");
        std::unique_ptr<HierarchyAveragedDataManager> f_avg_manager(new HierarchyAveragedDataManager(
            "FaceAvgManager", f_var, avg_manager_db, patch_hierarchy, time_pts, t_start, t_end, grid_geometry, false));
        f_avg_manager->setSteadyStateThreshold(0.2);

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        int c_idx = var_db->registerVariableAndContext(c_var, var_db->getContext("ctx"), 1);
        int n_idx = var_db->registerVariableAndContext(n_var, var_db->getContext("ctx"), 1);
        int s_idx = var_db->registerVariableAndContext(s_var, var_db->getContext("ctx"), 1);
        int e_idx = var_db->registerVariableAndContext(e_var, var_db->getContext("ctx"), 1);
        int f_idx = var_db->registerVariableAndContext(f_var, var_db->getContext("ctx"), 1);

        // Do tests
        test("Cell",
             c_idx,
             c_avg_manager,
             patch_hierarchy,
             time_pts,
             max_periods,
             "CONSERVATIVE_LINEAR_REFINE",
             noise_max);
        test("Node", n_idx, n_avg_manager, patch_hierarchy, time_pts, max_periods, "LINEAR_REFINE", noise_max);
        test("Side",
             s_idx,
             s_avg_manager,
             patch_hierarchy,
             time_pts,
             max_periods,
             "CONSERVATIVE_LINEAR_REFINE",
             noise_max);
        test("Edge",
             e_idx,
             e_avg_manager,
             patch_hierarchy,
             time_pts,
             max_periods,
             "CONSERVATIVE_LINEAR_REFINE",
             noise_max);
        test("Face",
             f_idx,
             f_avg_manager,
             patch_hierarchy,
             time_pts,
             max_periods,
             "CONSERVATIVE_LINEAR_REFINE",
             noise_max);
    }
} // main()

void
test(const std::string& test_name,
     const int idx,
     const std::unique_ptr<HierarchyAveragedDataManager>& avg_manager,
     Pointer<PatchHierarchy<NDIM> > hierarchy,
     const std::set<double>& times,
     const int max_periods,
     const std::string& refine_type,
     const double noise_max)
{
    // Fill in data.
    std::mt19937 gen(1);
    uniform_double_distribution dis(-noise_max, noise_max);
    std::vector<std::pair<bool, int> > steady_state_vec(times.size(), std::make_pair(false, 0));
    pout << "Testing " << test_name << "\n";
    for (int period = 0; period < max_periods; ++period)
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
          Pointer<PatchHierarchy<NDIM> > hierarchy,
          const double time,
          std::mt19937& gen,
          uniform_double_distribution& dis)
{
    // Fill idx with the function value
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(idx)) level->allocatePatchData(idx, time);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
            const double* const dx = pgeom->getDx();
            const double* const xlow = pgeom->getXLower();
            const hier::Index<NDIM>& idx_low = patch->getBox().lower();
            Pointer<CellData<NDIM, double> > c_data = patch->getPatchData(idx);
            Pointer<NodeData<NDIM, double> > n_data = patch->getPatchData(idx);
            Pointer<EdgeData<NDIM, double> > e_data = patch->getPatchData(idx);
            Pointer<SideData<NDIM, double> > s_data = patch->getPatchData(idx);
            Pointer<FaceData<NDIM, double> > f_data = patch->getPatchData(idx);

            if (c_data)
            {
                for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
                {
                    const CellIndex<NDIM>& idx = ci();
                    VectorNd x;
                    for (unsigned int d = 0; d < NDIM; ++d)
                        x[d] = xlow[d] + dx[d] * (static_cast<double>(idx(d) - idx_low(d)) + 0.5);
                    (*c_data)(idx) = fcn(x, time, dis(gen));
                }
            }
            else if (n_data)
            {
                for (NodeIterator<NDIM> ni(patch->getBox()); ni; ni++)
                {
                    const NodeIndex<NDIM>& idx = ni();
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
                    for (EdgeIterator<NDIM> ei(patch->getBox(), axis); ei; ei++)
                    {
                        const EdgeIndex<NDIM>& idx = ei();
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
                    for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                    {
                        const SideIndex<NDIM>& idx = si();
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
                    for (FaceIterator<NDIM> fi(patch->getBox(), axis); fi; fi++)
                    {
                        const FaceIndex<NDIM>& idx = fi();
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
