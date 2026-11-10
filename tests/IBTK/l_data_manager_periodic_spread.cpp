// ---------------------------------------------------------------------
//
// Copyright (c) 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Test the spreading of node- and edge-centered data at points near a periodic boundary. The program spreads the
// value 1 of each point, half of them within two cells of a boundary, with nine kernels. The integral of the
// node-centered grid function over the nodes, each counted once, must be the number of points. In three dimensions,
// the integral of the edge-centered grid function of each axis over the edges must also be the number of points.

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LInitStrategy.h>
#include <ibtk/LMesh.h>
#include <ibtk/LNode.h>
#include <ibtk/LNodeSetData.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>

#include <tbox/Pointer.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellIterator.h>
#include <EdgeData.h>
#include <EdgeVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <NodeData.h>
#include <NodeVariable.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
std::vector<IBTK::Point>
build_points(const int N, const int num_points)
{
    std::mt19937 generator(20261103);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::vector<IBTK::Point> points;
    for (int k = 0; k < num_points; ++k)
    {
        IBTK::Point X;
        for (int d = 0; d < NDIM; ++d)
        {
            X[d] = uniform(generator);
            if (k % 2 == 0)
            {
                // Within two cells of a boundary, on either side.
                const double offset = 2.0 * uniform(generator) / N;
                X[d] = (generator() % 2 == 0) ? offset : std::nextafter(1.0 - offset, 0.0);
            }
        }
        points.push_back(X);
    }
    return points;
}

// Initialize the points at their positions.
class PointInitializer : public IBTK::LInitStrategy
{
public:
    PointInitializer(const std::vector<IBTK::Point>& points) : d_X(points)
    {
    }

    bool getLevelHasLagrangianData(int level_number, bool /*can_be_refined*/) const override
    {
        return level_number == 0;
    }

    bool getIsAllLagrangianDataInDomain(Pointer<PatchHierarchy<NDIM>> /*hierarchy*/) const override
    {
        return true;
    }

    unsigned int computeGlobalNodeCountOnPatchLevel(Pointer<PatchHierarchy<NDIM>> /*hierarchy*/,
                                                    int level_number,
                                                    double /*init_data_time*/,
                                                    bool /*can_be_refined*/,
                                                    bool /*initial_time*/) override
    {
        return level_number == 0 ? d_X.size() : 0;
    }

    unsigned int computeLocalNodeCountOnPatchLevel(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                   int level_number,
                                                   double /*init_data_time*/,
                                                   bool /*can_be_refined*/,
                                                   bool /*initial_time*/) override
    {
        if (level_number != 0) return 0;

        unsigned int count = 0;
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_number);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Box<NDIM>& patch_box = level->getPatch(p())->getBox();
            for (const IBTK::Point& X : d_X)
            {
                if (patch_box.contains(
                        IndexUtilities::getAssignedCellIndex(X, hierarchy->getGridGeometry(), level->getRatio())))
                {
                    ++count;
                }
            }
        }
        return count;
    }

    void initializeStructureIndexingOnPatchLevel(std::map<int, std::string>& strct_id_to_strct_name_map,
                                                 std::map<int, std::pair<int, int>>& strct_id_to_lag_idx_range_map,
                                                 int level_number,
                                                 double /*init_data_time*/,
                                                 bool /*can_be_refined*/,
                                                 bool /*initial_time*/,
                                                 IBTK::LDataManager* /*l_data_manager*/) override
    {
        if (level_number != 0) return;
        strct_id_to_strct_name_map[0] = "points";
        strct_id_to_lag_idx_range_map[0] = std::make_pair(0, static_cast<int>(d_X.size()));
    }

    unsigned int initializeDataOnPatchLevel(int lag_node_index_idx,
                                            unsigned int global_index_offset,
                                            unsigned int local_index_offset,
                                            Pointer<IBTK::LData> X_data,
                                            Pointer<IBTK::LData> /*U_data*/,
                                            Pointer<PatchHierarchy<NDIM>> hierarchy,
                                            int level_number,
                                            double /*init_data_time*/,
                                            bool /*can_be_refined*/,
                                            bool /*initial_time*/,
                                            IBTK::LDataManager* /*l_data_manager*/) override
    {
        if (level_number != 0) return 0;

        boost::multi_array_ref<double, 2>& X_array = *X_data->getLocalFormVecArray();
        int local_idx = -1;
        unsigned int local_node_count = 0;
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_number);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<IBTK::LNodeSetData> index_data = patch->getPatchData(lag_node_index_idx);
            for (std::size_t k = 0; k < d_X.size(); ++k)
            {
                const CellIndex<NDIM> cell_idx =
                    IndexUtilities::getAssignedCellIndex(d_X[k], hierarchy->getGridGeometry(), level->getRatio());
                if (!patch->getBox().contains(cell_idx)) continue;

                const int local_petsc_idx = ++local_idx + static_cast<int>(local_index_offset);
                const int global_petsc_idx = ++local_node_count - 1 + static_cast<int>(global_index_offset);
                for (int d = 0; d < NDIM; ++d) X_array[local_petsc_idx][d] = d_X[k][d];
                if (!index_data->isElement(cell_idx))
                {
                    index_data->appendItemPointer(cell_idx, new IBTK::LNodeSet());
                }
                index_data->getItem(cell_idx)->push_back(new IBTK::LNode(static_cast<int>(k),
                                                                         global_petsc_idx,
                                                                         local_petsc_idx,
                                                                         IntVector<NDIM>(0),
                                                                         IntVector<NDIM>(0),
                                                                         IBTK::Vector::Zero(),
                                                                         IBTK::Vector::Zero()));
            }
        }
        X_data->restoreArrays();
        return static_cast<unsigned int>(local_idx + 1);
    }

private:
    std::vector<IBTK::Point> d_X;
};

// Run the checks for one kernel with a manager of its own, because the manager fixes the spreading kernel.
void
check_kernel(const std::string& kernel,
             Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
             const int ghost_width,
             const int node_idx,
             const int edge_idx,
             const std::vector<IBTK::Point>& points)
{
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);
    const double num_points = static_cast<double>(points.size());
    const double tol = 1.0e-12 * num_points;
    double cell_volume = 1.0;
    for (int d = 0; d < NDIM; ++d)
        cell_volume *= Pointer<CartesianGridGeometry<NDIM>>(patch_hierarchy->getGridGeometry())->getDx()[d];

    IBTK::LDataManager* const l_data_manager = IBTK::LDataManager::getManager(
        "LDataManagerPeriodicSpread_" + kernel, kernel, kernel, false, IntVector<NDIM>(ghost_width), false);
    l_data_manager->registerLInitStrategy(new PointInitializer(points));
    l_data_manager->setPatchHierarchy(patch_hierarchy);
    l_data_manager->setPatchLevels(0, 0);
    l_data_manager->initializeLevelData(patch_hierarchy, 0, 0.0, false, true);
    l_data_manager->resetHierarchyConfiguration(patch_hierarchy, 0, 0);

    Pointer<IBTK::LData> X_data = l_data_manager->getLData(IBTK::LDataManager::POSN_DATA_NAME, 0);
    RobinPhysBdryPatchStrategy* const no_bdry_op = nullptr;

    // The value 1 for each point, spread to node-centered data.
    {
        Pointer<IBTK::LData> F_data = l_data_manager->createLData("F_node", 0, 1);
        VecSet(F_data->getVec(), 1.0);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<NodeData<NDIM, double>> data = level->getPatch(p())->getPatchData(node_idx);
            data->fillAll(0.0);
        }
        l_data_manager->spread(node_idx, F_data, X_data, no_bdry_op, 0, 0.0);
        double sum = 0.0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<NodeData<NDIM, double>> data = patch->getPatchData(node_idx);
            for (CellIterator<NDIM> ic(patch->getBox()); ic; ic++) sum += data->getArrayData()(ic(), 0);
        }
        sum = IBTK_MPI::sumReduction(sum) * cell_volume;
        if (!(std::abs(sum - num_points) <= tol))
        {
            TBOX_ERROR(kernel << ": the integral of the node-centered grid function is " << sum << " for " << num_points
                              << " points of value 1\n");
        }
    }
#if (NDIM == 3)
    // The value 1 in each component for each point, spread to edge-centered data.
    {
        Pointer<IBTK::LData> F_data = l_data_manager->createLData("F_edge", 0, NDIM);
        VecSet(F_data->getVec(), 1.0);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<EdgeData<NDIM, double>> data = level->getPatch(p())->getPatchData(edge_idx);
            data->fillAll(0.0);
        }
        l_data_manager->spread(edge_idx, F_data, X_data, no_bdry_op, 0, 0.0);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            double sum = 0.0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<EdgeData<NDIM, double>> data = patch->getPatchData(edge_idx);
                for (CellIterator<NDIM> ic(patch->getBox()); ic; ic++) sum += data->getArrayData(axis)(ic(), 0);
            }
            sum = IBTK_MPI::sumReduction(sum) * cell_volume;
            if (!(std::abs(sum - num_points) <= tol))
            {
                TBOX_ERROR(kernel << ": the integral of the edge-centered grid function of axis " << axis << " is "
                                  << sum << " for " << num_points << " points of value 1\n");
            }
        }
    }
#else
    NULL_USE(edge_idx);
#endif
    pout << kernel << ": passed\n";
} // check_kernel
} // namespace

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv);
    Logger::getInstance()->setWarning(false);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    // The kernels that set the stencil side by comparing the position with the cell center, then the controls.
    const std::vector<std::string> kernels = { "BSPLINE_4",
                                               "BSPLINE_6",
                                               "PIECEWISE_CUBIC",
                                               "COMPOSITE_BSPLINE_34",
                                               "COMPOSITE_BSPLINE_43",
                                               "COMPOSITE_BSPLINE_56",
                                               "COMPOSITE_BSPLINE_65",
                                               "IB_4",
                                               "PIECEWISE_LINEAR" };
    int ghost_width = 0;
    for (const std::string& kernel : kernels)
    {
        ghost_width = std::max(ghost_width, LEInteractor::getMinimumGhostWidth(kernel));
    }

    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector =
        new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize", nullptr, Pointer<Database>(nullptr));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

    // The grid functions that store the spread values.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("l_data_manager_periodic_spread");
    Pointer<NodeVariable<NDIM, double>> node_var = new NodeVariable<NDIM, double>("node");
    const int node_idx = var_db->registerVariableAndContext(node_var, ctx, IntVector<NDIM>(ghost_width));
    int edge_idx = IBTK::invalid_index;
#if (NDIM == 3)
    Pointer<EdgeVariable<NDIM, double>> edge_var = new EdgeVariable<NDIM, double>("edge");
    edge_idx = var_db->registerVariableAndContext(edge_var, ctx, IntVector<NDIM>(ghost_width));
#endif

    gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);
    level->allocatePatchData(node_idx, 0.0);
    if (edge_idx != IBTK::invalid_index) level->allocatePatchData(edge_idx, 0.0);

    const std::vector<IBTK::Point> points = build_points(input_db->getInteger("N"), NDIM == 2 ? 400 : 600);
    pout << "patches on the level = " << level->getNumberOfPatches() << '\n' << "points = " << points.size() << '\n';
    for (const std::string& kernel : kernels)
        check_kernel(kernel, patch_hierarchy, ghost_width, node_idx, edge_idx, points);
}
