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

// Test the single-level LDataManager::interp() and spread() on each level of a two-level hierarchy with markers on
// both levels. The second marker of the coarser level belongs to an inactive structure.

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LInitStrategy.h>
#include <ibtk/LMesh.h>
#include <ibtk/LNode.h>
#include <ibtk/LNodeSetData.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>

#include <tbox/Pointer.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <ibtk/app_namespaces.h>

namespace
{
// Two markers on each level, each in a structure of its own. The markers of level 1 lie inside the refined region.
std::vector<std::vector<IBTK::Point>>
get_marker_positions()
{
    return { { IBTK::Point(0.15, 0.85), IBTK::Point(0.85, 0.15) },
             { IBTK::Point(0.45, 0.50), IBTK::Point(0.55, 0.53) } };
}

class MarkerInitializer : public IBTK::LInitStrategy
{
public:
    bool getLevelHasLagrangianData(int level_number, bool /*can_be_refined*/) const override
    {
        return level_number < static_cast<int>(d_X.size());
    }

    bool getIsAllLagrangianDataInDomain(Pointer<PatchHierarchy<NDIM>> /*hierarchy*/) const override
    {
        return true;
    }

    unsigned int computeGlobalNodeCountOnPatchLevel(Pointer<PatchHierarchy<NDIM>> /*hierarchy*/,
                                                    int level_number,
                                                    double /*init_data_time*/,
                                                    bool can_be_refined,
                                                    bool /*initial_time*/) override
    {
        return getLevelHasLagrangianData(level_number, can_be_refined) ? d_X[level_number].size() : 0;
    }

    unsigned int computeLocalNodeCountOnPatchLevel(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                   int level_number,
                                                   double /*init_data_time*/,
                                                   bool can_be_refined,
                                                   bool /*initial_time*/) override
    {
        if (!getLevelHasLagrangianData(level_number, can_be_refined)) return 0;

        unsigned int count = 0;
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_number);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Box<NDIM>& patch_box = level->getPatch(p())->getBox();
            for (const IBTK::Point& X : d_X[level_number])
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
                                                 bool can_be_refined,
                                                 bool /*initial_time*/,
                                                 IBTK::LDataManager* /*l_data_manager*/) override
    {
        if (!getLevelHasLagrangianData(level_number, can_be_refined)) return;
        for (int k = 0; k < static_cast<int>(d_X[level_number].size()); ++k)
        {
            strct_id_to_strct_name_map[k] = "marker_" + std::to_string(k);
            strct_id_to_lag_idx_range_map[k] = std::make_pair(k, k + 1);
        }
    }

    unsigned int initializeDataOnPatchLevel(int lag_node_index_idx,
                                            unsigned int global_index_offset,
                                            unsigned int local_index_offset,
                                            Pointer<IBTK::LData> X_data,
                                            Pointer<IBTK::LData> /*U_data*/,
                                            Pointer<PatchHierarchy<NDIM>> hierarchy,
                                            int level_number,
                                            double /*init_data_time*/,
                                            bool can_be_refined,
                                            bool /*initial_time*/,
                                            IBTK::LDataManager* /*l_data_manager*/) override
    {
        if (!getLevelHasLagrangianData(level_number, can_be_refined)) return 0;

        boost::multi_array_ref<double, 2>& X_array = *X_data->getLocalFormVecArray();
        int local_idx = -1;
        unsigned int local_node_count = 0;
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_number);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<IBTK::LNodeSetData> index_data = patch->getPatchData(lag_node_index_idx);
            for (std::size_t k = 0; k < d_X[level_number].size(); ++k)
            {
                const IBTK::Point& X = d_X[level_number][k];
                const CellIndex<NDIM> cell_idx =
                    IndexUtilities::getAssignedCellIndex(X, hierarchy->getGridGeometry(), level->getRatio());
                if (!patch->getBox().contains(cell_idx)) continue;

                const int local_petsc_idx = ++local_idx + static_cast<int>(local_index_offset);
                const int global_petsc_idx = ++local_node_count - 1 + static_cast<int>(global_index_offset);
                for (int d = 0; d < NDIM; ++d) X_array[local_petsc_idx][d] = X[d];
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
    std::vector<std::vector<IBTK::Point>> d_X = get_marker_positions();
};

// A linear function, which the IB_4 kernel interpolates exactly.
double
linear_field(const double* x)
{
    return 1.0 + x[0] + 2.0 * x[1];
}

// Return the entries of a quantity with one value per marker of the level, indexed by the Lagrangian index.
std::vector<double>
gather_by_lagrangian_index(Pointer<IBTK::LData> Q_data, IBTK::LDataManager* l_data_manager, const int ln)
{
    std::vector<double> Q(l_data_manager->getNumberOfNodes(ln), 0.0);
    boost::multi_array_ref<double, 2>& Q_array = *Q_data->getLocalFormVecArray();
    for (const IBTK::LNode* node : l_data_manager->getLMesh(ln)->getLocalNodes())
    {
        Q[node->getLagrangianIndex()] = Q_array[node->getLocalPETScIndex()][0];
    }
    Q_data->restoreArrays();
    IBTK_MPI::sumReduction(Q.data(), static_cast<int>(Q.size()));
    return Q;
}

// Integrate a cell-centered grid function over the interior of the level.
double
integrate(Pointer<PatchLevel<NDIM>> level, const int idx)
{
    double integral = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, double>> q_data = patch->getPatchData(idx);
        const double* const dx = Pointer<CartesianPatchGeometry<NDIM>>(patch->getPatchGeometry())->getDx();
        double cell_volume = 1.0;
        for (int d = 0; d < NDIM; ++d) cell_volume *= dx[d];
        for (CellIterator<NDIM> ic(patch->getBox()); ic; ic++) integral += (*q_data)(ic()) * cell_volume;
    }
    return IBTK_MPI::sumReduction(integral);
}
} // namespace

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv);
    Logger::getInstance()->setWarning(false);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);

    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

    // The interpolated field and the spread values.
    IBTK::LDataManager* const l_data_manager = IBTK::LDataManager::getManager("LDataManagerTwoLevels", "IB_4", "IB_4");
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("l_data_manager_two_levels");
    Pointer<CellVariable<NDIM, double>> q_var = new CellVariable<NDIM, double>("q");
    Pointer<CellVariable<NDIM, double>> f_var = new CellVariable<NDIM, double>("f");
    const int q_idx = var_db->registerVariableAndContext(q_var, ctx, l_data_manager->getGhostCellWidth());
    const int f_idx = var_db->registerVariableAndContext(f_var, ctx, l_data_manager->getGhostCellWidth());

    gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
    for (int ln = 0; gridding_algorithm->levelCanBeRefined(ln); ++ln)
    {
        gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, /*tag_buffer*/ 1);
    }
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    if (finest_ln != 1)
    {
        TBOX_ERROR("the hierarchy has " << finest_ln + 1 << " levels, not 2\n");
    }

    l_data_manager->registerLInitStrategy(new MarkerInitializer());
    l_data_manager->setPatchHierarchy(patch_hierarchy);
    l_data_manager->setPatchLevels(0, finest_ln);
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        l_data_manager->initializeLevelData(patch_hierarchy, ln, 0.0, ln < finest_ln, true);
    }
    l_data_manager->resetHierarchyConfiguration(patch_hierarchy, 0, finest_ln);
    l_data_manager->inactivateLagrangianStructures({ 1 }, 0);

    // Set the field, ghost cells included, and zero the spread values.
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(q_idx, 0.0);
        level->allocatePatchData(f_idx, 0.0);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> q_data = patch->getPatchData(q_idx);
            Pointer<CellData<NDIM, double>> f_data = patch->getPatchData(f_idx);
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            const double* const x_lower = pgeom->getXLower();
            const double* const dx = pgeom->getDx();
            for (CellIterator<NDIM> ic(q_data->getGhostBox()); ic; ic++)
            {
                double x[NDIM];
                for (int d = 0; d < NDIM; ++d) x[d] = x_lower[d] + dx[d] * (ic()(d) - patch->getBox().lower(d) + 0.5);
                (*q_data)(ic()) = linear_field(x);
            }
            f_data->fillAll(0.0);
        }
    }

    // Interpolate on each level. The interpolated value at the inactive marker is zero.
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<IBTK::LData> X_data = l_data_manager->getLData(IBTK::LDataManager::POSN_DATA_NAME, ln);
        Pointer<IBTK::LData> Q_data = l_data_manager->createLData("Q", ln, 1);
        l_data_manager->interp(q_idx, Q_data, X_data, ln);
        const std::vector<double> Q = gather_by_lagrangian_index(Q_data, l_data_manager, ln);
        for (std::size_t k = 0; k < Q.size(); ++k)
        {
            plog << "level " << ln << ", marker " << k << ": interpolated value = " << Q[k] << '\n';
        }
    }

    // Spread the value 1 from the markers of one level at a time. The integral on that level increases by the number of
    // its active markers.
    RobinPhysBdryPatchStrategy* const no_bdry_op = nullptr;
    for (int spread_ln = 0; spread_ln <= finest_ln; ++spread_ln)
    {
        Pointer<IBTK::LData> X_data = l_data_manager->getLData(IBTK::LDataManager::POSN_DATA_NAME, spread_ln);
        Pointer<IBTK::LData> F_data = l_data_manager->createLData("F", spread_ln, 1);
        {
            boost::multi_array_ref<double, 2>& F = *F_data->getLocalFormVecArray();
            std::fill_n(F.data(), F.num_elements(), 1.0);
            F_data->restoreArrays();
        }
        l_data_manager->spread(f_idx, F_data, X_data, no_bdry_op, spread_ln);
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            plog << "after spreading on level " << spread_ln << ": integral on level " << ln << " = "
                 << integrate(patch_hierarchy->getPatchLevel(ln), f_idx) << '\n';
        }
    }
}
