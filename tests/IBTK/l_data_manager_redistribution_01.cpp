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

#include <ibamr/IBTargetPointForceSpec.h>

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

#include <tbox/Pointer.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <StandardTagAndInitialize.h>

#include <array>
#include <sstream>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
constexpr int N_MARKERS = 4;

std::array<std::array<IBTK::Vector, N_MARKERS>, 4>
get_step_translations()
{
    const IBTK::Vector dx_pp(0.005, 0.005);
    const IBTK::Vector dx_mp(-0.005, 0.005);
    const IBTK::Vector dx_pm(0.005, -0.005);
    const IBTK::Vector dx_mm(-0.005, -0.005);
    return { std::array<IBTK::Vector, N_MARKERS>{ dx_pp, dx_mp, dx_pm, dx_mm },
             std::array<IBTK::Vector, N_MARKERS>{ dx_pp, dx_mp, dx_pm, dx_mm },
             std::array<IBTK::Vector, N_MARKERS>{ dx_pp, dx_mp, dx_pm, dx_mm },
             std::array<IBTK::Vector, N_MARKERS>{ dx_pp, dx_mp, dx_pm, dx_mm } };
}

std::array<IBTK::Point, N_MARKERS>
get_initial_positions()
{
    return { IBTK::Point(0.49, 0.49), IBTK::Point(0.51, 0.49), IBTK::Point(0.49, 0.51), IBTK::Point(0.51, 0.51) };
}

std::array<IBTK::Point, N_MARKERS>
get_exact_positions(const int step_number)
{
    auto positions = get_initial_positions();
    const auto step_translations = get_step_translations();
    for (int step = 0; step < step_number; ++step)
    {
        for (int k = 0; k < N_MARKERS; ++k)
        {
            positions[k] += step_translations[step][k];
        }
    }
    return positions;
}

class FourPointInitializer : public IBTK::LInitStrategy
{
public:
    FourPointInitializer()
    {
        IBAMR::IBTargetPointForceSpec::registerWithStreamableManager();
        d_X = get_initial_positions();
        d_U = get_step_translations()[0];
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
        return level_number == 0 ? N_MARKERS : 0;
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
        const IntVector<NDIM>& ratio = level->getRatio();
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            for (int k = 0; k < N_MARKERS; ++k)
            {
                const CellIndex<NDIM> cell_idx =
                    IndexUtilities::getAssignedCellIndex(d_X[k], hierarchy->getGridGeometry(), ratio);
                if (patch_box.contains(cell_idx)) ++count;
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
        for (int k = 0; k < N_MARKERS; ++k)
        {
            strct_id_to_strct_name_map[k] = "marker_" + std::to_string(k);
            strct_id_to_lag_idx_range_map[k] = std::make_pair(k, k + 1);
        }
    }

    unsigned int initializeDataOnPatchLevel(int lag_node_index_idx,
                                            unsigned int global_index_offset,
                                            unsigned int local_index_offset,
                                            Pointer<IBTK::LData> X_data,
                                            Pointer<IBTK::LData> U_data,
                                            Pointer<PatchHierarchy<NDIM>> hierarchy,
                                            int level_number,
                                            double /*init_data_time*/,
                                            bool /*can_be_refined*/,
                                            bool /*initial_time*/,
                                            IBTK::LDataManager* /*l_data_manager*/) override
    {
        if (level_number != 0) return 0;

        boost::multi_array_ref<double, 2>& X_array = *X_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& U_array = *U_data->getLocalFormVecArray();
        int local_idx = -1;
        unsigned int local_node_count = 0;

        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_number);
        const IntVector<NDIM>& ratio = level->getRatio();
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<IBTK::LNodeSetData> index_data = patch->getPatchData(lag_node_index_idx);
            const Box<NDIM>& patch_box = patch->getBox();
            for (int k = 0; k < N_MARKERS; ++k)
            {
                const CellIndex<NDIM> cell_idx =
                    IndexUtilities::getAssignedCellIndex(d_X[k], hierarchy->getGridGeometry(), ratio);
                if (!patch_box.contains(cell_idx)) continue;

                const int local_petsc_idx = ++local_idx + static_cast<int>(local_index_offset);
                const int global_petsc_idx = ++local_node_count - 1 + static_cast<int>(global_index_offset);

                for (int d = 0; d < NDIM; ++d)
                {
                    X_array[local_petsc_idx][d] = d_X[k][d];
                    U_array[local_petsc_idx][d] = d_U[k][d];
                }

                if (!index_data->isElement(cell_idx))
                {
                    index_data->appendItemPointer(cell_idx, new IBTK::LNodeSet());
                }
                IBTK::LNodeSet* const node_set = index_data->getItem(cell_idx);
                std::vector<Pointer<IBTK::Streamable>> node_data;
                node_data.push_back(new IBAMR::IBTargetPointForceSpec(k, 1.0, 0.0, d_X[k]));
                node_set->push_back(new IBTK::LNode(k,
                                                    global_petsc_idx,
                                                    local_petsc_idx,
                                                    IntVector<NDIM>(0),
                                                    IntVector<NDIM>(0),
                                                    IBTK::Vector::Zero(),
                                                    IBTK::Vector::Zero(),
                                                    node_data));
            }
        }

        local_node_count = static_cast<unsigned int>(local_idx + 1);

        X_data->restoreArrays();
        U_data->restoreArrays();
        return local_node_count;
    }

private:
    std::array<IBTK::Point, N_MARKERS> d_X;
    std::array<IBTK::Vector, N_MARKERS> d_U;
};

void
print_patch_data(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                 IBTK::LDataManager* l_data_manager,
                 Pointer<IBTK::LData> X_data,
                 Pointer<IBTK::LData> U_data,
                 const std::string& label,
                 const int step_number)
{
    const int rank = IBTK_MPI::getRank();
    if (rank != 0) return;

    std::ostringstream out;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<IBTK::LMesh> mesh = l_data_manager->getLMesh(finest_ln);
    const std::vector<IBTK::LNode*>& local_nodes = mesh->getLocalNodes();
    boost::multi_array_ref<double, 2>& X_array = *X_data->getGhostedLocalFormVecArray();
    boost::multi_array_ref<double, 2>& U_array = *U_data->getGhostedLocalFormVecArray();
    const auto exact_positions = get_exact_positions(step_number);

    out << label << " rank = 0 local_nodes = " << local_nodes.size() << '\n';
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(finest_ln);
    const IntVector<NDIM>& ratio = level->getRatio();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        out << "patch = " << patch->getBox() << '\n';
        int patch_count = 0;
        for (IBTK::LNode* node : local_nodes)
        {
            const int local_idx = node->getLocalPETScIndex();
            IBTK::Point X;
            for (int d = 0; d < NDIM; ++d) X[d] = X_array[local_idx][d];
            const CellIndex<NDIM> cell_idx =
                IndexUtilities::getAssignedCellIndex(X, patch_hierarchy->getGridGeometry(), ratio);
            if (!patch->getBox().contains(cell_idx)) continue;
            const auto* target_spec = node->getNodeDataItem<IBAMR::IBTargetPointForceSpec>();
            TBOX_ASSERT(target_spec);
            const IBTK::Point& X_target = target_spec->getTargetPointPosition();
            const IBTK::Point& X_target_exact = exact_positions[node->getLagrangianIndex()];

            out << "  index = " << node->getLagrangianIndex() << " X = " << X[0] << ", " << X[1]
                << " U = " << U_array[local_idx][0] << ", " << U_array[local_idx][1] << " X_target = " << X_target[0]
                << ", " << X_target[1] << " X_target_exact = " << X_target_exact[0] << ", " << X_target_exact[1]
                << '\n';
            ++patch_count;
        }
        out << "  count = " << patch_count << '\n';
    }

    X_data->restoreArrays();
    U_data->restoreArrays();
    plog << out.str();
}
} // namespace

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv);
    Logger::getInstance()->setWarning(false);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);

    if (IBTK_MPI::getNodes() != 4)
    {
        TBOX_ERROR("This test must be run with exactly 4 MPI processes.\n");
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

    gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);

    IBTK::LDataManager* const l_data_manager =
        IBTK::LDataManager::getManager("LDataManagerRedistribution01", "IB_4", "IB_4");
    l_data_manager->registerLInitStrategy(new FourPointInitializer());
    l_data_manager->setPatchHierarchy(patch_hierarchy);
    l_data_manager->setPatchLevels(0, 0);
    l_data_manager->initializeLevelData(patch_hierarchy, 0, 0.0, false, true);
    l_data_manager->resetHierarchyConfiguration(patch_hierarchy, 0, 0);

    Pointer<IBTK::LData> X_data = l_data_manager->getLData(IBTK::LDataManager::POSN_DATA_NAME, 0);
    Pointer<IBTK::LData> U_data = l_data_manager->getLData(IBTK::LDataManager::VEL_DATA_NAME, 0);
    print_patch_data(patch_hierarchy, l_data_manager, X_data, U_data, "before redistribution", 0);
    const auto step_translations = get_step_translations();
    for (std::size_t step = 0; step < step_translations.size(); ++step)
    {
        Pointer<IBTK::LMesh> mesh = l_data_manager->getLMesh(0);
        const std::vector<IBTK::LNode*>& local_nodes = mesh->getLocalNodes();
        boost::multi_array_ref<double, 2>& X_array = *X_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& U_array = *U_data->getLocalFormVecArray();
        for (IBTK::LNode* node : local_nodes)
        {
            const int lag_idx = node->getLagrangianIndex();
            const int local_idx = node->getLocalPETScIndex();
            for (int d = 0; d < NDIM; ++d)
            {
                const double dX = step_translations[step][lag_idx][d];
                U_array[local_idx][d] = dX;
                X_array[local_idx][d] += dX;
            }
            auto* target_spec = node->getNodeDataItem<IBAMR::IBTargetPointForceSpec>();
            TBOX_ASSERT(target_spec);
            IBTK::Point& X_target = target_spec->getTargetPointPosition();
            for (int d = 0; d < NDIM; ++d) X_target[d] = X_array[local_idx][d];
        }
        X_data->restoreArrays();
        U_data->restoreArrays();

        if (IBTK_MPI::getRank() == 0) plog << "step " << step + 1 << ": begin redistribution\n";
        l_data_manager->beginDataRedistribution(0, 0);
        l_data_manager->endDataRedistribution(0, 0);
        if (IBTK_MPI::getRank() == 0) plog << "step " << step + 1 << ": end redistribution\n";

        print_patch_data(
            patch_hierarchy, l_data_manager, X_data, U_data, "after step " + std::to_string(step + 1), step + 1);
    }
}
