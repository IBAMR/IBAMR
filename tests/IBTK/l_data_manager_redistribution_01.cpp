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
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/IBOperator.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LInitStrategy.h>
#include <ibtk/LMesh.h>
#include <ibtk/LNode.h>
#include <ibtk/LNodeSetData.h>
#include <ibtk/ib_kernel_evaluators.h>

#include <tbox/Pointer.h>
#include <tbox/RestartManager.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellVariable.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenOperator.h>
#include <CoarsenSchedule.h>
#include <EdgeData.h>
#include <EdgeGeometry.h>
#include <EdgeVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <NodeGeometry.h>
#include <NodeVariable.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideIterator.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
constexpr int N_MARKERS = 4;

IBTK::Point
make_point(const double x, const double y, const double z = 0.5)
{
    IBTK::Point result;
    result[0] = x;
    result[1] = y;
#if NDIM == 3
    result[2] = z;
#else
    (void)z;
#endif
    return result;
}

std::array<std::array<IBTK::Vector, N_MARKERS>, 4>
get_step_translations()
{
    const IBTK::Vector dx_pp = make_point(0.005, 0.005, 0.0);
    const IBTK::Vector dx_mp = make_point(-0.005, 0.005, 0.0);
    const IBTK::Vector dx_pm = make_point(0.005, -0.005, 0.0);
    const IBTK::Vector dx_mm = make_point(-0.005, -0.005, 0.0);
    return { std::array<IBTK::Vector, N_MARKERS>{ dx_pp, dx_mp, dx_pm, dx_mm },
             std::array<IBTK::Vector, N_MARKERS>{ dx_pp, dx_mp, dx_pm, dx_mm },
             std::array<IBTK::Vector, N_MARKERS>{ dx_pp, dx_mp, dx_pm, dx_mm },
             std::array<IBTK::Vector, N_MARKERS>{ dx_pp, dx_mp, dx_pm, dx_mm } };
}

std::array<IBTK::Point, N_MARKERS>
get_initial_positions()
{
    return { make_point(0.49, 0.49), make_point(0.51, 0.49), make_point(0.49, 0.51), make_point(0.51, 0.51) };
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
    explicit FourPointInitializer(const bool coupling = false, const bool two_levels = false)
        : d_finest_level(two_levels ? 1 : 0)
    {
        IBAMR::IBTargetPointForceSpec::registerWithStreamableManager();
        d_X = get_initial_positions();
        d_U = get_step_translations()[0];
        if (coupling)
        {
            d_X = { make_point(0.49, 0.01, 0.02),
                    make_point(0.51, 0.99, 0.98),
                    make_point(0.01, 0.49, 0.49),
                    make_point(0.99, 0.51, 0.51) };
        }
        if (two_levels)
        {
            d_X = {
                make_point(0.255, 0.49), make_point(0.265, 0.51), make_point(0.735, 0.49), make_point(0.745, 0.51)
            };
        }
    }

    bool getLevelHasLagrangianData(int level_number, bool /*can_be_refined*/) const override
    {
        return level_number <= d_finest_level;
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
        return level_number <= d_finest_level ? N_MARKERS : 0;
    }

    unsigned int computeLocalNodeCountOnPatchLevel(Pointer<PatchHierarchy<NDIM>> hierarchy,
                                                   int level_number,
                                                   double /*init_data_time*/,
                                                   bool /*can_be_refined*/,
                                                   bool /*initial_time*/) override
    {
        if (level_number > d_finest_level) return 0;

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
        if (level_number > d_finest_level) return;
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
        if (level_number > d_finest_level) return 0;

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
    const int d_finest_level;
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
// Reuse the initializer and native executable with two managers on real redistributed nodes.
int
check_manager_coupling(Pointer<AppInitializer> app)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 2);
    const bool physical = app->getInputDatabase()->getBoolWithDefault("physical_boundary", false);
    const std::string spread_override = app->getInputDatabase()->getStringWithDefault("spread_override", "");
    const bool wide_override = app->getInputDatabase()->getBoolWithDefault("wide_override", false);
    const std::string kernel = app->getInputDatabase()->getStringWithDefault("kernel", "BSPLINE_3");
    const int requested_ghost_width = app->getInputDatabase()->getIntegerWithDefault("requested_ghost_width", 3);
    const bool custom = app->getInputDatabase()->getBoolWithDefault("custom", true);
    const bool weighted = app->getInputDatabase()->getBoolWithDefault("weighted", false);
    const bool restart_test = app->getInputDatabase()->getBoolWithDefault("restart_test", false);
    const bool from_restart = RestartManager::getManager()->isFromRestart();
    const std::string spreading_kernel = app->getInputDatabase()->getStringWithDefault("spreading_kernel", kernel);
    Logger::getInstance()->setAbortAppender(new TestAppender());
    LDataManager* const direct = [&]
    {
        if (!custom)
        {
            return LDataManager::getManager("DirectManager",
                                            kernel,
                                            spreading_kernel,
                                            false,
                                            IntVector<NDIM>(requested_ghost_width),
                                            restart_test,
                                            true);
        }
        const IBOperator interpolation{ IBKernelTensorProduct(kernel) };
        const IBOperator spreading{ IBKernelTensorProduct(spreading_kernel) };
        return LDataManager::getManager(
            "DirectManager", interpolation, spreading, false, IntVector<NDIM>(requested_ghost_width), restart_test);
    }();
    LDataManager* const legacy =
        LDataManager::getManager("LegacyManager",
                                 kernel,
                                 spreading_kernel,
                                 false,
                                 IntVector<NDIM>(requested_ghost_width),
                                 restart_test,
                                 kernel == "COMPOSITE_BSPLINE_12" || spread_override == "COMPOSITE_BSPLINE_12");
    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = variables->getContext("coupling");
    Pointer<SideVariable<NDIM, double>> u_variable = new SideVariable<NDIM, double>("u");
    Pointer<SideVariable<NDIM, double>> f_variable = new SideVariable<NDIM, double>("f");
    Pointer<SideVariable<NDIM, double>> g_variable = new SideVariable<NDIM, double>("g");
    const bool missing_ghost = app->getInputDatabase()->getBoolWithDefault("missing_ghost", false);
    const int ghosts = missing_ghost ? 1 : direct->getGhostCellWidth().max();
    const int u = variables->registerVariableAndContext(u_variable, context, ghosts);
    const int f = variables->registerVariableAndContext(f_variable, context, ghosts);
    const int g = variables->registerVariableAndContext(g_variable, context, ghosts);
    Pointer<CartesianGridGeometry<NDIM>> geometry =
        new CartesianGridGeometry<NDIM>("CartesianGeometry", app->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", geometry);
    Pointer<StandardTagAndInitialize<NDIM>> detector =
        new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize", nullptr, Pointer<Database>(nullptr));
    Pointer<BergerRigoutsos<NDIM>> generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> grid = new GriddingAlgorithm<NDIM>(
        "GriddingAlgorithm", app->getComponentDatabase("GriddingAlgorithm"), detector, generator, balancer);
    Pointer<FourPointInitializer> initializer = new FourPointInitializer(true);
    if (from_restart)
    {
        hierarchy->getFromRestart(1);
    }
    else
    {
        grid->makeCoarsestLevel(hierarchy, 0.0);
    }
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    level->allocatePatchData(u);
    level->allocatePatchData(f);
    level->allocatePatchData(g);
    for (LDataManager* manager : { direct, legacy })
    {
        manager->registerLInitStrategy(initializer);
        manager->setPatchHierarchy(hierarchy);
        manager->setPatchLevels(0, 0);
        manager->initializeLevelData(hierarchy, 0, 0.0, false, !from_restart, nullptr, !from_restart);
        manager->resetHierarchyConfiguration(hierarchy, 0, 0);
        manager->beginDataRedistribution(0, 0);
        manager->endDataRedistribution(0, 0);
    }
    if (from_restart)
    {
        const std::array<IBTK::Point, N_MARKERS> initial{ make_point(0.49, 0.01, 0.02),
                                                          make_point(0.51, 0.99, 0.98),
                                                          make_point(0.01, 0.49, 0.49),
                                                          make_point(0.99, 0.51, 0.51) };
        for (LDataManager* manager : { direct, legacy })
        {
            Pointer<LData> X = manager->getLData(LDataManager::POSN_DATA_NAME, 0);
            const auto& positions = *X->getLocalFormVecArray();
            for (LNode* node : manager->getLMesh(0)->getLocalNodes())
            {
                for (int d = 0; d < NDIM; ++d)
                {
                    const double moved = initial[node->getLagrangianIndex()][d] + 0.03;
                    TBOX_ASSERT(std::abs(positions[node->getLocalPETScIndex()][d] - (moved - std::floor(moved))) <
                                1.0e-12);
                }
            }
            X->restoreArrays();
        }
        plog << "restored nodes = " << direct->getNumberOfNodes(0) << '\n';
    }
    if (wide_override)
    {
        Logger::getInstance()->setAbortAppender(new TestAppender());
        direct->spread(f,
                       direct->getLData(LDataManager::VEL_DATA_NAME, 0),
                       direct->getLData(LDataManager::POSN_DATA_NAME, 0),
                       app->getInputDatabase()->getStringWithDefault("wide_override_kernel", "BSPLINE_6"),
                       nullptr,
                       0);
        return 0;
    }
    if (missing_ghost)
    {
        Logger::getInstance()->setAbortAppender(new TestAppender());
        direct->interp(
            u, direct->getLData(LDataManager::VEL_DATA_NAME, 0), direct->getLData(LDataManager::POSN_DATA_NAME, 0), 0);
        return 0;
    }
    TBOX_ASSERT(IBTK_MPI::sumReduction(static_cast<int>(direct->getLMesh(0)->getGhostNodes().size())) > 0);
    plog << std::scientific << std::setprecision(12);
    LocationIndexRobinBcCoefs<NDIM> coefficients;
    for (int side = 0; side < 2 * NDIM; ++side)
    {
        coefficients.setBoundarySlope(side, 0.0);
    }
    const std::vector<RobinBcCoefStrategy<NDIM>*> boundary_coefficients(NDIM, &coefficients);
    CartSideRobinPhysBdryOp boundary(f, boundary_coefficients, true);
    for (int step = from_restart ? 1 : 0; step < 2; ++step)
    {
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pg = patch->getPatchGeometry();
            Pointer<SideData<NDIM, double>> velocity = patch->getPatchData(u);
            Pointer<SideData<NDIM, double>> force = patch->getPatchData(f), other = patch->getPatchData(g);
            force->fillAll(0.125);
            other->fillAll(0.125);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> side(velocity->getGhostBox(), axis); side; side++)
                {
                    double value = 0.3 + axis;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        const double x =
                            pg->getXLower()[d] +
                            (side()(d) - patch->getBox().lower()(d) + (axis == d ? 0.0 : 0.5)) * pg->getDx()[d];
                        value += (d + 1) * std::sin(2.0 * 3.141592653589793 * x + 0.37);
                    }
                    (*velocity)(side()) = value;
                }
            }
            if (physical)
            {
                boundary.setPatchDataIndex(u);
                boundary.setPhysicalBoundaryConditions(*patch, 0.0, IntVector<NDIM>(3));
            }
        }
        std::array<std::array<double, NDIM>, N_MARKERS> gathered_direct{}, gathered_legacy{};
        for (LDataManager* manager : { direct, legacy })
        {
            Pointer<LData> X = manager->getLData(LDataManager::POSN_DATA_NAME, 0);
            Pointer<LData> U = manager->getLData(LDataManager::VEL_DATA_NAME, 0);
            manager->interp(u, U, X, 0);
            Pointer<LData> ds = manager->createLData("weights", 0, 1, false);
            boost::multi_array_ref<double, 1>& weights = *ds->getLocalFormArray();
            boost::multi_array_ref<double, 2>& values = *U->getLocalFormVecArray();
            for (LNode* node : manager->getLMesh(0)->getLocalNodes())
            {
                const int local = node->getLocalPETScIndex(), lag = node->getLagrangianIndex();
                weights[local] = 0.3 + 0.1 * lag;
                for (int d = 0; d < NDIM; ++d)
                {
                    TBOX_ASSERT(std::isfinite(values[local][d]));
                    (manager == direct ? gathered_direct : gathered_legacy)[lag][d] = values[local][d];
                    values[local][d] =
                        (0.2 * (lag + 1) + 0.1 * d) * (weighted && manager == legacy ? weights[local] : 1.0);
                }
            }
            U->restoreArrays();
            ds->restoreArrays();
            if (weighted && manager == direct)
            {
                if (spread_override.empty())
                {
                    manager->spread(f, U, X, ds, physical ? &boundary : nullptr, 0);
                }
                else
                {
                    manager->spread(f, U, X, ds, spread_override, physical ? &boundary : nullptr, 0);
                }
            }
            else if (spread_override.empty())
            {
                manager->spread(manager == direct ? f : g, U, X, physical ? &boundary : nullptr, 0);
            }
            else
            {
                manager->spread(manager == direct ? f : g, U, X, spread_override, physical ? &boundary : nullptr, 0);
            }
        }
        double gather_error = 0.0, spread_error = 0.0, mass = 0.0, pairing = 0.0, expected_mass = 0.0;
        for (int k = 0; k < N_MARKERS; ++k)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                gather_error = std::max(gather_error, std::abs(gathered_direct[k][d] - gathered_legacy[k][d]));
                const double force = (0.2 * (k + 1) + 0.1 * d) * (weighted ? 0.3 + 0.1 * k : 1.0);
                expected_mass += force;
                pairing += gathered_direct[k][d] * force;
            }
        }
        double grid_pairing = 0.0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pg = patch->getPatchGeometry();
            Pointer<SideData<NDIM, double>> velocity = patch->getPatchData(u), force = patch->getPatchData(f),
                                            other = patch->getPatchData(g);
            double volume = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                volume *= pg->getDx()[d];
            }
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> side(force->getGhostBox(), axis); side; side++)
                {
                    TBOX_ASSERT(std::isfinite((*force)(side())) && std::isfinite((*other)(side())));
                    spread_error = std::max(spread_error, volume * std::abs((*force)(side()) - (*other)(side())));
                }
                // Count each periodic side once by excluding the upper normal face of each patch.
                for (SideIterator<NDIM> side(patch->getBox(), axis); side; side++)
                {
                    if (side()(axis) > patch->getBox().upper()(axis))
                    {
                        continue;
                    }
                    const double contribution = volume * ((*force)(side()) - 0.125);
                    mass += contribution;
                    grid_pairing += contribution * (*velocity)(side());
                }
            }
        }
        gather_error = IBTK_MPI::maxReduction(gather_error);
        spread_error = IBTK_MPI::maxReduction(spread_error);
        mass = IBTK_MPI::sumReduction(mass);
        const double adjoint = std::abs(IBTK_MPI::sumReduction(grid_pairing - pairing));
        plog << "step " << step << " manager_errors " << gather_error << ' ' << spread_error;
        if (!physical)
        {
            plog << " grid_sum " << mass << " pairing_difference " << adjoint;
        }
        plog << '\n';
        TBOX_ASSERT(gather_error < 1.0e-11 && spread_error < 1.0e-11);
        if (!physical)
        {
            TBOX_ASSERT(std::abs(mass - expected_mass) < 1.0e-11 &&
                        (!spread_override.empty() || spreading_kernel != kernel || adjoint < 1.0e-11));
        }
        for (LDataManager* manager : { direct, legacy })
        {
            Pointer<LData> X = manager->getLData(LDataManager::POSN_DATA_NAME, 0);
            boost::multi_array_ref<double, 2>& coordinates = *X->getLocalFormVecArray();
            for (LNode* node : manager->getLMesh(0)->getLocalNodes())
            {
                const int local = node->getLocalPETScIndex();
                for (int d = 0; d < NDIM; ++d)
                {
                    coordinates[local][d] += physical ? (coordinates[local][d] < 0.5 ? 0.03 : -0.03) : 0.03;
                }
            }
            X->restoreArrays();
            manager->beginDataRedistribution(0, 0);
            manager->endDataRedistribution(0, 0);
        }
        if (restart_test && step == 0)
        {
            RestartManager::getManager()->writeRestartFile("restart", 1);
        }
    }
    return 0;
}
// Exercise component layouts and real hierarchy schedules with independent managers.
template <DataCentering C>
int
check_component_manager(Pointer<AppInitializer> app)
{
    static_assert(C == DataCentering::CELL || C == DataCentering::NODE || C == DataCentering::EDGE);
    using Data = typename CartesianCentering<C>::template Data<double>;
    TBOX_ASSERT(IBTK_MPI::getNodes() == 2);
    Pointer<Database> input = app->getInputDatabase();
    const bool amr = input->getBoolWithDefault("two_levels", false);
    const bool custom = input->getBoolWithDefault("custom", true);
    const bool missing_ghost = input->getBoolWithDefault("missing_ghost", false);
    const bool mismatch = input->getBoolWithDefault("mismatched_depth", false);
    const int depth = C == DataCentering::EDGE ? NDIM : input->getIntegerWithDefault("component_depth", 1);
    const bool wrong_positions = input->getBoolWithDefault("wrong_positions", false);
    const bool wrong_field = input->getBoolWithDefault("wrong_field", false);
    const bool dimension_error = input->getBoolWithDefault("dimension_error", false);
    const bool wide_override = input->getBoolWithDefault("wide_override", false);
    const int registered_ghosts = input->getIntegerWithDefault("registered_ghosts", 5);
    const std::string kernel = input->getStringWithDefault("kernel", "COMPOSITE_BSPLINE_23");
    const std::string spreading_kernel = input->getStringWithDefault("spreading_kernel", kernel);
    const std::string override_kernel = input->getStringWithDefault("spread_override", "");
    Logger::getInstance()->setAbortAppender(new TestAppender());
    LDataManager* const direct = [&]
    {
        if (custom)
        {
            const IBOperator interpolation =
                C == DataCentering::EDGE && kernel == "COMPOSITE_BSPLINE_56" ?
                    IBOperator{ IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::BSpline<5>{},
                                                                IBKernelEvaluators::BSpline<6>{} } } :
                    IBOperator{ IBKernelTensorProduct(kernel) };
            const IBOperator spreading =
                spreading_kernel == kernel ? interpolation : IBOperator{ IBKernelTensorProduct(spreading_kernel) };
            return LDataManager::getManager(
                "ComponentDirect", interpolation, spreading, false, IntVector<NDIM>(registered_ghosts), false);
        }
        return LDataManager::getManager(
            "ComponentDirect", kernel, spreading_kernel, false, IntVector<NDIM>(registered_ghosts), false, true);
    }();
    LDataManager* const legacy = LDataManager::getManager(
        "ComponentLegacy", kernel, spreading_kernel, false, IntVector<NDIM>(registered_ghosts), false);
    const auto entry = [](Data& field, const hier::Index<NDIM>& i, const int c) -> double&
    {
        if constexpr (C == DataCentering::EDGE)
        {
            return field.getArrayData(c)(i, 0);
        }
        else
        {
            return field.getArrayData()(i, c);
        }
    };
    const auto interior_box = [](const Box<NDIM>& box, const int c)
    {
        if constexpr (C == DataCentering::EDGE)
        {
            return EdgeGeometry<NDIM>::toEdgeBox(box, c);
        }
        else
        {
            return C == DataCentering::NODE ? NodeGeometry<NDIM>::toNodeBox(box) : box;
        }
    };
    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<hier::Variable<NDIM>> variable;
    if constexpr (C == DataCentering::CELL)
    {
        variable = new CellVariable<NDIM, double>("component_field", depth);
    }
    else if constexpr (C == DataCentering::EDGE)
    {
        variable = new EdgeVariable<NDIM, double>("edge_field", wrong_field ? 2 : 1);
    }
    else
    {
        variable = new NodeVariable<NDIM, double>("component_field", depth);
    }
    const int ghosts = missing_ghost ? 1 : direct->getGhostCellWidth().max();
    const int u = variables->registerVariableAndContext(variable, variables->getContext("component_u"), ghosts);
    const int f = variables->registerVariableAndContext(variable, variables->getContext("component_f"), ghosts);
    const int g = variables->registerVariableAndContext(variable, variables->getContext("component_g"), ghosts);
    Pointer<CartesianGridGeometry<NDIM>> geometry =
        new CartesianGridGeometry<NDIM>("CartesianGeometry", app->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", geometry);
    Pointer<StandardTagAndInitialize<NDIM>> detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize",
        nullptr,
        amr ? app->getComponentDatabase("StandardTagAndInitialize") : Pointer<Database>(nullptr));
    Pointer<BergerRigoutsos<NDIM>> generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> grid = new GriddingAlgorithm<NDIM>(
        "GriddingAlgorithm", app->getComponentDatabase("GriddingAlgorithm"), detector, generator, balancer);
    grid->makeCoarsestLevel(hierarchy, 0.0);
    if (amr)
    {
        grid->makeFinerLevel(hierarchy, 0.0, 0.0, 0);
    }
    const int finest = hierarchy->getFinestLevelNumber();
    TBOX_ASSERT(finest == (amr ? 1 : 0));
    Pointer<FourPointInitializer> initializer = new FourPointInitializer(true, amr);
    for (LDataManager* manager : { direct, legacy })
    {
        manager->registerLInitStrategy(initializer);
        manager->setPatchHierarchy(hierarchy);
        manager->setPatchLevels(0, finest);
        for (int ln = 0; ln <= finest; ++ln)
        {
            manager->initializeLevelData(hierarchy, ln, 0.0, ln < finest, true);
        }
        manager->resetHierarchyConfiguration(hierarchy, 0, finest);
        manager->beginDataRedistribution(0, finest);
        manager->endDataRedistribution(0, finest);
        for (int ln = 0; ln <= finest; ++ln)
        {
            TBOX_ASSERT(manager->getNumberOfNodes(ln) == N_MARKERS);
        }
    }
    CartExtrapPhysBdryOp boundary(u, "LINEAR");
    RefineAlgorithm<NDIM> fill, prolong_f, prolong_g;
    CoarsenAlgorithm<NDIM> synchronize;
    Pointer<RefineOperator<NDIM>> refine_op = geometry->lookupRefineOperator(
        variable, C == DataCentering::EDGE ? "CONSERVATIVE_LINEAR_REFINE" : "LINEAR_REFINE");
    Pointer<CoarsenOperator<NDIM>> coarsen_op = geometry->lookupCoarsenOperator(
        variable, C == DataCentering::NODE ? "CONSTANT_COARSEN" : "CONSERVATIVE_COARSEN");
    TBOX_ASSERT(refine_op && coarsen_op);
    fill.registerRefine(u, u, u, refine_op);
    prolong_f.registerRefine(f, f, f, refine_op);
    prolong_g.registerRefine(g, g, g, refine_op);
    synchronize.registerCoarsen(u, u, coarsen_op);
    std::vector<Pointer<RefineSchedule<NDIM>>> fill_schedules(finest + 1), f_schedules(finest + 1),
        g_schedules(finest + 1);
    std::vector<Pointer<CoarsenSchedule<NDIM>>> sync_schedules(finest + 1);
    for (int ln = 0; ln <= finest; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (int idx : { u, f, g })
        {
            level->allocatePatchData(idx);
        }
        fill_schedules[ln] =
            fill.createSchedule(level, ln - 1, hierarchy, C == DataCentering::EDGE ? nullptr : &boundary);
        if (ln > 0)
        {
            sync_schedules[ln] = synchronize.createSchedule(hierarchy->getPatchLevel(ln - 1), level);
            f_schedules[ln] = prolong_f.createSchedule(level, nullptr, ln - 1, hierarchy);
            g_schedules[ln] = prolong_g.createSchedule(level, nullptr, ln - 1, hierarchy);
        }
    }
    if (amr)
    {
        const Pointer<PatchLevel<NDIM>> fine_level = hierarchy->getPatchLevel(1);
        const BoxArray<NDIM>& fine_boxes = fine_level->getBoxes();
        const Box<NDIM> expected_region(hier::Index<NDIM>(8), hier::Index<NDIM>(23));
        int cells = 0;
        for (int k = 0; k < fine_boxes.size(); ++k)
        {
            TBOX_ASSERT((fine_boxes[k] * expected_region) == fine_boxes[k]);
            cells += fine_boxes[k].size();
        }
        TBOX_ASSERT(cells == expected_region.size());
        TBOX_ASSERT(cells > 0 &&
                    cells < Box<NDIM>::refine(geometry->getPhysicalDomain()[0], IntVector<NDIM>(2)).size());
        const boost::multi_array_ref<double, 2>& X =
            *direct->getLData(LDataManager::POSN_DATA_NAME, 1)->getLocalFormVecArray();
        int near_interface = 0;
        for (LNode* node : direct->getLMesh(1)->getLocalNodes())
        {
            const double x = X[node->getLocalPETScIndex()][0];
            if (std::min(std::abs(x - 0.25), std::abs(x - 0.75)) < 1.0 / 32.0)
            {
                ++near_interface;
            }
        }
        direct->getLData(LDataManager::POSN_DATA_NAME, 1)->restoreArrays();
        near_interface = IBTK_MPI::sumReduction(near_interface);
        TBOX_ASSERT(near_interface == N_MARKERS);
        plog << "fine_cells " << cells << " interface_markers " << near_interface << '\n';
    }
    if constexpr (C == DataCentering::EDGE)
    {
        if (amr)
        {
            for (int ln = 0; ln <= finest; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Data> field = level->getPatch(p())->getPatchData(f);
                    field->fillAll(ln == 0 ? 2.0 : 0.0);
                }
            }
            f_schedules[1]->fillData(0.0);
            double prolongation = 0.0;
            Pointer<PatchLevel<NDIM>> fine = hierarchy->getPatchLevel(1);
            for (PatchLevel<NDIM>::Iterator p(fine); p; p++)
            {
                Pointer<Patch<NDIM>> patch = fine->getPatch(p());
                Pointer<Data> field = patch->getPatchData(f);
                for (int c = 0; c < depth; ++c)
                {
                    for (Box<NDIM>::Iterator i(interior_box(patch->getBox(), c)); i; i++)
                    {
                        const double value = entry(*field, i(), c);
                        if (!std::isfinite(value))
                        {
                            TBOX_ERROR("Nonfinite edge prolongation.\n");
                        }
                        prolongation = std::max(prolongation, std::abs(value));
                    }
                }
            }
            prolongation = IBTK_MPI::maxReduction(prolongation);
            TBOX_ASSERT(std::abs(prolongation - 2.0) < 1.0e-12);
            plog << "edge_prolongation " << prolongation << '\n';
        }
    }
    plog << std::scientific << std::setprecision(12);
    for (int step = 0; step < 2; ++step)
    {
        std::vector<std::vector<double>> gathered[2];
        gathered[0].resize(finest + 1, std::vector<double>(N_MARKERS * depth));
        gathered[1].resize(finest + 1, std::vector<double>(N_MARKERS * depth));
        for (int ln = 0; ln <= finest; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Data> force = level->getPatch(p())->getPatchData(f),
                              other = level->getPatch(p())->getPatchData(g);
                force->fillAll(0.125);
                other->fillAll(0.125);
            }
        }
        int m = 0;
        for (LDataManager* manager : { direct, legacy })
        {
            std::vector<Pointer<LData>> X(finest + 1), Q(finest + 1), weights(finest + 1);
            for (int ln = 0; ln <= finest; ++ln)
            {
                X[ln] = wrong_positions ? manager->createLData("bad_positions", ln, NDIM + 1, false) :
                                          manager->getLData(LDataManager::POSN_DATA_NAME, ln);
                Q[ln] = manager->createLData("components", ln, depth + (mismatch ? 1 : 0), false);
                weights[ln] = manager->createLData("quadrature", ln, 1, false);
                IBTK_CHKERRQ(VecSet(Q[ln]->getVec(), -123.0));
                // Restore independently: interpolation synchronizes coarse data from the fine level.
                Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<CartesianPatchGeometry<NDIM>> pg = patch->getPatchGeometry();
                    Pointer<Data> field = patch->getPatchData(u);
                    field->fillAll(std::numeric_limits<double>::quiet_NaN());
                    for (int c = 0; c < depth; ++c)
                    {
                        for (Box<NDIM>::Iterator i(interior_box(patch->getBox(), c)); i; i++)
                        {
                            double value = 0.3 + c + 0.2 * ln;
                            for (int d = 0; d < NDIM; ++d)
                            {
                                const double x =
                                    pg->getXLower()[d] +
                                    (i()(d) - patch->getBox().lower()(d) +
                                     (C == DataCentering::CELL || (C == DataCentering::EDGE && d == c) ? 0.5 : 0.0)) *
                                        pg->getDx()[d];
                                value += (d + 1) * std::sin(2.0 * 3.141592653589793 * x + 0.37);
                            }
                            entry(*field, i(), c) = value;
                        }
                    }
                }
            }
            manager->interp(u, Q, X, sync_schedules, fill_schedules, 0.0);
            if (amr)
            {
                // A fine-level offset must reach the covered coarse node through synchronization.
                const hier::Index<NDIM> center(8);
                double synchronization_change = 0.0, initial_value = 0.3;
                for (int d = 0; d < NDIM; ++d)
                {
                    initial_value +=
                        (d + 1) * std::sin(2.0 * 3.141592653589793 *
                                               (0.5 + (C == DataCentering::EDGE && d == 0 ? 1.0 / 32.0 : 0.0)) +
                                           0.37);
                }
                const Pointer<PatchLevel<NDIM>> coarse_level = hierarchy->getPatchLevel(0);
                for (PatchLevel<NDIM>::Iterator p(coarse_level); p; p++)
                {
                    const Pointer<Patch<NDIM>> patch = coarse_level->getPatch(p());
                    if (patch->getBox().contains(center))
                    {
                        const Pointer<Data> field = patch->getPatchData(u);
                        const double value = entry(*field, center, 0);
                        if (!std::isfinite(value))
                        {
                            TBOX_ERROR("Nonfinite synchronized manager field.\n");
                        }
                        synchronization_change = std::max(synchronization_change, std::abs(value - initial_value));
                    }
                }
                synchronization_change = IBTK_MPI::maxReduction(synchronization_change);
                TBOX_ASSERT(synchronization_change > 0.1);
                if (m == 0)
                {
                    plog << "step " << step << " coarse_sync_change " << synchronization_change << '\n';
                }
            }
            if (mismatch || missing_ghost || wrong_positions || wrong_field || dimension_error)
            {
                return 0;
            }
            for (int ln = 0; ln <= finest; ++ln)
            {
                boost::multi_array_ref<double, 2>& values = *Q[ln]->getLocalFormVecArray();
                boost::multi_array_ref<double, 1>& ds = *weights[ln]->getLocalFormArray();
                for (LNode* node : manager->getLMesh(ln)->getLocalNodes())
                {
                    const int local = node->getLocalPETScIndex(), lag = node->getLagrangianIndex();
                    ds[local] = 0.3 + 0.1 * lag;
                    for (int c = 0; c < depth; ++c)
                    {
                        if (!std::isfinite(values[local][c]) || values[local][c] == -123.0)
                        {
                            TBOX_ERROR("Manager gather did not overwrite a finite value.\n");
                        }
                        gathered[m][ln][lag * depth + c] = values[local][c];
                        values[local][c] = (0.2 * (lag + 1) + 0.1 * c) * (m == 1 ? ds[local] : 1.0);
                    }
                }
                Q[ln]->restoreArrays();
                weights[ln]->restoreArrays();
            }
            if (m == 0)
            {
                if (override_kernel.empty())
                {
                    manager->spread(f, Q, X, weights, nullptr, f_schedules, 0.0);
                }
                else
                {
                    manager->spread(f, Q, X, weights, override_kernel, nullptr, f_schedules, 0.0);
                }
            }
            else
            {
                if (override_kernel.empty())
                {
                    manager->spread(g, Q, X, nullptr, g_schedules, 0.0);
                }
                else
                {
                    manager->spread(g, Q, X, override_kernel, nullptr, g_schedules, 0.0);
                }
            }
            if (wide_override)
            {
                return 0;
            }
            ++m;
        }
        if constexpr (C == DataCentering::EDGE)
        {
            for (int ln = 0; ln <= finest; ++ln)
            {
                Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                const Box<NDIM> domain = Box<NDIM>::refine(geometry->getPhysicalDomain()[0], level->getRatio());
                const int n = domain.numberCells()(0), count = NDIM * domain.size();
                std::vector<double> lower(count, std::numeric_limits<double>::infinity()),
                    upper(count, -std::numeric_limits<double>::infinity());
                std::vector<int> copies(count, 0);
                for (PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(p());
                    Pointer<Data> force = patch->getPatchData(f);
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        for (Box<NDIM>::Iterator i(interior_box(patch->getBox(), axis)); i; i++)
                        {
                            int slot = axis;
                            for (int d = NDIM - 1; d >= 0; --d)
                            {
                                TBOX_ASSERT(domain.numberCells()(d) == n);
                                slot = slot * n + (i()(d) - domain.lower()(d)) % n;
                            }
                            const double value = entry(*force, i(), axis);
                            if (!std::isfinite(value))
                            {
                                TBOX_ERROR("Nonfinite shared edge copy.\n");
                            }
                            lower[slot] = std::min(lower[slot], value);
                            upper[slot] = std::max(upper[slot], value);
                            ++copies[slot];
                        }
                    }
                }
                IBTK_MPI::minReduction(lower.data(), count);
                IBTK_MPI::maxReduction(upper.data(), count);
                IBTK_MPI::sumReduction(copies.data(), count);
                double copy_error = 0.0;
                int shared = 0;
                for (int k = 0; k < count; ++k)
                {
                    if (copies[k] > 1)
                    {
                        ++shared;
                        copy_error = std::max(copy_error, std::abs(upper[k] - lower[k]));
                    }
                }
                TBOX_ASSERT(shared > 0 && copy_error < 1.0e-10);
                plog << "step " << step << " level " << ln << " shared_edges " << shared << " copy_error " << copy_error
                     << '\n';
            }
        }
        for (int ln = 0; ln <= finest; ++ln)
        {
            double gather_error = 0.0, spread_error = 0.0, mass = 0.0, grid_pairing = 0.0, pairing = 0.0,
                   expected_mass = 0.0;
            for (int k = 0; k < N_MARKERS; ++k)
            {
                for (int c = 0; c < depth; ++c)
                {
                    const int idx = k * depth + c;
                    gather_error = std::max(gather_error, std::abs(gathered[0][ln][idx] - gathered[1][ln][idx]));
                    const double force = (0.2 * (k + 1) + 0.1 * c) * (0.3 + 0.1 * k);
                    expected_mass += force;
                    pairing += force * gathered[0][ln][idx];
                }
            }
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CartesianPatchGeometry<NDIM>> pg = patch->getPatchGeometry();
                Pointer<Data> velocity = patch->getPatchData(u), force = patch->getPatchData(f),
                              other = patch->getPatchData(g);
                double volume = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    volume *= pg->getDx()[d];
                }
                for (int c = 0; c < depth; ++c)
                {
                    for (Box<NDIM>::Iterator i(interior_box(force->getGhostBox(), c)); i; i++)
                    {
                        const double a = entry(*force, i(), c), b = entry(*other, i(), c);
                        if (!std::isfinite(a) || !std::isfinite(b))
                        {
                            TBOX_ERROR("Nonfinite manager spread comparison.\n");
                        }
                        spread_error = std::max(spread_error, volume * std::abs(a - b));
                    }
                }
                // The cell range keeps all tangent edges and excludes upper transverse copies.
                // It also counts shared/periodic nodes once for NODE data.
                for (Box<NDIM>::Iterator i(patch->getBox()); i; i++)
                {
                    for (int c = 0; c < depth; ++c)
                    {
                        const double contribution = volume * (entry(*force, i(), c) - 0.125);
                        mass += contribution;
                        grid_pairing += contribution * entry(*velocity, i(), c);
                    }
                }
            }
            if (!std::isfinite(mass) || !std::isfinite(grid_pairing) || !std::isfinite(pairing))
            {
                TBOX_ERROR("Nonfinite manager diagnostic.\n");
            }
            gather_error = IBTK_MPI::maxReduction(gather_error);
            spread_error = IBTK_MPI::maxReduction(spread_error);
            mass = IBTK_MPI::sumReduction(mass);
            const double adjoint = std::abs(IBTK_MPI::sumReduction(grid_pairing - pairing));
            plog << "step " << step << " level " << ln << " depth " << depth << " errors " << gather_error << ' '
                 << spread_error;
            if (!amr)
            {
                plog << " mass " << mass << " pairing_difference " << adjoint;
            }
            plog << '\n';
            TBOX_ASSERT(gather_error < 1.0e-11 && spread_error < 1.0e-11);
            if (!amr)
            {
                TBOX_ASSERT(
                    std::abs(mass - expected_mass) < 1.0e-11 &&
                    ((override_kernel.empty() ? spreading_kernel : override_kernel) != kernel || adjoint < 1.0e-11));
            }
        }
        for (LDataManager* manager : { direct, legacy })
        {
            for (int ln = 0; ln <= finest; ++ln)
            {
                Pointer<LData> X = manager->getLData(LDataManager::POSN_DATA_NAME, ln);
                boost::multi_array_ref<double, 2>& positions = *X->getLocalFormVecArray();
                for (LNode* node : manager->getLMesh(ln)->getLocalNodes())
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        double& x = positions[node->getLocalPETScIndex()][d];
                        x += amr ? (x < 0.5 ? 0.01 : -0.01) : 0.03;
                    }
                }
                X->restoreArrays();
            }
            manager->beginDataRedistribution(0, finest);
            manager->endDataRedistribution(0, finest);
        }
    }
    return 0;
}
} // namespace

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv);
    Logger::getInstance()->setWarning(false);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);

    const std::string centering = app_initializer->getInputDatabase()->getStringWithDefault("component_centering", "");
    if (centering == "CELL")
    {
        return check_component_manager<DataCentering::CELL>(app_initializer);
    }
    if (centering == "EDGE")
    {
        return check_component_manager<DataCentering::EDGE>(app_initializer);
    }
    if (centering == "NODE")
    {
        return check_component_manager<DataCentering::NODE>(app_initializer);
    }

    if (app_initializer->getInputDatabase()->getBoolWithDefault("matrix_free_coupling", false))
    {
        return check_manager_coupling(app_initializer);
    }

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
