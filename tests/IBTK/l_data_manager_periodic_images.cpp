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

// Test LDataManager::interp() and spread() at points that have crossed a periodic boundary since the last
// redistribution. The first set of points is moved by up to 0.9 grid cells. The second set contains the moved points
// translated by one period into the domain. For every kernel, the program checks that the weights of each point sum to
// 1, that interpolation of a smooth periodic field gives equal values at the two points of each pair, and that
// spreading the forces of either set gives the same grid function, whose integral is the sum of the forces.

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
#include <CellData.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
// Motions of a point along one axis, in grid cells from the lower end of the domain: stay inside (1), cross a boundary
// (4), move near the lower boundary (1), end 6e-18 beyond it or on the upper one (2), or lie one ulp below a face (1).
struct Motion
{
    double x0;
    double x1;
    bool below_face = false;
};

constexpr int NUM_MOTIONS = 9;

std::array<Motion, NUM_MOTIONS>
get_motions(const int N, const int axis)
{
    return { Motion{ 0.5 * N + 0.41 + 0.07 * axis, 0.5 * N + 0.81 + 0.07 * axis },
             Motion{ N - 0.30, N + 0.60 },
             Motion{ N - 0.45, N + 0.25 },
             Motion{ 0.30, -0.60 },
             Motion{ 0.45, -0.25 },
             Motion{ 2.30, 1.90 },
             Motion{ 0.30, -1.0e-16 },
             Motion{ N - 0.30, static_cast<double>(N) },
             Motion{ 0.5 * N, 0.5 * N, true } };
}

// The positions of the points, with N cells along each axis of the domain. The first set is stored at "stored" and
// moved to "moved". The second set is at "image", the moved positions translated by one period into the domain.
struct Points
{
    std::vector<IBTK::Point> stored;
    std::vector<IBTK::Point> moved;
    std::vector<IBTK::Point> image;
    int num_crossing = 0;
};

Points
build_points(const int N, const double* const x_lo, const double* const x_up)
{
    Points points;
    const double h = 1.0 / N;
    int num_points = 1;
    for (int d = 0; d < NDIM; ++d) num_points *= NUM_MOTIONS;
    for (int k = 0; k < num_points; ++k)
    {
        IBTK::Point X0, X1, X2;
        int remaining = k;
        bool crossing = false;
        for (int d = 0; d < NDIM; ++d)
        {
            const Motion m = get_motions(N, d)[remaining % NUM_MOTIONS];
            remaining /= NUM_MOTIONS;
            const double length = x_up[d] - x_lo[d];
            X0[d] = x_lo[d] + m.x0 * h * length;
            X1[d] = x_lo[d] + m.x1 * h * length;
            if (m.below_face) X0[d] = X1[d] = std::nextafter(x_lo[d] + 0.5 * length, x_lo[d]);
            X2[d] = X1[d] < x_lo[d] ? std::min(X1[d] + length, std::nextafter(x_up[d], x_lo[d])) :
                                      (X1[d] >= x_up[d] ? X1[d] - length : X1[d]);
            crossing = crossing || X1[d] < x_lo[d] || X1[d] >= x_up[d];
        }
        points.stored.push_back(X0);
        points.moved.push_back(X1);
        points.image.push_back(X2);
        if (crossing) ++points.num_crossing;
    }
    return points;
}

// Initialize the points of the first set at their stored positions and the points of the second set at their images.
class PointInitializer : public IBTK::LInitStrategy
{
public:
    PointInitializer(const Points& points)
    {
        d_X = points.stored;
        d_X.insert(d_X.end(), points.image.begin(), points.image.end());
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

// Move the first set to its moved positions without assigning the points to patches again, as when points move
// between redistributions.
void
move_points(Pointer<IBTK::LData> X_data, IBTK::LDataManager* l_data_manager, const Points& points)
{
    boost::multi_array_ref<double, 2>& X = *X_data->getLocalFormVecArray();
    for (const IBTK::LNode* node : l_data_manager->getLMesh(0)->getLocalNodes())
    {
        const int lag_idx = node->getLagrangianIndex();
        if (lag_idx >= static_cast<int>(points.moved.size())) continue;
        for (int d = 0; d < NDIM; ++d) X[node->getLocalPETScIndex()][d] = points.moved[lag_idx][d];
    }
    X_data->restoreArrays();
}

// Return the entries of a quantity with one value per point, indexed by the Lagrangian index.
std::vector<double>
gather_by_lagrangian_index(Pointer<IBTK::LData> Q_data,
                           IBTK::LDataManager* l_data_manager,
                           const std::size_t num_points)
{
    std::vector<double> Q(num_points, 0.0);
    boost::multi_array_ref<double, 2>& Q_array = *Q_data->getLocalFormVecArray();
    for (const IBTK::LNode* node : l_data_manager->getLMesh(0)->getLocalNodes())
    {
        Q[node->getLagrangianIndex()] = Q_array[node->getLocalPETScIndex()][0];
    }
    Q_data->restoreArrays();
    IBTK_MPI::sumReduction(Q.data(), static_cast<int>(Q.size()));
    return Q;
}

// Set the forces of the points of one set to a function of the index within the set, and the others to zero.
void
set_forces(Pointer<IBTK::LData> F_data,
           IBTK::LDataManager* l_data_manager,
           const std::size_t num_in_set,
           const bool first_set)
{
    boost::multi_array_ref<double, 2>& F = *F_data->getLocalFormVecArray();
    for (const IBTK::LNode* node : l_data_manager->getLMesh(0)->getLocalNodes())
    {
        const std::size_t lag_idx = static_cast<std::size_t>(node->getLagrangianIndex());
        const bool in_first_set = lag_idx < num_in_set;
        F[node->getLocalPETScIndex()][0] =
            in_first_set == first_set ? 1.0 + 0.01 * static_cast<double>(lag_idx % num_in_set) : 0.0;
    }
    F_data->restoreArrays();
}

// The four-point IB kernel as a user-defined kernel with an even stencil size.
double
four_point_kernel(const double r)
{
    const double abs_r = std::abs(r);
    if (abs_r < 1.0) return 0.125 * (3.0 - 2.0 * abs_r + std::sqrt(1.0 + 4.0 * abs_r - 4.0 * abs_r * abs_r));
    if (abs_r < 2.0) return 0.125 * (5.0 - 2.0 * abs_r - std::sqrt(-7.0 + 12.0 * abs_r - 4.0 * abs_r * abs_r));
    return 0.0;
}

double
smooth_field(const double* x)
{
    double value = 1.0;
    for (int d = 0; d < NDIM; ++d) value += 0.3 * std::sin(2.0 * M_PI * (1.0 + d) * x[d] + 0.4 * d);
    return value;
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

// Run the checks for one kernel with a manager of its own, because the manager fixes the interpolation kernel.
void
check_kernel(const std::string& kernel,
             Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
             const int ghost_width,
             const std::array<int, 3>& idx,
             const Points& points)
{
    const int q_idx = idx[0], f_first_idx = idx[1], f_second_idx = idx[2];
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);
    const std::size_t num_in_set = points.stored.size();
    // The wider kernels have larger round-off.
    const double tol = LEInteractor::getStencilSize(kernel) > 4 ? 1.0e-10 : 1.0e-12;

    IBTK::LDataManager* const l_data_manager = IBTK::LDataManager::getManager(
        "LDataManagerPeriodicImages_" + kernel, kernel, kernel, false, IntVector<NDIM>(ghost_width), false);
    l_data_manager->registerLInitStrategy(new PointInitializer(points));
    l_data_manager->setPatchHierarchy(patch_hierarchy);
    l_data_manager->setPatchLevels(0, 0);
    l_data_manager->initializeLevelData(patch_hierarchy, 0, 0.0, false, true);
    l_data_manager->resetHierarchyConfiguration(patch_hierarchy, 0, 0);

    // Check the weights of the field 1 at the stored positions of the first set and the positions of the second set,
    // including a point within round-off of a face shared by patches.
    {
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<CellData<NDIM, double>> q_data = level->getPatch(p())->getPatchData(q_idx);
            q_data->fillAll(1.0);
        }
        Pointer<IBTK::LData> X_initial = l_data_manager->getLData(IBTK::LDataManager::POSN_DATA_NAME, 0);
        Pointer<IBTK::LData> U_initial = l_data_manager->createLData("U_initial", 0, 1);
        l_data_manager->interp(q_idx, U_initial, X_initial, 0, {}, {}, 0.0);
        const std::vector<double> U_stored = gather_by_lagrangian_index(U_initial, l_data_manager, 2 * num_in_set);
        for (std::size_t k = 0; k < U_stored.size(); ++k)
        {
            if (!(std::abs(U_stored[k] - 1.0) <= tol))
            {
                TBOX_ERROR(kernel << ": interpolation of the field 1 at the stored positions gives " << U_stored[k]
                                  << " at point " << k << ", which is at "
                                  << (k < num_in_set ? points.stored[k] : points.image[k - num_in_set]).transpose()
                                  << "\n");
            }
        }
    }

    // Move the first set.
    Pointer<IBTK::LData> X_data = l_data_manager->getLData(IBTK::LDataManager::POSN_DATA_NAME, 0);
    move_points(X_data, l_data_manager, points);
    Pointer<IBTK::LData> U_data = l_data_manager->createLData("U", 0, 1);
    Pointer<IBTK::LData> F_data = l_data_manager->createLData("F", 0, 1);

    Pointer<RefineAlgorithm<NDIM>> ghost_fill_alg = new RefineAlgorithm<NDIM>();
    ghost_fill_alg->registerRefine(q_idx, q_idx, q_idx, nullptr);
    const std::vector<Pointer<RefineSchedule<NDIM>>> ghost_fill_scheds = { ghost_fill_alg->createSchedule(level) };

    // The weights sum to 1.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<CellData<NDIM, double>> q_data = level->getPatch(p())->getPatchData(q_idx);
        q_data->fillAll(1.0);
    }
    l_data_manager->interp(q_idx, U_data, X_data, 0, {}, {}, 0.0);
    const std::vector<double> U_one = gather_by_lagrangian_index(U_data, l_data_manager, 2 * num_in_set);
    for (std::size_t k = 0; k < U_one.size(); ++k)
    {
        if (!(std::abs(U_one[k] - 1.0) <= tol))
        {
            TBOX_ERROR(kernel << ": interpolation of the field 1 gives " << U_one[k] << " at point " << k
                              << ", which is at "
                              << (k < num_in_set ? points.moved[k] : points.image[k - num_in_set]).transpose() << "\n");
        }
    }

    // A periodic field has the same value at the two points of a pair.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, double>> q_data = patch->getPatchData(q_idx);
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const x_lower = pgeom->getXLower();
        const double* const dx = pgeom->getDx();
        q_data->fillAll(0.0);
        for (CellIterator<NDIM> ic(patch->getBox()); ic; ic++)
        {
            double x[NDIM];
            for (int d = 0; d < NDIM; ++d) x[d] = x_lower[d] + dx[d] * (ic()(d) - patch->getBox().lower(d) + 0.5);
            (*q_data)(ic()) = smooth_field(x);
        }
    }
    ghost_fill_scheds[0]->fillData(0.0);
    l_data_manager->interp(q_idx, U_data, X_data, 0, {}, {}, 0.0);
    const std::vector<double> U = gather_by_lagrangian_index(U_data, l_data_manager, 2 * num_in_set);
    for (std::size_t k = 0; k < num_in_set; ++k)
    {
        if (!(std::abs(U[k] - U[k + num_in_set]) <= tol))
        {
            TBOX_ERROR(kernel << ": the interpolated values at the point " << points.moved[k].transpose()
                              << " and at its translate by a period are " << U[k] << " and " << U[k + num_in_set]
                              << ", not equal\n");
        }
    }

    // Spreading the forces of either set gives the same grid function, whose integral is the total force.
    RobinPhysBdryPatchStrategy* const no_bdry_op = nullptr;
    double force_sum = 0.0;
    for (std::size_t k = 0; k < num_in_set; ++k) force_sum += 1.0 + 0.01 * static_cast<double>(k);
    for (int set = 0; set < 2; ++set)
    {
        const int f_idx = set == 0 ? f_first_idx : f_second_idx;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<CellData<NDIM, double>> f_data = level->getPatch(p())->getPatchData(f_idx);
            f_data->fillAll(0.0);
        }
        set_forces(F_data, l_data_manager, num_in_set, set == 0);
        l_data_manager->spread(f_idx, F_data, X_data, no_bdry_op, 0, {}, 0.0);
        const double integral = integrate(level, f_idx);
        if (!(std::abs(integral - force_sum) <= tol * force_sum))
        {
            TBOX_ERROR(kernel << ": spreading the " << (set == 0 ? "first" : "second") << " set gives a total force of "
                              << integral << ", not " << force_sum << "\n");
        }
    }
    double max_difference = 0.0, max_value = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, double>> f_first = patch->getPatchData(f_first_idx);
        Pointer<CellData<NDIM, double>> f_second = patch->getPatchData(f_second_idx);
        for (CellIterator<NDIM> ic(patch->getBox()); ic; ic++)
        {
            max_difference = std::max(max_difference, std::abs((*f_first)(ic()) - (*f_second)(ic())));
            max_value = std::max(max_value, std::abs((*f_second)(ic())));
        }
    }
    max_difference = IBTK_MPI::maxReduction(max_difference);
    max_value = IBTK_MPI::maxReduction(max_value);
    if (!(max_difference <= tol * max_value))
    {
        TBOX_ERROR(kernel << ": the spread forces of the two sets differ by " << max_difference
                          << ", with a largest value of " << max_value << "\n");
    }
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
    LEInteractor::s_kernel_fcn = &four_point_kernel;
    LEInteractor::s_kernel_fcn_stencil_size = 4;
    const std::vector<std::string> kernels = { "USER_DEFINED",
                                               "BSPLINE_4",
                                               "BSPLINE_6",
                                               "PIECEWISE_CUBIC",
                                               "COMPOSITE_BSPLINE_34",
                                               "COMPOSITE_BSPLINE_43",
                                               "COMPOSITE_BSPLINE_56",
                                               "COMPOSITE_BSPLINE_65",
                                               "BSPLINE_3",
                                               "BSPLINE_5",
                                               "COMPOSITE_BSPLINE_23",
                                               "COMPOSITE_BSPLINE_32",
                                               "COMPOSITE_BSPLINE_45",
                                               "COMPOSITE_BSPLINE_54",
                                               "IB_3",
                                               "IB_4",
                                               "IB_5",
                                               "IB_6",
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

    // The field, and the grid functions that store the spread forces of the two sets.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("l_data_manager_periodic_images");
    std::array<int, 3> idx;
    for (int k = 0; k < 3; ++k)
    {
        Pointer<CellVariable<NDIM, double>> var = new CellVariable<NDIM, double>("q_" + std::to_string(k));
        idx[k] = var_db->registerVariableAndContext(var, ctx, IntVector<NDIM>(ghost_width));
    }

    gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);
    for (int k = 0; k < 3; ++k) level->allocatePatchData(idx[k], 0.0);

    Pointer<Database> geometry_db = app_initializer->getComponentDatabase("CartesianGeometry");
    double x_lo[NDIM], x_up[NDIM];
    geometry_db->getDoubleArray("x_lo", x_lo, NDIM);
    geometry_db->getDoubleArray("x_up", x_up, NDIM);
    const Points points = build_points(input_db->getInteger("N"), x_lo, x_up);
    pout << "patches on the level = " << level->getNumberOfPatches() << '\n'
         << "points in each set = " << points.stored.size() << '\n'
         << "points of the first set that cross a periodic boundary = " << points.num_crossing << '\n';
    for (const std::string& kernel : kernels) check_kernel(kernel, patch_hierarchy, ghost_width, idx, points);
}
