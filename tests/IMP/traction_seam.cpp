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

// Check that IMPMethod velocity interpolation and force spreading are independent of the division of the level into
// patches, for material points next to a traction boundary that have crossed a patch boundary without being assigned
// to patches again. The program shifts the points by marker_shift grid cells, interpolates a velocity field, spreads a
// stress, and prints sums and maxima for comparison with the expected output.
//
// With the input key layouts, the program builds one complete stack for each named layout of patches and moves the
// points of every stack to the same designed positions: on and within a few ulps of cell faces and centers, beside and
// beyond the physical boundaries, and at periodic seams. It stops unless U and Grad U are exactly equal at every point
// and the spread force agrees to 1e-12 on every face.

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IMPInitializer.h>
#include <ibamr/IMPMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LMesh.h>
#include <ibtk/LNode.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <LoadBalancer.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
// One complete stack: the fluid solver, the IB method and integrator, and the grid, with object names made distinct by
// a suffix so that several stacks can coexist.
struct Stack
{
    explicit Stack(const libMesh::Parallel::Communicator& comm) : mesh(comm, NDIM)
    {
    }

    std::string layout;
    Mesh mesh;
    Pointer<INSHierarchyIntegrator> ins_integrator;
    Pointer<IMPMethod> ib_method_ops;
    Pointer<IBHierarchyIntegrator> ib_integrator;
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry;
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy;
    Pointer<StandardTagAndInitialize<NDIM>> error_detector;
    Pointer<BergerRigoutsos<NDIM>> box_generator;
    Pointer<LoadBalancer<NDIM>> load_balancer;
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm;
    std::vector<std::unique_ptr<muParserRobinBcCoefs>> u_bc_coef_owners;
    std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs;
};

// The interpolated values at one material point.
struct PointData
{
    std::array<double, NDIM> X{}, U{};
    std::array<double, NDIM * NDIM> Grad_U{};
};

// The results of one stack: the point data by Lagrangian index and the spread force on each face of the domain. The
// force is stored by axis and global face index, with the flag set for the faces that are counted.
struct StackResult
{
    std::map<int, PointData> points;
    std::array<std::vector<double>, NDIM> f;
    std::array<std::vector<char>, NDIM> f_counted;
};

// Build a stack whose object names end with suffix, up to the initialization of the patch hierarchy. The existing
// inputs use the empty suffix, so their names and output do not change.
std::unique_ptr<Stack>
build_stack(const LibMeshInit& init, Pointer<AppInitializer> app_initializer, const std::string& layout)
{
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    const std::string suffix = layout.empty() ? "" : "_" + layout;
    auto stack = std::make_unique<Stack>(init.comm());
    stack->layout = layout;

    // A strip of elements along the lower boundary, across the boundary between two patches.
    std::array<double, 4> strip{};
    input_db->getDoubleArray("strip", strip.data(), 4);
    MeshTools::Generation::build_square(stack->mesh,
                                        static_cast<int>(input_db->getInteger("num_elements_x")),
                                        static_cast<int>(input_db->getInteger("num_elements_y")),
                                        strip[0],
                                        strip[1],
                                        strip[2],
                                        strip[3],
                                        QUAD4);
    stack->mesh.prepare_for_use();

    stack->ins_integrator =
        new INSStaggeredHierarchyIntegrator("INSStaggeredHierarchyIntegrator" + suffix,
                                            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    stack->ib_method_ops = new IMPMethod("IMPMethod" + suffix, app_initializer->getComponentDatabase("IMPMethod"));
    stack->ib_integrator =
        new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator" + suffix,
                                          app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                          stack->ib_method_ops,
                                          stack->ins_integrator);
    stack->grid_geometry = new CartesianGridGeometry<NDIM>("CartesianGeometry" + suffix,
                                                           app_initializer->getComponentDatabase("CartesianGeometry"));
    stack->patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy" + suffix, stack->grid_geometry);
    stack->error_detector =
        new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize" + suffix,
                                           stack->ib_integrator,
                                           app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    stack->box_generator = new BergerRigoutsos<NDIM>();
    stack->load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer" + suffix, app_initializer->getComponentDatabase("LoadBalancer"));
    stack->gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm" + suffix,
                                    app_initializer->getComponentDatabase("GriddingAlgorithm" + suffix),
                                    stack->error_detector,
                                    stack->box_generator,
                                    stack->load_balancer);

    Pointer<IMPInitializer> ib_initializer = new IMPInitializer("IMPInitializer" + suffix,
                                                                app_initializer->getComponentDatabase("IMPInitializer"),
                                                                stack->patch_hierarchy,
                                                                stack->gridding_algorithm);
    ib_initializer->registerMesh(&stack->mesh);
    stack->ib_method_ops->registerLInitStrategy(ib_initializer);

    for (unsigned int d = 0; d < NDIM; ++d)
    {
        stack->u_bc_coef_owners.push_back(std::make_unique<muParserRobinBcCoefs>(
            "u_bc_coefs_" + std::to_string(d) + suffix,
            app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
            stack->grid_geometry));
        stack->u_bc_coefs.push_back(stack->u_bc_coef_owners.back().get());
    }
    stack->ins_integrator->registerPhysicalBoundaryConditions(stack->u_bc_coefs);

    stack->ib_integrator->initializePatchHierarchy(stack->patch_hierarchy, stack->gridding_algorithm);
    stack->ib_method_ops->freeLInitStrategy();
    ib_initializer.setNull();
    return stack;
} // build_stack

// Preprocessing sets up the velocity boundary conditions.
void
preprocess(Stack& stack)
{
    const double current_time = stack.ib_integrator->getIntegratorTime();
    stack.ib_integrator->preprocessIntegrateHierarchy(
        current_time, current_time + stack.ib_integrator->getMaximumTimeStepSize(), 1);
} // preprocess

// Set each coordinate of each point to new_value(Lagrangian index, axis, current value), without assigning the points
// to patches again, and update the ghost positions.
void
set_positions(Stack& stack, const std::function<double(int, unsigned int, double)>& new_value)
{
    const int ln = stack.patch_hierarchy->getFinestLevelNumber();
    LDataManager* l_data_manager = stack.ib_method_ops->getLDataManager();
    Pointer<LData> X_data = l_data_manager->getLData(LDataManager::POSN_DATA_NAME, ln);
    {
        boost::multi_array_ref<double, 2>& X = *X_data->getLocalFormVecArray();
        for (const auto& node : l_data_manager->getLMesh(ln)->getLocalNodes())
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                double& x = X[node->getLocalPETScIndex()][d];
                x = new_value(node->getLagrangianIndex(), d, x);
            }
        }
        X_data->restoreArrays();
    }
    X_data->beginGhostUpdate();
    X_data->endGhostUpdate();
} // set_positions

// The positions along one axis: each face and center and the doubles one, two, and six ulps on either side of it,
// eight golden-ratio positions, and positions beyond each end of the domain.
std::vector<double>
axis_positions(const double x_lo, const double x_up, const int N)
{
    const double h = (x_up - x_lo) / N;
    std::vector<double> positions;
    const auto add_with_neighbors = [&positions](const double x)
    {
        positions.push_back(x);
        double below = x, above = x;
        for (int k = 1; k <= 6; ++k)
        {
            below = std::nextafter(below, -std::numeric_limits<double>::infinity());
            above = std::nextafter(above, std::numeric_limits<double>::infinity());
            if (k == 1 || k == 2 || k == 6)
            {
                positions.push_back(below);
                positions.push_back(above);
            }
        }
    };
    for (int j = 0; j <= N; ++j)
    {
        add_with_neighbors(j == 0 ? x_lo : (j == N ? x_up : x_lo + j * h));
        if (j < N)
        {
            add_with_neighbors(x_lo + (j + 0.5) * h);
        }
    }
    for (int k = 1; k <= 8; ++k)
    {
        const double fraction = k * 0.6180339887498949 - std::floor(k * 0.6180339887498949);
        positions.push_back(x_lo + fraction * (x_up - x_lo));
    }
    for (const double cells : { 0.3, 0.9 })
    {
        positions.push_back(x_lo - cells * h);
        positions.push_back(x_up + cells * h);
    }
    std::vector<double> kept;
    for (const double x : positions)
    {
        if (x >= x_lo - h && x <= x_up + h)
        {
            kept.push_back(x);
        }
    }
    std::sort(kept.begin(), kept.end());
    kept.erase(std::unique(kept.begin(), kept.end()), kept.end());
    return kept;
} // axis_positions

// The designed points, NDIM values each: for each axis, every position along it, combined with positions along the
// other axes in a fixed order.
std::vector<double>
designed_points(const Stack& stack)
{
    std::array<std::vector<double>, NDIM> along;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        along[d] = axis_positions(stack.grid_geometry->getXLower()[d],
                                  stack.grid_geometry->getXUpper()[d],
                                  stack.grid_geometry->getPhysicalDomain()[0].numberCells(d));
    }
    std::vector<double> X;
    for (unsigned int lead = 0; lead < NDIM; ++lead)
    {
        for (std::size_t k = 0; k < along[lead].size(); ++k)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::size_t n = along[d].size();
                X.push_back(d == lead ? along[d][k] : along[d][(7 * k + 13 * d + 5) % n]);
            }
        }
    }
    return X;
} // designed_points

// Move the points toward the designed positions, point n to position n modulo the number of positions, in steps of at
// most one cell along each axis with a regrid after each step. A regrid clamps a point beyond a non-periodic boundary,
// or too close to the upper one, with clampToDomain(), and wraps one beyond a periodic boundary, so the steps end at
// the clamped position or, along a periodic axis, at the nearest position inside the domain. Return the number of
// regrids.
int
move_points_in_steps(Stack& stack, const std::vector<double>& design)
{
    const int ln = stack.patch_hierarchy->getFinestLevelNumber();
    LDataManager* l_data_manager = stack.ib_method_ops->getLDataManager();
    const double* const x_lo = stack.grid_geometry->getXLower();
    const double* const x_up = stack.grid_geometry->getXUpper();
    const double* const dx = stack.grid_geometry->getDx();
    const IntVector<NDIM> periodic_shift = stack.grid_geometry->getPeriodicShift();
    const int num_design = static_cast<int>(design.size() / NDIM);

    // The start positions and the positions that the steps end at.
    std::map<int, std::array<double, NDIM>> start, end;
    Pointer<LData> X_data = l_data_manager->getLData(LDataManager::POSN_DATA_NAME, ln);
    {
        boost::multi_array_ref<double, 2>& X = *X_data->getLocalFormVecArray();
        for (const auto& node : l_data_manager->getLMesh(ln)->getLocalNodes())
        {
            const int lag_idx = node->getLagrangianIndex();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                start[lag_idx][d] = X[node->getLocalPETScIndex()][d];
            }
        }
        X_data->restoreArrays();
    }
    if (static_cast<int>(start.size()) < num_design)
    {
        TBOX_ERROR("traction_seam: the mesh has " << start.size() << " points, fewer than the " << num_design
                                                  << " designed positions\n");
    }
    double max_cells = 0.0;
    for (const auto& [lag_idx, X0] : start)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const double target = design[NDIM * (lag_idx % num_design) + d];
            if (periodic_shift(d))
            {
                end[lag_idx][d] =
                    target < x_lo[d] ? x_lo[d] : (target < x_up[d] ? target : std::nextafter(x_up[d], x_lo[d]));
            }
            else
            {
                end[lag_idx][d] = IndexUtilities::clampToDomain(target, x_lo[d], x_up[d]);
            }
            max_cells = std::max(max_cells, std::abs(end[lag_idx][d] - X0[d]) / dx[d]);
        }
    }
    const int num_steps = std::max(1, static_cast<int>(std::ceil(max_cells)));
    for (int step = 1; step <= num_steps; ++step)
    {
        set_positions(stack,
                      [&](const int lag_idx, const unsigned int d, const double)
                      {
                          // The last step assigns the end position, so that it is exact.
                          if (step == num_steps)
                          {
                              return end.at(lag_idx)[d];
                          }
                          return start.at(lag_idx)[d] +
                                 (end.at(lag_idx)[d] - start.at(lag_idx)[d]) * (static_cast<double>(step) / num_steps);
                      });
        stack.ib_integrator->regridHierarchy();
    }
    return num_steps;
} // move_points_in_steps

// Move the points to their designed positions, at most one cell from where the steps ended, without a regrid. Points
// beyond a periodic boundary are then reached through their periodic images.
void
move_points_to_design(Stack& stack, const std::vector<double>& design)
{
    const int num_design = static_cast<int>(design.size() / NDIM);
    set_positions(stack,
                  [&](const int lag_idx, const unsigned int d, const double)
                  { return design[NDIM * (lag_idx % num_design) + d]; });
} // move_points_to_design

// Interpolate a velocity field and spread a stress. Print the sums and maxima of the results if print_values is true,
// and store the point data and the spread force in result if it is given. The velocity is a
// function of the position in the frame of the patch, or in the frame of the level if level_frame_field is true, where
// it is periodic along the periodic axes.
void
check_partition_independence(Stack& stack,
                             const bool print_values,
                             const bool level_frame_field,
                             StackResult* const result)
{
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = stack.patch_hierarchy;
    Pointer<INSHierarchyIntegrator> ins_integrator = stack.ins_integrator;
    Pointer<IBHierarchyIntegrator> ib_integrator = stack.ib_integrator;
    Pointer<IMPMethod> ib_method_ops = stack.ib_method_ops;
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = stack.grid_geometry;
    const int ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
    const double time = ib_integrator->getIntegratorTime();
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<Variable<NDIM>> u_var = ins_integrator->getVelocityVariable();
    const int u_idx = var_db->registerVariableAndContext(
        u_var, var_db->getContext("traction_seam"), ib_method_ops->getMinimumGhostCellWidth());
    const int f_idx = var_db->registerClonedPatchDataIndex(u_var, u_idx);
    level->allocatePatchData(u_idx, time);
    level->allocatePatchData(f_idx, time);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<SideData<NDIM, double>> f_data = level->getPatch(p())->getPatchData(f_idx);
        f_data->fillAll(0.0);
    }

    // The velocity is a function of position. The operator with the boundary conditions of the fluid solver sets its
    // ghost values outside the domain.
    RobinPhysBdryPatchStrategy* bdry_op = ib_integrator->getVelocityPhysBdryOp();
    bdry_op->setPatchDataIndex(u_idx);
    bdry_op->setHomogeneousBc(false);
    Pointer<RefineAlgorithm<NDIM>> ghost_fill_alg = new RefineAlgorithm<NDIM>();
    ghost_fill_alg->registerRefine(u_idx, u_idx, u_idx, nullptr);
    std::vector<Pointer<RefineSchedule<NDIM>>> ghost_fill_scheds(ln + 1);
    ghost_fill_scheds[ln] = ghost_fill_alg->createSchedule(level, bdry_op);
    const IntVector<NDIM> periodic_shift = grid_geom->getPeriodicShift();
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const x_lower = level_frame_field ? grid_geom->getXLower() : pgeom->getXLower();
        const double* const dx = level_frame_field ? grid_geom->getDx() : pgeom->getDx();
        const hier::Index<NDIM>& index_lower =
            level_frame_field ? grid_geom->getPhysicalDomain()[0].lower() : patch_box.lower();
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        u_data->fillAll(0.0);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                double phase = 1.7 * axis;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    if (level_frame_field && periodic_shift(d))
                    {
                        // The index modulo the number of cells gives the two copies of a periodic face equal values.
                        const int n = grid_geom->getPhysicalDomain()[0].numberCells(d);
                        const int k = ((b()(d) - index_lower(d)) % n + n) % n;
                        const double x = x_lower[d] + dx[d] * (k + (d == axis ? 0.0 : 0.5));
                        phase += 2.0 * M_PI * (x - x_lower[d]) / (grid_geom->getXUpper()[d] - x_lower[d]);
                    }
                    else
                    {
                        const double x = x_lower[d] + dx[d] * (b()(d) - index_lower(d) + (d == axis ? 0.0 : 0.5));
                        phase += (3.1 + 2.2 * d) * x;
                    }
                }
                u_data->getArrayData(axis)(b(), 0) = std::sin(phase);
            }
        }
    }
    ghost_fill_scheds[ln]->fillData(time);

    // Interpolate u, and set a stress that depends only on the Lagrangian index.
    ib_method_ops->interpolateVelocity(u_idx, {}, ghost_fill_scheds, time);
    LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    Pointer<LData> U_data = l_data_manager->getLData(LDataManager::VEL_DATA_NAME, ln);
    Pointer<LData> Grad_U_data = l_data_manager->getLData("Grad_U", ln);
    Pointer<LData> tau_data = l_data_manager->getLData("tau", ln);
    Pointer<LData> X_data = l_data_manager->getLData(LDataManager::POSN_DATA_NAME, ln);
    double num_points = 0.0, num_not_finite = 0.0;
    double U_sum_squares = 0.0, U_max = 0.0, Grad_U_sum_squares = 0.0, Grad_U_max = 0.0;
    {
        boost::multi_array_ref<double, 2>& U = *U_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& Grad_U = *Grad_U_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& tau = *tau_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& X = *X_data->getLocalFormVecArray();
        for (const auto& node : l_data_manager->getLMesh(ln)->getLocalNodes())
        {
            const int lag_idx = node->getLagrangianIndex();
            const int petsc_idx = node->getLocalPETScIndex();
            num_points += 1.0;
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                if (result)
                {
                    PointData& point = result->points[lag_idx];
                    point.X[i] = X[petsc_idx][i];
                    point.U[i] = U[petsc_idx][i];
                    for (unsigned int j = 0; j < NDIM; ++j)
                    {
                        point.Grad_U[NDIM * i + j] = Grad_U[petsc_idx][NDIM * i + j];
                    }
                }
                if (!std::isfinite(U[petsc_idx][i]))
                {
                    num_not_finite += 1.0;
                }
                U_sum_squares += U[petsc_idx][i] * U[petsc_idx][i];
                U_max = std::max(U_max, std::abs(U[petsc_idx][i]));
                for (unsigned int j = 0; j < NDIM; ++j)
                {
                    const double g = Grad_U[petsc_idx][NDIM * i + j];
                    if (!std::isfinite(g))
                    {
                        num_not_finite += 1.0;
                    }
                    Grad_U_sum_squares += g * g;
                    Grad_U_max = std::max(Grad_U_max, std::abs(g));
                    tau[petsc_idx][NDIM * i + j] = std::cos(1.3 * lag_idx + 0.7 * (NDIM * i + j) + 0.2);
                }
            }
        }
        U_data->restoreArrays();
        Grad_U_data->restoreArrays();
        tau_data->restoreArrays();
        X_data->restoreArrays();
    }

    if (IBTK_MPI::sumReduction(num_not_finite) != 0.0)
    {
        TBOX_ERROR("traction_seam: the interpolated values of layout " << stack.layout << " are not finite\n");
    }
    if (print_values)
    {
        pout << std::setprecision(14);
        pout << "points = " << IBTK_MPI::sumReduction(num_points) << '\n';
        pout << "U_sum_squares = " << IBTK_MPI::sumReduction(U_sum_squares) << '\n';
        pout << "U_max = " << IBTK_MPI::maxReduction(U_max) << '\n';
        pout << "Grad_U_sum_squares = " << IBTK_MPI::sumReduction(Grad_U_sum_squares) << '\n';
        pout << "Grad_U_max = " << IBTK_MPI::maxReduction(Grad_U_max) << '\n';
    }

    // Spread the stress. The spread force is zero outside the domain. Each face is counted once.
    ib_method_ops->spreadForce(f_idx, bdry_op, {}, time);
    double f_sum_squares = 0.0, f_max = 0.0;
    if (result)
    {
        const Box<NDIM>& domain_box = grid_geom->getPhysicalDomain()[0];
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            std::size_t num_faces = 1;
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                num_faces *= static_cast<std::size_t>(domain_box.numberCells(d) + (d == axis ? 1 : 0));
            }
            result->f[axis].assign(num_faces, 0.0);
            result->f_counted[axis].assign(num_faces, 0);
        }
    }
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); b; b++)
            {
                // Faces on a patch boundary that is not a physical boundary are shared with another patch.
                const bool on_lower_face = b()(axis) == patch->getBox().lower(axis);
                const bool shared = on_lower_face && !Pointer<CartesianPatchGeometry<NDIM>>(patch->getPatchGeometry())
                                                          ->getTouchesRegularBoundary(axis, 0);
                const double v = f_data->getArrayData(axis)(b(), 0);
                if (shared)
                {
                    continue;
                }
                f_sum_squares += v * v;
                f_max = std::max(f_max, std::abs(v));
                if (result)
                {
                    const Box<NDIM>& domain_box = grid_geom->getPhysicalDomain()[0];
                    std::size_t face = 0, stride = 1;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        face += stride * static_cast<std::size_t>(b()(d) - domain_box.lower(d));
                        stride *= static_cast<std::size_t>(domain_box.numberCells(d) + (d == axis ? 1 : 0));
                    }
                    result->f[axis][face] = v;
                    result->f_counted[axis][face] = 1;
                }
            }
        }
    }
    level->deallocatePatchData(u_idx);
    level->deallocatePatchData(f_idx);

    if (print_values)
    {
        pout << "f_sum_squares = " << IBTK_MPI::sumReduction(f_sum_squares) << '\n';
        pout << "f_max = " << IBTK_MPI::maxReduction(f_max) << '\n';
        plog << std::flush;
    }
    return;
} // check_partition_independence

// Stop unless U and Grad U are exactly equal at every point of the two results and the spread force agrees to 1e-12
// relative to its largest value on every face. Print the number of points compared and one sum per quantity.
void
compare_layouts(const std::string& name_a, const StackResult& a, const std::string& name_b, const StackResult& b)
{
    const auto fail = [&](const std::string& what, const int lag_idx, const PointData& pa, const PointData& pb)
    {
        std::ostringstream message;
        message << std::setprecision(17) << "traction_seam: " << what << " differs between layouts " << name_a
                << " and " << name_b << " at Lagrangian index " << lag_idx << "\n  X =";
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            message << ' ' << pa.X[d];
        }
        message << "\n  X of " << name_b << " =";
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            message << ' ' << pb.X[d];
        }
        message << "\n  U of " << name_a << " =";
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            message << ' ' << pa.U[d];
        }
        message << "\n  U of " << name_b << " =";
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            message << ' ' << pb.U[d];
        }
        message << "\n  Grad_U of " << name_a << " =";
        for (unsigned int k = 0; k < NDIM * NDIM; ++k)
        {
            message << ' ' << pa.Grad_U[k];
        }
        message << "\n  Grad_U of " << name_b << " =";
        for (unsigned int k = 0; k < NDIM * NDIM; ++k)
        {
            message << ' ' << pb.Grad_U[k];
        }
        message << '\n';
        TBOX_ERROR(message.str());
    };
    if (a.points.size() != b.points.size())
    {
        TBOX_ERROR("traction_seam: layout " << name_a << " has " << a.points.size() << " points, and layout " << name_b
                                            << " has " << b.points.size() << "\n");
    }
    double U_sum_squares = 0.0, Grad_U_sum_squares = 0.0;
    for (const auto& [lag_idx, pa] : a.points)
    {
        const auto it = b.points.find(lag_idx);
        if (it == b.points.end())
        {
            TBOX_ERROR("traction_seam: layout " << name_b << " has no point with Lagrangian index " << lag_idx << "\n");
        }
        const PointData& pb = it->second;
        if (pa.X != pb.X)
        {
            fail("the position", lag_idx, pa, pb);
        }
        if (pa.U != pb.U)
        {
            fail("U", lag_idx, pa, pb);
        }
        if (pa.Grad_U != pb.Grad_U)
        {
            fail("Grad_U", lag_idx, pa, pb);
        }
        for (const double u : pa.U)
        {
            U_sum_squares += u * u;
        }
        for (const double g : pa.Grad_U)
        {
            Grad_U_sum_squares += g * g;
        }
    }

    double f_sum_squares = 0.0, f_max = 0.0, f_max_difference = 0.0;
    std::size_t num_faces = 0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (std::size_t face = 0; face < a.f[axis].size(); ++face)
        {
            if (a.f_counted[axis][face] != b.f_counted[axis][face])
            {
                TBOX_ERROR("traction_seam: the layouts count different faces of axis " << axis << "\n");
            }
            if (!a.f_counted[axis][face])
            {
                continue;
            }
            if (!std::isfinite(a.f[axis][face]) || !std::isfinite(b.f[axis][face]))
            {
                TBOX_ERROR("traction_seam: the spread force is not finite on face " << face << " of axis " << axis
                                                                                    << "\n");
            }
            ++num_faces;
            f_sum_squares += a.f[axis][face] * a.f[axis][face];
            f_max = std::max({ f_max, std::abs(a.f[axis][face]), std::abs(b.f[axis][face]) });
            f_max_difference = std::max(f_max_difference, std::abs(a.f[axis][face] - b.f[axis][face]));
        }
    }
    if (f_max_difference > 1.0e-12 * f_max)
    {
        std::ostringstream message;
        message << std::setprecision(17) << "traction_seam: the spread force differs between layouts " << name_a
                << " and " << name_b << " by " << f_max_difference << ", more than 1e-12 times its largest value "
                << f_max << "\n";
        TBOX_ERROR(message.str());
    }

    pout << std::setprecision(10);
    pout << "points compared = " << a.points.size() << '\n';
    pout << "faces compared = " << num_faces << '\n';
    pout << "U_sum_squares = " << U_sum_squares << '\n';
    pout << "Grad_U_sum_squares = " << Grad_U_sum_squares << '\n';
    pout << "f_sum_squares = " << f_sum_squares << '\n';
} // compare_layouts
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const LibMeshInit& init = ibtk_init.getLibMeshInit();

    // Suppress warnings caused by running without Silo.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        if (!input_db->keyExists("layouts"))
        {
            std::unique_ptr<Stack> stack = build_stack(init, app_initializer, "");
            // Regrid, as advanceHierarchy() does at the initial time, to set up the Lagrangian indices in ghost cells.
            stack->ib_integrator->regridHierarchy();
            preprocess(*stack);

            // Move the points by marker_shift grid cells along each axis, without assigning them to patches again.
            std::array<double, NDIM> marker_shift{};
            input_db->getDoubleArray("marker_shift", marker_shift.data(), NDIM);
            const double* const dx = stack->grid_geometry->getDx();
            set_positions(*stack,
                          [&](const int, const unsigned int d, const double x) { return x + marker_shift[d] * dx[d]; });

            check_partition_independence(*stack, /*print_values*/ true, false, nullptr);
            return 0;
        }

        // One stack for each layout, with the points moved to the same designed positions.
        std::vector<std::string> layouts(input_db->getArraySize("layouts"));
        input_db->getStringArray("layouts", layouts.data(), static_cast<int>(layouts.size()));
        if (layouts.size() != 2)
        {
            TBOX_ERROR("traction_seam: the input key layouts must name two layouts\n");
        }
        std::vector<std::unique_ptr<Stack>> stacks;
        std::vector<StackResult> results(layouts.size());
        for (std::size_t n = 0; n < layouts.size(); ++n)
        {
            stacks.push_back(build_stack(init, app_initializer, layouts[n]));
            Stack& stack = *stacks.back();
            const std::vector<double> design = designed_points(stack);
            const int num_regrids = move_points_in_steps(stack, design);
            preprocess(stack);
            move_points_to_design(stack, design);
            pout << "layout " << layouts[n] << ": regrids = " << num_regrids << '\n';
            check_partition_independence(stack, /*print_values*/ false, true, &results[n]);
        }
        compare_layouts(layouts[0], results[0], layouts[1], results[1]);
        // Destroy the stacks in the reverse order of construction.
        while (!stacks.empty())
        {
            stacks.pop_back();
        }
    }
    return 0;
} // main
