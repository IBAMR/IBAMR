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

// Interpolate a discretely divergence-free velocity with a divergence-preserving kernel at points near the physical
// boundaries, using the velocity boundary operator of the IB integrator, and check that the divergence of the
// interpolated velocity is small relative to max |u| / h, where h is the cell width.

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/LMesh.h>
#include <ibtk/LNode.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/RobinPhysBdryPatchStrategy.h>
#include <ibtk/ibtk_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ArrayData.h>
#include <BergerRigoutsos.h>
#include <BoundaryBox.h>
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
#include <iomanip>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
struct MarkerParameters
{
    int finest_ln;
    int num_cells;
};

// The divergence of the interpolated velocity is a centered difference over EPS cell widths, evaluated at these
// depths (in cell widths) from each boundary. No evaluation point is within EPS of a kink of the kernel.
constexpr double EPS = 1.0e-3;
constexpr std::array<double, 3> DEPTHS = { 0.2, 0.7, 1.3 };

std::vector<IBTK::Point>
generate_evaluation_points(const int num_cells)
{
    const double h = 1.0 / num_cells;
    const std::array<double, 3> along = { 0.3 * h, (0.5 * num_cells + 0.2) * h, 1.0 - 0.3 * h };
    std::vector<IBTK::Point> evaluation_points;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (const bool upper : { false, true })
        {
            for (const double depth : DEPTHS)
            {
                // Each axis along the boundary takes one of the positions along it.
                for (int combination = 0; combination < (NDIM == 2 ? 3 : 9); ++combination)
                {
                    IBTK::Point X;
                    int digits = combination;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            X[d] = upper ? 1.0 - depth * h : depth * h;
                        }
                        else
                        {
                            X[d] = along[digits % 3];
                            digits /= 3;
                        }
                    }
                    evaluation_points.push_back(X);
                }
            }
        }
    }
    return evaluation_points;
} // generate_evaluation_points

// Place markers around each evaluation point, so that the divergence of the interpolated velocity can be computed
// there.
void
generate_markers(const unsigned int& /*strct_num*/,
                 const int& ln,
                 int& num_vertices,
                 std::vector<IBTK::Point>& vertex_posn,
                 void* ctx)
{
    const MarkerParameters& params = *static_cast<const MarkerParameters*>(ctx);
    vertex_posn.clear();
    if (ln == params.finest_ln)
    {
        // For each evaluation point, the points displaced by +EPS and -EPS cell widths along each axis.
        const double h = 1.0 / params.num_cells;
        for (const IBTK::Point& evaluation_point : generate_evaluation_points(params.num_cells))
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                for (const double sign : { 1.0, -1.0 })
                {
                    IBTK::Point X = evaluation_point;
                    X[d] += sign * EPS * h;
                    vertex_posn.push_back(X);
                }
            }
        }
    }
    num_vertices = static_cast<int>(vertex_posn.size());
    return;
} // generate_markers

// Register and allocate velocity data on the level, with the ghost width that the IB integrator registers: the ghost
// width that the IB method needs, plus the stencil width of the velocity boundary operator if the divergence-free
// velocity extension is used. Return its patch data index.
int
allocate_velocity_data(Pointer<PatchLevel<NDIM>> level,
                       Pointer<INSHierarchyIntegrator> ins_integrator,
                       Pointer<IBHierarchyIntegrator> ib_integrator,
                       Pointer<IBMethod> ib_method_ops,
                       const bool divergence_free_extension,
                       const double time)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<Variable<NDIM>> u_var = ins_integrator->getVelocityVariable();
    IntVector<NDIM> ghost_width = ib_method_ops->getMinimumGhostCellWidth();
    if (divergence_free_extension)
    {
        ghost_width += ib_integrator->getVelocityPhysBdryOp()->getRefineOpStencilWidth();
    }
    const int u_idx =
        var_db->registerVariableAndContext(u_var, var_db->getContext("interpolated_velocity_divergence"), ghost_width);
    level->allocatePatchData(u_idx, time);
    return u_idx;
} // allocate_velocity_data

// Return, for each axis, whether the velocity condition on the lower boundary normal to it leaves the normal velocity
// unprescribed (b = 1), from the velocity boundary coefficients of the fluid solver. The conditions are constant
// along each boundary.
std::array<bool, NDIM>
get_open_lower_boundaries(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                          Pointer<INSHierarchyIntegrator> ins_integrator,
                          const double time)
{
    const std::vector<RobinBcCoefStrategy<NDIM>*>& bc_coefs = ins_integrator->getPhysicalBoundaryConditions();
    std::array<int, NDIM> is_open;
    is_open.fill(0);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(patch_hierarchy->getFinestLevelNumber());
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const tbox::Array<BoundaryBox<NDIM>> boundary_boxes =
            PhysicalBoundaryUtilities::getPhysicalBoundaryCodim1Boxes(*patch);
        for (int k = 0; k < boundary_boxes.size(); ++k)
        {
            const unsigned int location_index = boundary_boxes[k].getLocationIndex();
            if (location_index % 2 != 0)
            {
                continue;
            }
            const unsigned int axis = location_index / 2;
            const Box<NDIM> coef_box = PhysicalBoundaryUtilities::makeSideBoundaryCodim1Box(boundary_boxes[k]);
            Pointer<ArrayData<NDIM, double>> a_data = new ArrayData<NDIM, double>(coef_box, 1);
            Pointer<ArrayData<NDIM, double>> b_data = new ArrayData<NDIM, double>(coef_box, 1);
            Pointer<ArrayData<NDIM, double>> g_data = new ArrayData<NDIM, double>(coef_box, 1);
            bc_coefs[axis]->setBcCoefs(
                a_data, b_data, g_data, Pointer<Variable<NDIM>>(), *patch, boundary_boxes[k], time);
            if ((*a_data)(coef_box.lower(), 0) == 0.0 && (*b_data)(coef_box.lower(), 0) == 1.0)
            {
                is_open[axis] = 1;
            }
        }
    }
    IBTK_MPI::maxReduction(is_open.data(), NDIM);
    std::array<bool, NDIM> open_lower;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        open_lower[d] = is_open[d] != 0;
    }
    return open_lower;
} // get_open_lower_boundaries

// A smooth potential (stream function in 2D, vector potential in 3D) whose discrete curl is the test velocity. It
// vanishes on the boundaries where the normal velocity is prescribed, so that the velocity has no normal component
// there.
double
potential(const unsigned int component, const std::array<double, NDIM>& x, const std::array<bool, NDIM>& open_lower)
{
    double value = 1.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (NDIM == 2 || d != component)
        {
            value *= open_lower[d] ? 1.0 - x[d] : x[d] * (1.0 - x[d]);
        }
    }
#if (NDIM == 2)
    return value * (1.0 + 0.6 * x[0] + 0.4 * std::sin(3.0 * x[1] + 0.5));
#else
    return value * (0.5 + std::cos((1.1 + 0.4 * component) * x[0] + (0.7 + 0.5 * component) * x[1] +
                                   (0.3 + 0.8 * component) * x[2] + 0.2 * component));
#endif
} // potential

// The difference of the potential over a cell width in direction d, centered at x.
double
potential_difference(const unsigned int component,
                     const std::array<double, NDIM>& x,
                     const unsigned int d,
                     const double dx,
                     const std::array<bool, NDIM>& open_lower)
{
    std::array<double, NDIM> x_upper = x, x_lower = x;
    x_upper[d] += 0.5 * dx;
    x_lower[d] -= 0.5 * dx;
    return (potential(component, x_upper, open_lower) - potential(component, x_lower, open_lower)) / dx;
} // potential_difference

// Interpolate a velocity that is discretely divergence-free in the domain, as
// the differences of a node-centered stream function in 2D and as the discrete
// curl of an edge-centered vector potential in 3D, near the physical
// boundaries with the IB integrator's velocity boundary operator and the IB
// kernel, and compare the centered-difference divergence of the interpolated
// velocity with zero.
void
check_interpolated_divergence(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                              Pointer<INSHierarchyIntegrator> ins_integrator,
                              Pointer<IBHierarchyIntegrator> ib_integrator,
                              Pointer<IBMethod> ib_method_ops,
                              const bool divergence_free_extension,
                              Pointer<Database> input_db)
{
    const int ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
    const double time = ib_integrator->getIntegratorTime();
    const int u_idx =
        allocate_velocity_data(level, ins_integrator, ib_integrator, ib_method_ops, divergence_free_extension, time);
    const std::array<bool, NDIM> open_lower = get_open_lower_boundaries(patch_hierarchy, ins_integrator, time);
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = patch_hierarchy->getGridGeometry();
    const double* const domain_x_lower = grid_geom->getXLower();

    RobinPhysBdryPatchStrategy* bdry_op = ib_integrator->getVelocityPhysBdryOp();
    bdry_op->setPatchDataIndex(u_idx);
    bdry_op->setHomogeneousBc(false);
    Pointer<RefineAlgorithm<NDIM>> ghost_fill_alg = new RefineAlgorithm<NDIM>();
    ghost_fill_alg->registerRefine(u_idx, u_idx, u_idx, nullptr);
    std::vector<Pointer<RefineSchedule<NDIM>>> ghost_fill_scheds(ln + 1);
    ghost_fill_scheds[ln] = ghost_fill_alg->createSchedule(level, bdry_op);

    // Set the velocity on the faces of each patch, and find its largest
    // divergence in the cells of the domain.
    double max_grid_div = 0.0;
    double max_abs_u = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        u_data->fillAll(0.0);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                std::array<double, NDIM> x;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    x[d] = domain_x_lower[d] + dx[d] * (b()(d) + (d == axis ? 0.0 : 0.5));
                }
#if (NDIM == 2)
                const double sign = axis == 0 ? 1.0 : -1.0;
                u_data->getArrayData(axis)(b(), 0) =
                    sign * potential_difference(0, x, 1 - axis, dx[1 - axis], open_lower);
#else
                const unsigned int next = (axis + 1) % 3, previous = (axis + 2) % 3;
                u_data->getArrayData(axis)(b(), 0) = potential_difference(previous, x, next, dx[next], open_lower) -
                                                     potential_difference(next, x, previous, dx[previous], open_lower);
#endif
                max_abs_u = std::max(max_abs_u, std::abs(u_data->getArrayData(axis)(b(), 0)));
            }
        }
        for (Box<NDIM>::Iterator b(patch_box); b; b++)
        {
            double div = 0.0;
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                hier::Index<NDIM> upper = b();
                upper(axis) += 1;
                div += (u_data->getArrayData(axis)(upper, 0) - u_data->getArrayData(axis)(b(), 0)) / dx[axis];
            }
            if (!std::isfinite(div))
            {
                TBOX_ERROR("the divergence of the velocity in the domain is not finite in the cell " << b() << '\n');
            }
            max_grid_div = std::max(max_grid_div, std::abs(div));
        }
    }
    max_grid_div = IBTK_MPI::maxReduction(max_grid_div);
    max_abs_u = IBTK_MPI::maxReduction(max_abs_u);
    constexpr double GRID_DIV_BOUND = 1.0e-9;
    if (!(max_grid_div < GRID_DIV_BOUND))
    {
        TBOX_ERROR("max |div u| of the velocity in the domain is " << max_grid_div << '\n');
    }

    // Interpolate u at the points around each evaluation point.
    LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    Pointer<LData> X_data = l_data_manager->getLData("X", ln);
    Pointer<LData> U_data = l_data_manager->createLData("interpolated_velocity_divergence_U", ln, NDIM);
    l_data_manager->interp(u_idx, U_data, X_data, ln, {}, ghost_fill_scheds, time);
    const int num_evaluation_points = static_cast<int>(generate_evaluation_points(input_db->getInteger("N")).size());
    const int num_points = 2 * NDIM * num_evaluation_points;
    std::vector<double> U(num_points * NDIM, 0.0);
    {
        boost::multi_array_ref<double, 2>& U_local = *U_data->getLocalFormVecArray();
        for (const auto& node : l_data_manager->getLMesh(ln)->getLocalNodes())
        {
            const int lag_idx = node->getLagrangianIndex();
            TBOX_ASSERT(lag_idx >= 0 && lag_idx < num_points);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                U[lag_idx * NDIM + d] = U_local[node->getLocalPETScIndex()][d];
            }
        }
        U_data->restoreArrays();
    }
    IBTK_MPI::sumReduction(U.data(), static_cast<int>(U.size()));

    // The centered-difference divergence of the interpolated velocity at each evaluation point.
    const double eps_length = EPS / input_db->getInteger("N");
    double max_interpolated_div = 0.0;
    for (int point = 0; point < num_evaluation_points; ++point)
    {
        double div = 0.0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const int plus = 2 * NDIM * point + 2 * d;
            div += (U[plus * NDIM + d] - U[(plus + 1) * NDIM + d]) / (2.0 * eps_length);
        }
        if (!std::isfinite(div))
        {
            TBOX_ERROR("the divergence of the interpolated velocity is not finite at evaluation point " << point
                                                                                                        << '\n');
        }
        max_interpolated_div = std::max(max_interpolated_div, std::abs(div));
    }
    level->deallocatePatchData(u_idx);

    // The divergence is scaled by h / max |u| so that the bound does not depend on the magnitude of the velocity or on
    // the cell width h. The centered difference over EPS cell widths amplifies roundoff in the interpolated velocity
    // by 1/EPS, so the bound is well above roundoff and far below the divergence that the option off produces.
    const double h = 1.0 / input_db->getInteger("N");
    const double relative_div = max_interpolated_div * h / max_abs_u;
    constexpr double RELATIVE_DIV_BOUND = 1.0e-9;
    if (!(relative_div < RELATIVE_DIV_BOUND))
    {
        TBOX_ERROR("max |div| of the interpolated velocity / (max |u| / h) is " << relative_div << ", which exceeds "
                                                                                << RELATIVE_DIV_BOUND << '\n');
    }
    pout << "patches on the level = " << level->getNumberOfPatches() << '\n';
    pout << "evaluation points = " << num_evaluation_points << '\n';
    pout << std::setprecision(8) << "max |u| in the domain = " << max_abs_u << '\n';
    return;
} // check_interpolated_divergence
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

#ifndef IBTK_HAVE_SILO
    // Suppress warnings caused by running without Silo.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        Pointer<INSHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        Pointer<IBHierarchyIntegrator> ib_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              ins_integrator);
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               ib_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        const bool divergence_free_extension = app_initializer->getComponentDatabase("IBHierarchyIntegrator")
                                                   ->getBoolWithDefault("divergence_free_velocity_extension", false);
        MarkerParameters marker_params{ input_db->getInteger("MAX_LEVELS") - 1, input_db->getInteger("N") };
        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        ib_initializer->setStructureNamesOnLevel(marker_params.finest_ln, { "markers" });
        ib_initializer->registerInitStructureFunction(generate_markers, &marker_params);
        ib_method_ops->registerLInitStrategy(ib_initializer);
        ib_method_ops->registerIBLagrangianForceFunction(new IBStandardForceGen());

        Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
            "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
        ins_integrator->registerVelocityInitialConditions(u_init);
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            u_bc_coefs[d] =
                new muParserRobinBcCoefs("u_bc_coefs_" + std::to_string(d),
                                         app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                                         grid_geometry);
        }
        ins_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

        ib_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();

        // As advanceHierarchy() does at the initial time, regrid first, which
        // also sets up the Lagrangian indices in the ghost cells that spreading
        // uses; the velocity boundary conditions are set up at the start of a
        // time step.
        ib_integrator->regridHierarchy();
        const double current_time = ib_integrator->getIntegratorTime();
        ib_integrator->preprocessIntegrateHierarchy(
            current_time, current_time + ib_integrator->getMaximumTimeStepSize(), 1);
        check_interpolated_divergence(
            patch_hierarchy, ins_integrator, ib_integrator, ib_method_ops, divergence_free_extension, input_db);

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
