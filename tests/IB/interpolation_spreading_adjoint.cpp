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

// Check that IB force spreading is the adjoint of IB velocity interpolation,
// including the velocity boundary conditions at physical boundaries, on a level
// with several patches: (F, J u) = (S F, u).
//
// With check_interpolated_divergence = TRUE in the input, instead interpolate a
// discretely divergence-free velocity with a divergence-preserving kernel at
// points near the physical boundaries, and report the divergence of the
// interpolated velocity.

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
#include <ibtk/RobinPhysBdryPatchStrategy.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

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
    bool probes;
};

// The divergence of the interpolated velocity is the centered difference over
// EPS cell widths. The velocity is interpolated at points at these depths (in
// cell widths) from a boundary and, along the boundary, near its lower end, near
// its middle, and near its upper end. The interpolation weights are not smooth
// at odd multiples of half a cell width, and none of the points is within EPS
// cell widths of one.
constexpr double EPS = 1.0e-3;
constexpr std::array<double, 3> PROBE_DEPTHS = { 0.2, 0.7, 1.3 };

std::vector<IBTK::Point>
generate_probe_points(const int num_cells)
{
    const double h = 1.0 / num_cells;
    const std::array<double, 3> along = { 0.3 * h, (0.5 * num_cells + 0.2) * h, 1.0 - 0.3 * h };
    std::vector<IBTK::Point> probes;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        for (const bool upper : { false, true })
        {
            for (const double depth : PROBE_DEPTHS)
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
                    probes.push_back(X);
                }
            }
        }
    }
    return probes;
} // generate_probe_points

// Place markers on lines parallel to each coordinate axis, a fraction of a grid
// cell away from each pair of boundaries that meet along an edge of the unit
// domain, so that the IB kernels overlap the physical boundaries, their
// corners (and edges), and the boundaries between patches. When probing the
// divergence of the interpolated velocity, place markers around the probe points.
void
generate_markers(const unsigned int& /*strct_num*/,
                 const int& ln,
                 int& num_vertices,
                 std::vector<IBTK::Point>& vertex_posn,
                 void* ctx)
{
    const MarkerParameters& params = *static_cast<const MarkerParameters*>(ctx);
    vertex_posn.clear();
    if (ln == params.finest_ln && params.probes)
    {
        // For each probe point, the points displaced by +EPS and -EPS cell widths along each axis.
        const double h = 1.0 / params.num_cells;
        for (const IBTK::Point& probe : generate_probe_points(params.num_cells))
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                for (const double sign : { 1.0, -1.0 })
                {
                    IBTK::Point X = probe;
                    X[d] += sign * EPS * h;
                    vertex_posn.push_back(X);
                }
            }
        }
    }
    else if (ln == params.finest_ln)
    {
        const double h = 1.0 / params.num_cells;
        const int num_along_line = 2 * params.num_cells;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            // Each other axis contributes a lower or an upper offset.
            for (int sides = 0; sides < (1 << (NDIM - 1)); ++sides)
            {
                for (int k = 0; k < num_along_line; ++k)
                {
                    IBTK::Point X;
                    unsigned int bit = 0;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (d == axis)
                        {
                            X[d] = (k + 0.5) / num_along_line;
                        }
                        else
                        {
                            const bool upper = (sides >> bit++) & 1;
                            X[d] = upper ? 1.0 - (0.35 + 0.05 * d) * h : (0.3 + 0.05 * d) * h;
                        }
                    }
                    vertex_posn.push_back(X);
                }
            }
        }
    }
    num_vertices = static_cast<int>(vertex_posn.size());
    return;
} // generate_markers

// Register and allocate velocity data on the level that has the ghost width of
// the IB integrator's velocity data, which is wider when its velocity boundary
// operator requires it, and return its patch data index.
int
allocate_velocity_data(Pointer<PatchLevel<NDIM>> level,
                       Pointer<INSHierarchyIntegrator> ins_integrator,
                       Pointer<IBHierarchyIntegrator> ib_integrator,
                       const double time)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<Variable<NDIM>> u_var = ins_integrator->getVelocityVariable();
    const int ib_u_idx =
        var_db->mapVariableAndContextToIndex(u_var, var_db->getContext(ib_integrator->getName() + "::IB"));
    const IntVector<NDIM> ghost_width = level->getPatchDescriptor()->getPatchDataFactory(ib_u_idx)->getGhostCellWidth();
    const int u_idx =
        var_db->registerVariableAndContext(u_var, var_db->getContext("interpolation_spreading_adjoint"), ghost_width);
    level->allocatePatchData(u_idx, time);
    return u_idx;
} // allocate_velocity_data

// With the IB integrator's velocity boundary conditions made homogeneous,
// velocity interpolation J and force spreading S must satisfy
// (F, J u) = (S F, u), where the Eulerian inner product counts each degree of
// freedom once and weights it by the cell volume.
void
check_interpolation_spreading_adjoint(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                                      Pointer<INSHierarchyIntegrator> ins_integrator,
                                      Pointer<IBHierarchyIntegrator> ib_integrator,
                                      Pointer<IBMethod> ib_method_ops)
{
    const int ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
    const double time = ib_integrator->getIntegratorTime();
    const int u_idx = allocate_velocity_data(level, ins_integrator, ib_integrator, time);
    const int f_idx = VariableDatabase<NDIM>::getDatabase()->registerClonedPatchDataIndex(
        ins_integrator->getVelocityVariable(), u_idx);
    level->allocatePatchData(f_idx, time);

    RobinPhysBdryPatchStrategy* bdry_op = ib_integrator->getVelocityPhysBdryOp();
    bdry_op->setPatchDataIndex(u_idx);
    bdry_op->setHomogeneousBc(true);
    Pointer<RefineAlgorithm<NDIM>> ghost_fill_alg = new RefineAlgorithm<NDIM>();
    ghost_fill_alg->registerRefine(u_idx, u_idx, u_idx, nullptr);
    std::vector<Pointer<RefineSchedule<NDIM>>> ghost_fill_scheds(ln + 1);
    ghost_fill_scheds[ln] = ghost_fill_alg->createSchedule(level, bdry_op);

    // Set an arbitrary velocity that depends only on position, and fill its
    // ghost values with the homogeneous boundary conditions.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const x_lower = pgeom->getXLower();
        const double* const dx = pgeom->getDx();
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        u_data->fillAll(0.0);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                double phase = 1.7 * axis;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    const double x = x_lower[d] + dx[d] * (b()(d) - patch_box.lower(d) + (d == axis ? 0.0 : 0.5));
                    phase += (3.1 + 2.2 * d) * x;
                }
                u_data->getArrayData(axis)(b(), 0) = std::sin(phase);
            }
        }
    }
    ghost_fill_scheds[ln]->fillData(time);

    // Interpolate u, and set an arbitrary force that depends only on the
    // Lagrangian index.
    LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    Pointer<LData> X_data = l_data_manager->getLData("X", ln);
    Pointer<LData> U_data = l_data_manager->createLData("interpolation_spreading_adjoint_U", ln, NDIM);
    Pointer<LData> F_data = l_data_manager->createLData("interpolation_spreading_adjoint_F", ln, NDIM);
    l_data_manager->interp(u_idx, U_data, X_data, ln, {}, ghost_fill_scheds, time);
    double F_dot_Ju = 0.0;
    {
        boost::multi_array_ref<double, 2>& U = *U_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& F = *F_data->getLocalFormVecArray();
        for (const auto& node : l_data_manager->getLMesh(ln)->getLocalNodes())
        {
            const int lag_idx = node->getLagrangianIndex();
            const int petsc_idx = node->getLocalPETScIndex();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F[petsc_idx][d] = std::cos(1.3 * lag_idx + 0.7 * d + 0.2);
                F_dot_Ju += F[petsc_idx][d] * U[petsc_idx][d];
            }
        }
        U_data->restoreArrays();
        F_data->restoreArrays();
    }

    // Spread F with the same boundary operator.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<SideData<NDIM, double>> f_data = level->getPatch(p())->getPatchData(f_idx);
        f_data->fillAll(0.0);
    }
    l_data_manager->spread(f_idx, F_data, X_data, bdry_op, ln, {}, time);
    double SF_dot_u = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        const double* const dx = pgeom->getDx();
        double cell_volume = 1.0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            cell_volume *= dx[d];
        }
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch_box, axis)); b; b++)
            {
                // Faces on the boundary of a patch that are not on a physical
                // boundary are shared with another patch.
                const bool on_lower_face = b()(axis) == patch_box.lower(axis);
                const bool on_upper_face = b()(axis) == patch_box.upper(axis) + 1;
                const bool shared = (on_lower_face && !pgeom->getTouchesRegularBoundary(axis, 0)) ||
                                    (on_upper_face && !pgeom->getTouchesRegularBoundary(axis, 1));
                const double weight = shared ? 0.5 * cell_volume : cell_volume;
                SF_dot_u += weight * f_data->getArrayData(axis)(b(), 0) * u_data->getArrayData(axis)(b(), 0);
            }
        }
    }
    level->deallocatePatchData(u_idx);
    level->deallocatePatchData(f_idx);
    pout << "patches on the level = " << level->getNumberOfPatches() << '\n';
    pout << std::setprecision(12) << "(F, J u) = " << IBTK_MPI::sumReduction(F_dot_Ju) << '\n'
         << "(S F, u) = " << IBTK_MPI::sumReduction(SF_dot_u) << '\n';
    return;
} // check_interpolation_spreading_adjoint

// A smooth potential, from which a velocity that is discretely divergence-free
// in the domain is derived: a stream function in 2D (only component 0), and the
// components of a vector potential in 3D. It vanishes on a boundary where it
// would give a velocity normal to the boundary, which is every boundary except
// the lower boundaries when they are open.
double
potential(const unsigned int component, const std::array<double, NDIM>& x, const bool open_lower_faces)
{
    double value = 1.0;
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        if (NDIM == 2 || d != component)
        {
            value *= open_lower_faces ? 1.0 - x[d] : x[d] * (1.0 - x[d]);
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
                     const bool open_lower_faces)
{
    std::array<double, NDIM> x_upper = x, x_lower = x;
    x_upper[d] += 0.5 * dx;
    x_lower[d] -= 0.5 * dx;
    return (potential(component, x_upper, open_lower_faces) - potential(component, x_lower, open_lower_faces)) / dx;
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
                              Pointer<Database> input_db)
{
    const int ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
    const double time = ib_integrator->getIntegratorTime();
    const int u_idx = allocate_velocity_data(level, ins_integrator, ib_integrator, time);
    const bool open_lower_faces = input_db->getBoolWithDefault("open_lower_faces", false);
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
                    sign * potential_difference(0, x, 1 - axis, dx[1 - axis], open_lower_faces);
#else
                const unsigned int next = (axis + 1) % 3, previous = (axis + 2) % 3;
                u_data->getArrayData(axis)(b(), 0) =
                    potential_difference(previous, x, next, dx[next], open_lower_faces) -
                    potential_difference(next, x, previous, dx[previous], open_lower_faces);
#endif
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
    constexpr double GRID_DIV_BOUND = 1.0e-9;
    if (!(max_grid_div < GRID_DIV_BOUND))
    {
        TBOX_ERROR("max |div u| of the velocity in the domain is " << max_grid_div << '\n');
    }
    pout << std::setprecision(3) << "max |div u| of the velocity in the domain = " << max_grid_div << '\n';

    // Interpolate u at the points around each probe point.
    LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    Pointer<LData> X_data = l_data_manager->getLData("X", ln);
    Pointer<LData> U_data = l_data_manager->createLData("interpolation_spreading_adjoint_U", ln, NDIM);
    l_data_manager->interp(u_idx, U_data, X_data, ln, {}, ghost_fill_scheds, time);
    const int num_probes = static_cast<int>(generate_probe_points(input_db->getInteger("N")).size());
    const int num_points = 2 * NDIM * num_probes;
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

    // The centered-difference divergence of the interpolated velocity at each probe point.
    const double eps_length = EPS / input_db->getInteger("N");
    double max_interpolated_div = 0.0;
    for (int probe = 0; probe < num_probes; ++probe)
    {
        double div = 0.0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            const int plus = 2 * NDIM * probe + 2 * d;
            div += (U[plus * NDIM + d] - U[(plus + 1) * NDIM + d]) / (2.0 * eps_length);
        }
        if (!std::isfinite(div))
        {
            TBOX_ERROR("the divergence of the interpolated velocity is not finite at probe point " << probe << '\n');
        }
        max_interpolated_div = std::max(max_interpolated_div, std::abs(div));
    }
    level->deallocatePatchData(u_idx);
    pout << "patches on the level = " << level->getNumberOfPatches() << '\n';
    pout << "probe points = " << num_probes << '\n';
    pout << std::setprecision(3) << "max |div| of the interpolated velocity = " << max_interpolated_div << '\n';
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

        MarkerParameters marker_params{ input_db->getInteger("MAX_LEVELS") - 1,
                                        input_db->getInteger("N"),
                                        input_db->keyExists("check_interpolated_divergence") };
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
        if (marker_params.probes)
        {
            check_interpolated_divergence(patch_hierarchy, ins_integrator, ib_integrator, ib_method_ops, input_db);
        }
        else
        {
            check_interpolation_spreading_adjoint(patch_hierarchy, ins_integrator, ib_integrator, ib_method_ops);
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
