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
// Before the check, the markers are moved, without assigning them to patches
// again, by marker_displacement grid cells along each axis (see
// displace_markers).

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
#include <ibtk/ibtk_utilities.h>
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

// Place markers on lines parallel to each coordinate axis, a fraction of a grid
// cell away from each pair of boundaries that meet along an edge of the unit
// domain, so that the IB kernels overlap the physical boundaries, their
// corners (and edges), and the boundaries between patches.
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

// Markers may move up to one grid cell between the times at which they are
// assigned to patches, and then lie outside the patches that hold them. Move
// each marker by the given number of grid cells along each axis toward the
// center of the domain, without assigning the markers to patches again, and
// update the ghost values of the positions.
void
displace_markers(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                 Pointer<IBMethod> ib_method_ops,
                 const double displacement)
{
    const int ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<CartesianGridGeometry<NDIM>> grid_geom = patch_hierarchy->getGridGeometry();
    const IntVector<NDIM>& ratio = patch_hierarchy->getPatchLevel(ln)->getRatio();
    const double* const x_lower = grid_geom->getXLower();
    const double* const x_upper = grid_geom->getXUpper();
    const double* const dx_coarsest = grid_geom->getDx();
    LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    Pointer<LData> X_data = l_data_manager->getLData("X", ln);
    {
        boost::multi_array_ref<double, 2>& X = *X_data->getLocalFormVecArray();
        for (const auto& node : l_data_manager->getLMesh(ln)->getLocalNodes())
        {
            const int petsc_idx = node->getLocalPETScIndex();
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const double center = 0.5 * (x_lower[d] + x_upper[d]);
                const double h = dx_coarsest[d] / ratio(d);
                X[petsc_idx][d] += (X[petsc_idx][d] < center ? 1.0 : -1.0) * displacement * h;
            }
        }
        X_data->restoreArrays();
    }
    X_data->beginGhostUpdate();
    X_data->endGhostUpdate();
    return;
} // displace_markers

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
                                      Pointer<IBMethod> ib_method_ops,
                                      const bool divergence_free_extension)
{
    const int ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
    const double time = ib_integrator->getIntegratorTime();
    const int u_idx =
        allocate_velocity_data(level, ins_integrator, ib_integrator, ib_method_ops, divergence_free_extension, time);
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
    F_dot_Ju = IBTK_MPI::sumReduction(F_dot_Ju);
    SF_dot_u = IBTK_MPI::sumReduction(SF_dot_u);
    if (!IBTK::rel_equal_eps(F_dot_Ju, SF_dot_u, 1.0e-12))
    {
        TBOX_ERROR("(F, J u) = " << F_dot_Ju << " and (S F, u) = " << SF_dot_u << " do not agree\n");
    }
    pout << std::setprecision(12) << "(F, J u) = " << F_dot_Ju << '\n' << "(S F, u) = " << SF_dot_u << '\n';
    return;
} // check_interpolation_spreading_adjoint
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // Suppress warnings caused by running without Silo and by allowing patches smaller than the ghost cell width.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);

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
        const double marker_displacement = input_db->getDouble("marker_displacement");
        if (!(std::abs(marker_displacement) < 1.0))
        {
            TBOX_ERROR("marker_displacement = " << marker_displacement
                                                << " must be less than one grid cell in magnitude\n");
        }
        displace_markers(patch_hierarchy, ib_method_ops, marker_displacement);
        check_interpolation_spreading_adjoint(
            patch_hierarchy, ins_integrator, ib_integrator, ib_method_ops, divergence_free_extension);

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            delete u_bc_coefs[d];
        }
    }
} // main
