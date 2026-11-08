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
// stress, and stops unless the printed sums and maxima equal the reference values of the input to 1e-12.

#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IMPInitializer.h>
#include <ibamr/IMPMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
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
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
// Print a quantity with enough digits to see a difference between partitions, and record whether it differs from the
// reference value of the input by more than 1e-12, relative to the larger of the reference value and 1.
bool
check_against_reference(Pointer<Database> reference_db, const std::string& name, const double value)
{
    const double reference = reference_db->getDouble(name);
    const bool equal = std::abs(value - reference) <= 1.0e-12 * std::max(1.0, std::abs(reference));
    pout << name << " = " << value << '\n';
    if (!equal) pout << "  differs from the reference value " << reference << '\n';
    return equal;
} // check_against_reference

void
check_partition_independence(Pointer<PatchHierarchy<NDIM>> patch_hierarchy,
                             Pointer<INSHierarchyIntegrator> ins_integrator,
                             Pointer<IBHierarchyIntegrator> ib_integrator,
                             Pointer<IMPMethod> ib_method_ops,
                             Pointer<Database> reference_db)
{
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

    // Interpolate u, and set a stress that depends only on the Lagrangian index.
    ib_method_ops->interpolateVelocity(u_idx, {}, ghost_fill_scheds, time);
    LDataManager* l_data_manager = ib_method_ops->getLDataManager();
    Pointer<LData> U_data = l_data_manager->getLData(LDataManager::VEL_DATA_NAME, ln);
    Pointer<LData> Grad_U_data = l_data_manager->getLData("Grad_U", ln);
    Pointer<LData> tau_data = l_data_manager->getLData("tau", ln);
    double num_points = 0.0, num_not_finite = 0.0;
    double U_sum_squares = 0.0, U_max = 0.0, Grad_U_sum_squares = 0.0, Grad_U_max = 0.0;
    {
        boost::multi_array_ref<double, 2>& U = *U_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& Grad_U = *Grad_U_data->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& tau = *tau_data->getLocalFormVecArray();
        for (const auto& node : l_data_manager->getLMesh(ln)->getLocalNodes())
        {
            const int lag_idx = node->getLagrangianIndex();
            const int petsc_idx = node->getLocalPETScIndex();
            num_points += 1.0;
            for (unsigned int i = 0; i < NDIM; ++i)
            {
                if (!std::isfinite(U[petsc_idx][i])) num_not_finite += 1.0;
                U_sum_squares += U[petsc_idx][i] * U[petsc_idx][i];
                U_max = std::max(U_max, std::abs(U[petsc_idx][i]));
                for (unsigned int j = 0; j < NDIM; ++j)
                {
                    const double g = Grad_U[petsc_idx][NDIM * i + j];
                    if (!std::isfinite(g)) num_not_finite += 1.0;
                    Grad_U_sum_squares += g * g;
                    Grad_U_max = std::max(Grad_U_max, std::abs(g));
                    tau[petsc_idx][NDIM * i + j] = std::cos(1.3 * lag_idx + 0.7 * (NDIM * i + j) + 0.2);
                }
            }
        }
        U_data->restoreArrays();
        Grad_U_data->restoreArrays();
        tau_data->restoreArrays();
    }

    pout << std::setprecision(14);
    bool all_equal = true;
    all_equal &= check_against_reference(reference_db, "points", IBTK_MPI::sumReduction(num_points));
    all_equal &= check_against_reference(reference_db, "not_finite", IBTK_MPI::sumReduction(num_not_finite));
    all_equal &= check_against_reference(reference_db, "U_sum_squares", IBTK_MPI::sumReduction(U_sum_squares));
    all_equal &= check_against_reference(reference_db, "U_max", IBTK_MPI::maxReduction(U_max));
    all_equal &=
        check_against_reference(reference_db, "Grad_U_sum_squares", IBTK_MPI::sumReduction(Grad_U_sum_squares));
    all_equal &= check_against_reference(reference_db, "Grad_U_max", IBTK_MPI::maxReduction(Grad_U_max));
    plog << std::flush;

    // Spread the stress. The spread force is zero outside the domain. Each face is counted once.
    ib_method_ops->spreadForce(f_idx, bdry_op, {}, time);
    double f_sum = 0.0, f_sum_squares = 0.0, f_max = 0.0;
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
                if (shared) continue;
                f_sum += v;
                f_sum_squares += v * v;
                f_max = std::max(f_max, std::abs(v));
            }
        }
    }
    level->deallocatePatchData(u_idx);
    level->deallocatePatchData(f_idx);

    all_equal &= check_against_reference(reference_db, "f_sum", IBTK_MPI::sumReduction(f_sum));
    all_equal &= check_against_reference(reference_db, "f_sum_squares", IBTK_MPI::sumReduction(f_sum_squares));
    all_equal &= check_against_reference(reference_db, "f_max", IBTK_MPI::maxReduction(f_max));
    if (!all_equal) TBOX_ERROR("the interpolated and spread quantities differ from the reference values\n");
    return;
} // check_partition_independence
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

        // A strip of elements along the lower boundary, across the boundary between two patches.
        Mesh mesh(init.comm(), NDIM);
        std::array<double, 4> strip;
        input_db->getDoubleArray("strip", strip.data(), 4);
        MeshTools::Generation::build_square(mesh,
                                            static_cast<int>(input_db->getInteger("num_elements_x")),
                                            static_cast<int>(input_db->getInteger("num_elements_y")),
                                            strip[0],
                                            strip[1],
                                            strip[2],
                                            strip[3],
                                            QUAD4);
        mesh.prepare_for_use();

        Pointer<INSHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        Pointer<IMPMethod> ib_method_ops =
            new IMPMethod("IMPMethod", app_initializer->getComponentDatabase("IMPMethod"));
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

        Pointer<IMPInitializer> ib_initializer =
            new IMPInitializer("IMPInitializer",
                               app_initializer->getComponentDatabase("IMPInitializer"),
                               patch_hierarchy,
                               gridding_algorithm);
        ib_initializer->registerMesh(&mesh);
        ib_method_ops->registerLInitStrategy(ib_initializer);

        std::vector<std::unique_ptr<muParserRobinBcCoefs>> u_bc_coef_owners;
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            u_bc_coef_owners.push_back(std::make_unique<muParserRobinBcCoefs>(
                "u_bc_coefs_" + std::to_string(d),
                app_initializer->getComponentDatabase("VelocityBcCoefs_" + std::to_string(d)),
                grid_geometry));
            u_bc_coefs.push_back(u_bc_coef_owners.back().get());
        }
        ins_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);

        ib_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();

        // Regrid, as advanceHierarchy() does at the initial time, to set up the Lagrangian indices in ghost cells.
        // Preprocessing sets up the velocity boundary conditions.
        ib_integrator->regridHierarchy();
        const double current_time = ib_integrator->getIntegratorTime();
        ib_integrator->preprocessIntegrateHierarchy(
            current_time, current_time + ib_integrator->getMaximumTimeStepSize(), 1);

        // Move the points by marker_shift grid cells along each axis, without assigning them to patches again.
        std::array<double, NDIM> marker_shift;
        input_db->getDoubleArray("marker_shift", marker_shift.data(), NDIM);
        const int ln = patch_hierarchy->getFinestLevelNumber();
        LDataManager* l_data_manager = ib_method_ops->getLDataManager();
        Pointer<LData> X_data = l_data_manager->getLData(LDataManager::POSN_DATA_NAME, ln);
        {
            boost::multi_array_ref<double, 2>& X = *X_data->getLocalFormVecArray();
            for (const auto& node : l_data_manager->getLMesh(ln)->getLocalNodes())
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    X[node->getLocalPETScIndex()][d] += marker_shift[d] * grid_geometry->getDx()[d];
                }
            }
            X_data->restoreArrays();
        }
        X_data->beginGhostUpdate();
        X_data->endGhostUpdate();

        check_partition_independence(
            patch_hierarchy, ins_integrator, ib_integrator, ib_method_ops, input_db->getDatabase("Reference"));
    }
} // main
