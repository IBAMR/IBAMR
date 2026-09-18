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

// IBMethod's finest-level LE coupling position accessor and the implicit
// Jacobian's interpolation matrix, built with IB kernel evaluators.
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/ib_kernel_evaluators.h>

#include <tbox/Logger.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>

#include <petscvec.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <array>
#include <cmath>
#include <fstream>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

namespace
{
constexpr int n_markers = 3;
const std::array<double, n_markers> marker_x = { 0.10, 0.47, 0.72 };
const std::array<double, n_markers> marker_y = { 0.20, 0.83, 0.35 };
// Distinct, deterministic per-entry velocities so that a transposed or
// mis-indexed accessor result is caught by the position check.
const std::array<double, NDIM* n_markers> marker_velocity = { 0.020, 0.025, 0.030, 0.035, 0.040, 0.045 };

void
generate_markers(const unsigned int& structure,
                 const int& level,
                 int& count,
                 std::vector<IBTK::Point>& positions,
                 void*)
{
    TBOX_ASSERT(structure == 0 && level == 0);
    count = n_markers;
    positions.resize(count);
    for (int i = 0; i < n_markers; ++i)
    {
        positions[i](0) = marker_x[i];
        positions[i](1) = marker_y[i];
    }
}

// Independent one-dimensional weight formulas: the symmetric radial form of
// the IB_4 kernel (rather than the moment/square-root recurrence used by
// IBKernelEvaluators::IB4), and the elementary linear-interpolation formula.
double
ib4_weight_1d(int offset, double r)
{
    const double x = std::abs(r - offset);
    if (x <= 1.0) return (3.0 - 2.0 * x + std::sqrt(1.0 + 4.0 * x - 4.0 * x * x)) / 8.0;
    return (5.0 - 2.0 * x - std::sqrt(-7.0 + 12.0 * x - 4.0 * x * x)) / 8.0;
}

double
linear_weight_1d(int offset, double r)
{
    return offset == 0 ? 1.0 - r : r;
}

// Check J's columns and weights against an independent calculation that
// reuses only the documented stencil-placement formula (not the production
// weight evaluators). width is the (isotropic) per-axis stencil width.
template <class Weight1D>
void
check_interp_matrix(Mat J,
                    Vec X_vec,
                    Pointer<SideData<NDIM, int>> dofs,
                    Pointer<PatchLevel<NDIM>> level,
                    int width,
                    Weight1D weight_1d,
                    int& column_mismatches,
                    double& max_weight_error)
{
    column_mismatches = 0;
    max_weight_error = 0.0;
    Pointer<CartesianGridGeometry<NDIM>> geom = level->getGridGeometry();
    const double* const x_lower = geom->getXLower();
    const double* const dx = geom->getDx();
    const hier::Index<NDIM>& domain_lower = level->getPhysicalDomain()[0].lower();

    PetscErrorCode ierr;
    PetscInt begin, end;
    ierr = MatGetOwnershipRange(J, &begin, &end);
    IBTK_CHKERRQ(ierr);
    const PetscScalar* X_arr;
    ierr = VecGetArrayRead(X_vec, &X_arr);
    IBTK_CHKERRQ(ierr);
    for (PetscInt row = begin; row < end; ++row)
    {
        const int point = static_cast<int>(row - begin) / NDIM;
        const int axis = static_cast<int>(row - begin) % NDIM;
        std::array<double, NDIM> q, r;
        std::array<int, NDIM> lower;
        for (int d = 0; d < NDIM; ++d)
        {
            const double centering_offset = d == axis ? 0.0 : 0.5;
            q[d] = (X_arr[NDIM * point + d] - x_lower[d]) / dx[d] + domain_lower(d) - centering_offset;
            // Documented placement formula for even width: normal axis uses
            // floor(q) - N/2 + 1; other axes use ceil(q) - N/2.
            lower[d] = d == axis ? static_cast<int>(std::floor(q[d])) - width / 2 + 1 :
                                   static_cast<int>(std::ceil(q[d])) - width / 2;
            r[d] = q[d] - lower[d];
        }
        std::map<PetscInt, double> expected;
        for (int j = 0; j < width; ++j)
            for (int i = 0; i < width; ++i)
            {
                SAMRAI::hier::Index<NDIM> index;
                index(0) = lower[0] + i;
                index(1) = lower[1] + j;
                const SideIndex<NDIM> side(index, axis, SideIndex<NDIM>::Lower);
                const int column = (*dofs)(side);
                expected[column] = weight_1d(i, r[0]) * weight_1d(j, r[1]);
            }
        PetscInt count;
        const PetscInt* columns;
        const PetscScalar* values;
        ierr = MatGetRow(J, row, &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
        if (count != static_cast<PetscInt>(expected.size())) ++column_mismatches;
        for (PetscInt k = 0; k < count; ++k)
        {
            const auto found = expected.find(columns[k]);
            if (found == expected.end())
            {
                ++column_mismatches;
            }
            else
            {
                max_weight_error = std::max(max_weight_error, std::abs(PetscRealPart(values[k]) - found->second));
            }
        }
        ierr = MatRestoreRow(J, row, &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(X_vec, &X_arr);
    IBTK_CHKERRQ(ierr);
}

// Build the fixture (geometry, hierarchy, IBMethod with probe markers, and a
// side-centered DOF index variable), run the two-step position/Jacobian
// checks, and accumulate results into the caller's diagnostics.
void
run_fixture(Pointer<AppInitializer> app,
            bool use_fixed_ops,
            double& max_position_error,
            int& column_mismatches,
            double& max_weight_error,
            bool trigger_bad_accessor_time = false)
{
    // Every named object below registers itself with a process-global
    // registry (RestartManager, VariableDatabase); since run_fixture is
    // called more than once in the same process, give each call's objects
    // distinct names.
    const std::string suffix = use_fixed_ops ? "_fixed" : "_unfixed";
    Pointer<IBMethod> method = new IBMethod("IBMethod" + suffix, app->getComponentDatabase("IBMethod"));
    method->setUseFixedLEOperators(use_fixed_ops);
    // preprocessIntegrateData() needs a registered IBHierarchyIntegrator (for
    // getStartTime()); the explicit integrator is otherwise unused here.
    Pointer<INSStaggeredHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator" + suffix, app->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    Pointer<IBExplicitHierarchyIntegrator> ib_integrator =
        new IBExplicitHierarchyIntegrator("IBExplicitHierarchyIntegrator" + suffix,
                                          app->getComponentDatabase("IBExplicitHierarchyIntegrator"),
                                          method,
                                          ins_integrator);
    Pointer<CartesianGridGeometry<NDIM>> geometry =
        new CartesianGridGeometry<NDIM>("CartesianGeometry" + suffix, app->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy" + suffix, geometry);
    Pointer<StandardTagAndInitialize<NDIM>> tagger = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize" + suffix, method, app->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> balancer =
        new LoadBalancer<NDIM>("LoadBalancer" + suffix, app->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
        "GriddingAlgorithm" + suffix, app->getComponentDatabase("GriddingAlgorithm"), tagger, boxes, balancer);
    Pointer<IBRedundantInitializer> initializer = new IBRedundantInitializer(
        "IBRedundantInitializer" + suffix, app->getComponentDatabase("IBRedundantInitializer"));
    initializer->setStructureNamesOnLevel(0, { "probes" });
    initializer->registerInitStructureFunction(generate_markers);
    method->registerLInitStrategy(initializer);
    gridding->makeCoarsestLevel(hierarchy, 0.0);
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    if (IBTK_MPI::getNodes() != 1 || level->getNumberOfPatches() != 1)
        TBOX_ERROR("implicit_stokes_ib_solver_components_01 requires one patch on one rank\n");

    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = variables->getContext("implicit_stokes_ib_solver_components_01" + suffix);
    Pointer<SideVariable<NDIM, double>> velocity = new SideVariable<NDIM, double>("velocity" + suffix);
    Pointer<SideVariable<NDIM, int>> indices = new SideVariable<NDIM, int>("indices" + suffix);
    const int u = variables->registerVariableAndContext(velocity, context, method->getMinimumGhostCellWidth());
    const int dof = variables->registerVariableAndContext(indices, context, method->getMinimumGhostCellWidth());
    level->allocatePatchData(u, 0.0);
    level->allocatePatchData(dof, 0.0);
    Pointer<SideData<NDIM, double>> u_data = level->getPatch(0)->getPatchData(u);
    u_data->fillAll(0.0);
    std::vector<int> counts;
    PETScVecUtilities::constructPatchLevelDOFIndices(counts, dof, level);
    Pointer<SideData<NDIM, int>> dofs = level->getPatch(0)->getPatchData(dof);
    std::vector<Pointer<CoarsenSchedule<NDIM>>> synch(1);
    std::vector<Pointer<RefineSchedule<NDIM>>> ghost_fill(1);
    method->initializePatchHierarchy(hierarchy, gridding, u, synch, ghost_fill, 0, 0.0, true);
    method->freeLInitStrategy();
    initializer.setNull();

    // Snapshot the initial marker positions once; every later "true" position
    // is then known in closed form as X0 + marker_velocity * time.
    std::vector<double> X0(NDIM * n_markers);
    {
        Vec X0_vec = method->getLDataManager()->getLData(LDataManager::POSN_DATA_NAME, 0)->getVec();
        const PetscScalar* arr;
        int ierr = VecGetArrayRead(X0_vec, &arr);
        IBTK_CHKERRQ(ierr);
        std::copy(arr, arr + NDIM * n_markers, X0.begin());
        ierr = VecRestoreArrayRead(X0_vec, &arr);
        IBTK_CHKERRQ(ierr);
    }

    max_position_error = 0.0;
    column_mismatches = 0;
    max_weight_error = 0.0;
    for (int step = 0; step < 2; ++step)
    {
        const double current_time = step * 0.125, new_time = current_time + 0.125, half_time = current_time + 0.0625;
        method->preprocessIntegrateData(current_time, new_time, 1);

        Vec X_new_vec = nullptr;
        method->createSolverVecs(&X_new_vec, nullptr);
        {
            PetscScalar* arr;
            int ierr = VecGetArray(X_new_vec, &arr);
            IBTK_CHKERRQ(ierr);
            for (int k = 0; k < NDIM * n_markers; ++k) arr[k] = X0[k] + marker_velocity[k] * new_time;
            ierr = VecRestoreArray(X_new_vec, &arr);
            IBTK_CHKERRQ(ierr);
        }
        method->setUpdatedPosition(X_new_vec);
        method->updateFixedLEOperators();

        if (trigger_bad_accessor_time)
        {
            // Neither the current, half, nor new time: the accessor must
            // reject this rather than dereference an unset pointer.
            method->getFinestLevelLECouplingPositions(current_time + 0.03);
        }

        Vec new_vec = method->getFinestLevelLECouplingPositions(new_time);
        Vec half_vec = method->getFinestLevelLECouplingPositions(half_time);
        const PetscScalar *new_arr, *half_arr;
        int ierr = VecGetArrayRead(new_vec, &new_arr);
        IBTK_CHKERRQ(ierr);
        ierr = VecGetArrayRead(half_vec, &half_arr);
        IBTK_CHKERRQ(ierr);
        for (int k = 0; k < NDIM * n_markers; ++k)
        {
            const double expected_new = X0[k] + marker_velocity[k] * new_time;
            const double expected_half = X0[k] + marker_velocity[k] * half_time;
            max_position_error = std::max(max_position_error, std::abs(PetscRealPart(new_arr[k]) - expected_new));
            max_position_error = std::max(max_position_error, std::abs(PetscRealPart(half_arr[k]) - expected_half));
        }
        ierr = VecRestoreArrayRead(new_vec, &new_arr);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreArrayRead(half_vec, &half_arr);
        IBTK_CHKERRQ(ierr);

        // The integrator's two jacobian_delta_fcn kernels.
        Mat J_ib4 = nullptr;
        PETScMatUtilities::constructPatchLevelSCInterpOp(
            J_ib4,
            IBTK::IBKernelEvaluatorTensorProduct{ IBTK::IBKernelEvaluators::IB4{} },
            new_vec,
            counts,
            dof,
            level);
        int ib4_mismatches;
        double ib4_error;
        check_interp_matrix(J_ib4, new_vec, dofs, level, 4, ib4_weight_1d, ib4_mismatches, ib4_error);
        column_mismatches += ib4_mismatches;
        max_weight_error = std::max(max_weight_error, ib4_error);
        ierr = MatDestroy(&J_ib4);
        IBTK_CHKERRQ(ierr);

        Mat J_linear = nullptr;
        PETScMatUtilities::constructPatchLevelSCInterpOp(
            J_linear,
            IBTK::IBKernelEvaluatorTensorProduct{ IBTK::IBKernelEvaluators::BSpline<2>{} },
            new_vec,
            counts,
            dof,
            level);
        int linear_mismatches;
        double linear_error;
        check_interp_matrix(J_linear, new_vec, dofs, level, 2, linear_weight_1d, linear_mismatches, linear_error);
        column_mismatches += linear_mismatches;
        max_weight_error = std::max(max_weight_error, linear_error);
        ierr = MatDestroy(&J_linear);
        IBTK_CHKERRQ(ierr);

        ierr = VecDestroy(&X_new_vec);
        IBTK_CHKERRQ(ierr);
        method->postprocessIntegrateData(current_time, new_time, 1);
    }
    method->postprocessData();
    level->deallocatePatchData(u);
    level->deallocatePatchData(dof);
    variables->removePatchDataIndex(u);
    variables->removePatchDataIndex(dof);
}

// Construct an IBImplicitStaggeredHierarchyIntegrator with the given
// jacobian_delta_fcn value; used only to check construction-time validation.
void
try_construct_implicit_integrator(Pointer<AppInitializer> app)
{
    Pointer<INSStaggeredHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator", app->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    Pointer<IBMethod> method = new IBMethod("IBMethod", app->getComponentDatabase("IBMethod"));
    Pointer<IBImplicitStaggeredHierarchyIntegrator> integrator =
        new IBImplicitStaggeredHierarchyIntegrator("IBImplicitStaggeredHierarchyIntegrator",
                                                   app->getComponentDatabase("IBImplicitStaggeredHierarchyIntegrator"),
                                                   method,
                                                   ins_integrator,
                                                   false);
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "IB.log");
    const std::string input_file = argc > 1 ? argv[1] : "";

    SAMRAI::tbox::Pointer<SAMRAI::tbox::Logger::Appender> abort_appender = new TestAppender();
    SAMRAI::tbox::Logger::getInstance()->setAbortAppender(abort_appender);

    if (input_file.find("bad_jacobian_delta_fcn") != std::string::npos)
    {
        SAMRAI::tbox::PIO::logOnlyNodeZero("output");
        try_construct_implicit_integrator(app);
        return 0;
    }
    if (input_file.find("bad_accessor_time") != std::string::npos)
    {
        SAMRAI::tbox::PIO::logOnlyNodeZero("output");
        double max_position_error;
        int column_mismatches;
        double max_weight_error;
        run_fixture(app,
                    /*use_fixed_ops*/ true,
                    max_position_error,
                    column_mismatches,
                    max_weight_error,
                    /*trigger_bad_accessor_time*/ true);
        return 0;
    }

    double overall_max_position_error = 0.0;
    int overall_column_mismatches = 0;
    double overall_max_weight_error = 0.0;
    for (bool use_fixed_ops : { false, true })
    {
        double max_position_error;
        int column_mismatches;
        double max_weight_error;
        run_fixture(app, use_fixed_ops, max_position_error, column_mismatches, max_weight_error);
        overall_max_position_error = std::max(overall_max_position_error, max_position_error);
        overall_column_mismatches += column_mismatches;
        overall_max_weight_error = std::max(overall_max_weight_error, max_weight_error);
    }
    if (!(overall_max_position_error <= 1.0e-10))
    {
        TBOX_ERROR("Failed check: !(overall_max_position_error <= 1.0e-10).\n");
    }
    if (overall_column_mismatches != 0)
    {
        TBOX_ERROR("Failed check: overall_column_mismatches != 0.\n");
    }
    if (!(overall_max_weight_error <= 1.0e-10))
    {
        TBOX_ERROR("Failed check: !(overall_max_weight_error <= 1.0e-10).\n");
    }

    if (IBTK_MPI::getRank() == 0)
    {
        std::ofstream out("output");
        out << "max_position_error = " << overall_max_position_error << '\n';
        out << "column_mismatches = " << overall_column_mismatches << '\n';
        out << "max_weight_error = " << overall_max_weight_error << '\n';
    }
    return 0;
}
