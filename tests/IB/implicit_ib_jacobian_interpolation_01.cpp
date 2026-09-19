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
#include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBKernelConcepts.h>
#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/IBOperatorBuilder.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>

#include <tbox/Logger.h>
#include <tbox/MemoryDatabase.h>
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

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <fstream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

namespace
{
constexpr int N_MARKERS = 3;
constexpr std::size_t WIDE_TENT_WIDTH = 8;
constexpr std::array<double, N_MARKERS> MARKER_X = { 0.10, 0.47, 0.72 };
constexpr std::array<double, N_MARKERS> MARKER_Y = { 0.20, 0.83, 0.35 };
// Distinct per-entry velocities, so that a transposed or mis-indexed accessor
// result is caught by the position check.
constexpr std::array<double, NDIM* N_MARKERS> MARKER_VELOCITY = { 0.020, 0.025, 0.030, 0.035, 0.040, 0.045 };

void
generate_markers(const unsigned int& structure,
                 const int& level,
                 int& count,
                 std::vector<IBTK::Point>& positions,
                 void*)
{
    TBOX_ASSERT(structure == 0 && level == 0);
    count = N_MARKERS;
    positions.resize(count);
    for (int i = 0; i < N_MARKERS; ++i)
    {
        positions[i](0) = MARKER_X[i];
        positions[i](1) = MARKER_Y[i];
    }
}

// std::max keeps its first argument when the second is NaN, so a NaN
// difference must make the accumulated error infinite instead.
void
accumulate_error(double& maximum, const double error)
{
    maximum = std::isfinite(error) ? std::max(maximum, error) : std::numeric_limits<double>::infinity();
}

// One-dimensional kernels as functions of the distance from the evaluation
// point to a stencil point, in grid spacings. The IB_4 kernel is written in its
// symmetric radial form rather than the form used by IBKernelEvaluators::IB4.
double
ib4_weight(const double distance)
{
    const double x = std::abs(distance);
    if (x <= 1.0) return (3.0 - 2.0 * x + std::sqrt(1.0 + 4.0 * x - 4.0 * x * x)) / 8.0;
    return (5.0 - 2.0 * x - std::sqrt(-7.0 + 12.0 * x - 4.0 * x * x)) / 8.0;
}

double
bspline2_weight(const double distance)
{
    return std::max(0.0, 1.0 - std::abs(distance));
}

double
bspline3_weight(const double distance)
{
    const double x = std::abs(distance);
    return x < 0.5 ? 0.75 - x * x : (x < 1.5 ? 0.5 * (1.5 - x) * (1.5 - x) : 0.0);
}

double
bspline4_weight(const double distance)
{
    const double x = std::abs(distance);
    return x < 1.0 ? 2.0 / 3.0 - x * x + 0.5 * x * x * x : (x < 2.0 ? (2.0 - x) * (2.0 - x) * (2.0 - x) / 6.0 : 0.0);
}

// A kernel that the library does not define: a tent of half-width 3/2.
double
tent_weight(const double distance)
{
    return std::max(0.0, 1.0 - std::abs(distance) / 1.5);
}

template <std::size_t Width>
struct TentKernel
{
    static constexpr std::size_t get_stencil_width()
    {
        return Width;
    }

    template <IBTK::IBKernelWeights Output, std::floating_point Input>
    Output evaluate(const Input r) const
    {
        Output w{};
        for (std::size_t i = 0; i < Width; ++i) w[i] = tent_weight(static_cast<double>(r) - static_cast<double>(i));
        return w;
    }
};

struct KernelCase
{
    std::string name;
    int normal_width, transverse_width;
    double (*normal_weight)(double);
    double (*transverse_weight)(double);
};

// The kernels selected by jacobian_delta_fcn, and the width and weights that
// they must produce, in the normal and transverse directions.
const std::vector<KernelCase> BUILT_IN_CASES = { { "IB_4", 4, 4, ib4_weight, ib4_weight },
                                                 { "PIECEWISE_LINEAR", 2, 2, bspline2_weight, bspline2_weight },
                                                 { "BSPLINE_3", 3, 3, bspline3_weight, bspline3_weight },
                                                 { "BSPLINE_4", 4, 4, bspline4_weight, bspline4_weight },
                                                 { "COMPOSITE_BSPLINE_2_3", 2, 3, bspline2_weight, bspline3_weight },
                                                 { "COMPOSITE_BSPLINE_3_4", 3, 4, bspline3_weight, bspline4_weight },
                                                 { "COMPOSITE_BSPLINE_4_3", 4, 3, bspline4_weight, bspline3_weight } };

// First stencil index in one direction, for grid coordinate q that includes
// the data centering.
int
first_stencil_index(const double q, const int width, const bool normal)
{
    if (width % 2 != 0) return static_cast<int>(std::floor(q + 0.5)) - (width - 1) / 2;
    return normal ? static_cast<int>(std::floor(q)) - width / 2 + 1 : static_cast<int>(std::ceil(q)) - width / 2;
}

// Check J's columns and weights against the stencil placement documented for
// PETScMatUtilities::constructPatchLevelSCInterpOp() and the one-dimensional
// weights of kernel.
void
check_interp_matrix(Mat J,
                    Vec X_vec,
                    Pointer<SideData<NDIM, int>> dofs,
                    Pointer<PatchLevel<NDIM>> level,
                    const KernelCase& kernel,
                    int& column_mismatches,
                    double& max_weight_error)
{
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
        std::array<double, NDIM> q;
        std::array<int, NDIM> lower, width;
        for (int d = 0; d < NDIM; ++d)
        {
            const bool normal = d == axis;
            q[d] = (X_arr[NDIM * point + d] - x_lower[d]) / dx[d] + domain_lower(d) - (normal ? 0.0 : 0.5);
            width[d] = normal ? kernel.normal_width : kernel.transverse_width;
            lower[d] = first_stencil_index(q[d], width[d], normal);
        }
        std::map<PetscInt, double> expected;
        for (int j = 0; j < width[1]; ++j)
            for (int i = 0; i < width[0]; ++i)
            {
                SAMRAI::hier::Index<NDIM> index;
                index(0) = lower[0] + i;
                index(1) = lower[1] + j;
                const SideIndex<NDIM> side(index, axis, SideIndex<NDIM>::Lower);
                const auto weight_0 = axis == 0 ? kernel.normal_weight : kernel.transverse_weight;
                const auto weight_1 = axis == 1 ? kernel.normal_weight : kernel.transverse_weight;
                expected[(*dofs)(side)] = weight_0(q[0] - index(0)) * weight_1(q[1] - index(1));
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
                accumulate_error(max_weight_error, std::abs(PetscRealPart(values[k]) - found->second));
            }
        }
        ierr = MatRestoreRow(J, row, &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(X_vec, &X_arr);
    IBTK_CHKERRQ(ierr);
}

// The builder that an IBImplicitStaggeredHierarchyIntegrator constructed with
// jacobian_delta_fcn = kernel_name uses, after the optional application override.
IBTK::IBOperatorBuilder
make_jacobian_builder(Pointer<AppInitializer> app,
                      const std::string& tag,
                      const std::string& kernel_name,
                      const std::optional<IBTK::IBOperatorBuilder>& application_builder = std::nullopt)
{
    Pointer<INSStaggeredHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator" + tag, app->getComponentDatabase("INSStaggeredHierarchyIntegrator"), false);
    Pointer<IBMethod> method = new IBMethod("IBMethod" + tag, app->getComponentDatabase("IBMethod"), false);
    Pointer<Database> input_db = new MemoryDatabase("IBImplicitStaggeredHierarchyIntegrator" + tag);
    input_db->putBool("eliminate_eulerian_vars", true);
    input_db->putString("jacobian_delta_fcn", kernel_name);
    Pointer<IBImplicitStaggeredHierarchyIntegrator> integrator = new IBImplicitStaggeredHierarchyIntegrator(
        "IBImplicitStaggeredHierarchyIntegrator" + tag, input_db, method, ins_integrator, false);
    if (application_builder) integrator->setJacobianOperatorBuilder(*application_builder);
    return integrator->getJacobianOperatorBuilder();
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
            int& dof_ghost_mismatches,
            bool trigger_bad_accessor_time = false)
{
    // Every named object below registers itself with a process-global
    // registry (RestartManager, VariableDatabase); since run_fixture is
    // called more than once in the same process, give each call's objects
    // distinct names.
    const std::string suffix = use_fixed_ops ? "_fixed" : "_unfixed";
    Pointer<IBMethod> method = new IBMethod("IBMethod" + suffix, app->getComponentDatabase("IBMethod"));
    Pointer<INSStaggeredHierarchyIntegrator> ins_integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator" + suffix, app->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    // The integrator registers itself with method, which preprocessIntegrateData()
    // needs (for getStartTime()). One run uses the kernel of the input, and
    // the other an application kernel that is wider than the strategy's.
    Pointer<Database> integrator_db = new MemoryDatabase("IBImplicitStaggeredHierarchyIntegrator" + suffix);
    integrator_db->putBool("eliminate_eulerian_vars", true);
    if (use_fixed_ops) integrator_db->putString("jacobian_delta_fcn", "APPLICATION_KERNEL");
    Pointer<IBImplicitStaggeredHierarchyIntegrator> integrator = new IBImplicitStaggeredHierarchyIntegrator(
        "IBImplicitStaggeredHierarchyIntegrator" + suffix, integrator_db, method, ins_integrator, false);
    if (use_fixed_ops)
    {
        integrator->setJacobianOperatorBuilder(
            IBTK::IBOperatorBuilder(IBTK::IBKernelEvaluatorTensorProduct{ TentKernel<WIDE_TENT_WIDTH>{} }));
    }
    method->setUseFixedLEOperators(use_fixed_ops);
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
    // The integrator registers its variables in the process-global variable
    // database, so only one run can initialize one.
    if (use_fixed_ops)
    {
        integrator->initializeHierarchyIntegrator(hierarchy, gridding);
        // The integrator's DOF index data must cover the Jacobian's stencil as
        // well as the strategy's.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<hier::Variable<NDIM>> u_dof_var = var_db->getVariable(integrator->getName() + "::u_dof_index");
        const int u_dof_idx = var_db->mapVariableAndContextToIndex(u_dof_var, integrator->getScratchContext());
        const IntVector<NDIM> ghosts =
            var_db->getPatchDescriptor()->getPatchDataFactory(u_dof_idx)->getGhostCellWidth();
        IntVector<NDIM> expected = method->getMinimumGhostCellWidth();
        expected.max(IntVector<NDIM>(static_cast<int>((WIDE_TENT_WIDTH + 1) / 2 + 1)));
        dof_ghost_mismatches += !(ghosts == expected);
        dof_ghost_mismatches += !(ghosts.max() > method->getMinimumGhostCellWidth().max());
    }
    gridding->makeCoarsestLevel(hierarchy, 0.0);
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    if (IBTK_MPI::getNodes() != 1 || level->getNumberOfPatches() != 1)
        TBOX_ERROR("implicit_ib_jacobian_interpolation_01 requires one patch on one rank\n");

    // The Jacobian's kernels, each with the builder that an integrator uses for it.
    // The DOF index data must cover the stencil of every one of them.
    std::vector<std::pair<KernelCase, IBTK::IBOperatorBuilder>> jacobian_kernels;
    for (std::size_t k = 0; k < BUILT_IN_CASES.size(); ++k)
    {
        jacobian_kernels.emplace_back(
            BUILT_IN_CASES[k],
            make_jacobian_builder(app, suffix + "_built_in_" + std::to_string(k), BUILT_IN_CASES[k].name));
    }
    jacobian_kernels.emplace_back(
        KernelCase{ "APPLICATION_KERNEL", 4, 4, tent_weight, tent_weight },
        make_jacobian_builder(app,
                              suffix + "_application",
                              "APPLICATION_KERNEL",
                              IBTK::IBOperatorBuilder(IBTK::IBKernelEvaluatorTensorProduct{ TentKernel<4>{} })));
    IntVector<NDIM> dof_ghosts = method->getMinimumGhostCellWidth();
    for (const auto& [kernel, builder] : jacobian_kernels)
    {
        dof_ghosts.max(IntVector<NDIM>(builder.getMinimumGhostWidth()));
    }

    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = variables->getContext("implicit_ib_jacobian_interpolation_01" + suffix);
    Pointer<SideVariable<NDIM, double>> velocity = new SideVariable<NDIM, double>("velocity" + suffix);
    Pointer<SideVariable<NDIM, int>> indices = new SideVariable<NDIM, int>("indices" + suffix);
    const int u = variables->registerVariableAndContext(velocity, context, method->getMinimumGhostCellWidth());
    const int dof = variables->registerVariableAndContext(indices, context, dof_ghosts);
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
    // is then known in closed form as X0 + MARKER_VELOCITY * time.
    std::vector<double> X0(NDIM * N_MARKERS);
    {
        Vec X0_vec = method->getLDataManager()->getLData(LDataManager::POSN_DATA_NAME, 0)->getVec();
        const PetscScalar* arr;
        int ierr = VecGetArrayRead(X0_vec, &arr);
        IBTK_CHKERRQ(ierr);
        std::copy(arr, arr + NDIM * N_MARKERS, X0.begin());
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
            for (int k = 0; k < NDIM * N_MARKERS; ++k) arr[k] = X0[k] + MARKER_VELOCITY[k] * new_time;
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
        for (int k = 0; k < NDIM * N_MARKERS; ++k)
        {
            const double expected_new = X0[k] + MARKER_VELOCITY[k] * new_time;
            const double expected_half = X0[k] + MARKER_VELOCITY[k] * half_time;
            accumulate_error(max_position_error, std::abs(PetscRealPart(new_arr[k]) - expected_new));
            accumulate_error(max_position_error, std::abs(PetscRealPart(half_arr[k]) - expected_half));
        }
        ierr = VecRestoreArrayRead(new_vec, &new_arr);
        IBTK_CHKERRQ(ierr);
        ierr = VecRestoreArrayRead(half_vec, &half_arr);
        IBTK_CHKERRQ(ierr);

        for (const auto& [kernel, builder] : jacobian_kernels)
        {
            Mat J = nullptr;
            builder.constructInterpolationMatrixSide(J, new_vec, counts, dof, level);
            int mismatches = 0;
            double weight_error = 0.0;
            check_interp_matrix(J, new_vec, dofs, level, kernel, mismatches, weight_error);
            column_mismatches += mismatches;
            accumulate_error(max_weight_error, weight_error);
            ierr = MatDestroy(&J);
            IBTK_CHKERRQ(ierr);
        }

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

// Construct an IBImplicitStaggeredHierarchyIntegrator from the input database
// and, if requested, ask it for the Jacobian's builder.
void
try_construct_implicit_integrator(Pointer<AppInitializer> app, const bool query_builder)
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
    if (query_builder) integrator->getJacobianOperatorBuilder();
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
        try_construct_implicit_integrator(app, false);
        return 0;
    }
    if (input_file.find("unset_jacobian_builder") != std::string::npos)
    {
        SAMRAI::tbox::PIO::logOnlyNodeZero("output");
        try_construct_implicit_integrator(app, true);
        return 0;
    }
    if (input_file.find("bad_accessor_time") != std::string::npos)
    {
        SAMRAI::tbox::PIO::logOnlyNodeZero("output");
        double max_position_error;
        int column_mismatches;
        double max_weight_error;
        int dof_ghost_mismatches = 0;
        run_fixture(app,
                    /*use_fixed_ops*/ true,
                    max_position_error,
                    column_mismatches,
                    max_weight_error,
                    dof_ghost_mismatches,
                    /*trigger_bad_accessor_time*/ true);
        return 0;
    }

    // Built-in kernels beyond the ones with independent weights below.
    int names_accepted = 0;
    for (const char* name : { "IB_3", "IB_5", "IB_6", "BSPLINE_8", "COMPOSITE_BSPLINE_7_8", "DISCONTINUOUS_LINEAR" })
    {
        make_jacobian_builder(app, std::string("_name_") + name, name);
        ++names_accepted;
    }

    double overall_max_position_error = 0.0;
    int overall_column_mismatches = 0;
    int overall_dof_ghost_mismatches = 0;
    double overall_max_weight_error = 0.0;
    for (bool use_fixed_ops : { false, true })
    {
        double max_position_error;
        int column_mismatches;
        double max_weight_error;
        run_fixture(
            app, use_fixed_ops, max_position_error, column_mismatches, max_weight_error, overall_dof_ghost_mismatches);
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
    if (overall_dof_ghost_mismatches != 0)
    {
        TBOX_ERROR("Failed check: overall_dof_ghost_mismatches != 0.\n");
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
        out << "integrator dof ghost width mismatches = " << overall_dof_ghost_mismatches << '\n';
        out << "additional jacobian_delta_fcn names accepted = " << names_accepted << '\n';
    }
    return 0;
}
