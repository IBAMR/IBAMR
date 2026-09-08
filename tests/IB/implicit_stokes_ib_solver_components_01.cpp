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

// Direct interpolation contracts, including the standalone IBMethod lifecycle.
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBKernelEvaluators.h>
#include <ibtk/IBKernelTensorProductEvaluator.h>
#include <ibtk/IBOperatorRegistry.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>

#include <tbox/Logger.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <SideData.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <vector>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

namespace
{
constexpr double epsilon = 1.0 / 1024.0;
// Grid coordinates just below, at, and above cell- and side-centering ties.
const std::array<double, 6> probe = { 4.0 - epsilon, 4.0, 4.0 + epsilon, 4.5 - epsilon, 4.5, 4.5 + epsilon };
// Explicit expected indices, independent of the production rounding formula.
const std::array<int, 6> side_nearest = { 4, 4, 4, 4, 5, 5 };
const std::array<int, 6> cell_nearest = { 3, 4, 4, 4, 4, 4 };
const std::array<int, 6> side_even_lower = { 3, 4, 4, 4, 4, 4 };
const std::array<int, 6> cell_even_lower = { 3, 3, 3, 3, 3, 4 };

// Deliberately asymmetric, distance-dependent weights expose reversal and
// incorrect lower-stencil coordinates as well as misplaced columns.
// Application-defined, move-only evaluator state must survive registration.
struct ProbeEvaluator
{
    std::unique_ptr<double> slope = std::make_unique<double>(1.0 / 32.0);

    std::array<double, 3> operator()(double r) const
    {
        return { 0.125 + r * *slope, 0.375 - r * *slope, 0.5 };
    }
};

template <class Evaluator, std::size_t N>
double
sample_error(const Evaluator& evaluator, double r, const std::array<double, N>& expected)
{
    const auto weights = evaluator(r);
    static_assert(std::tuple_size<decltype(weights)>::value == N, "Natural stencil size changed");
    double error = 0.0;
    for (std::size_t i = 0; i < N; ++i)
    {
        if (!std::isfinite(weights[i])) TBOX_ERROR("Nonfinite kernel weight\n");
        error = std::max(error, std::abs(weights[i] - expected[i]));
    }
    return error;
}

template <class Evaluator>
double
moment_error(const Evaluator& evaluator)
{
    constexpr int width = std::tuple_size<decltype(evaluator(0.0))>::value;
    double error = 0.0;
    for (int k = 0; k < 64; ++k)
    {
        const double r = 0.5 * width - 1.0 + k / 64.0;
        const auto w = evaluator(r);
        double sum = 0.0, moment = 0.0;
        for (int i = 0; i < width; ++i)
        {
            sum += w[i];
            moment += i * w[i];
        }
        if (!std::isfinite(sum) || !std::isfinite(moment)) TBOX_ERROR("Nonfinite kernel moment\n");
        error = std::max({ error, std::abs(sum - 1.0), std::abs(moment - r) });
    }
    return error;
}

int
check_kernels()
{
    const double a = (2.0 - std::sqrt(2.0)) / 8.0, b = (2.0 + std::sqrt(2.0)) / 8.0;
    const double K6 = (59.0 - std::sqrt(261.0)) / 60.0;
    double error = std::max(
        { sample_error(IBKernelEvaluatorBSpline<1>{}, -0.25, std::array<double, 1>{ 1.0 }),
          sample_error(IBKernelEvaluatorBSpline<2>{}, 0.25, std::array<double, 2>{ 0.75, 0.25 }),
          sample_error(IBKernelEvaluatorBSpline<3>{}, 1.0, std::array<double, 3>{ 0.125, 0.75, 0.125 }),
          sample_error(
              IBKernelEvaluatorBSpline<4>{}, 1.5, std::array<double, 4>{ 1.0 / 48, 23.0 / 48, 23.0 / 48, 1.0 / 48 }),
          sample_error(
              IBKernelEvaluatorBSpline<5>{}, 1.5, std::array<double, 5>{ 1.0 / 24, 11.0 / 24, 11.0 / 24, 1.0 / 24, 0 }),
          sample_error(IBKernelEvaluatorBSpline<6>{},
                       2.5,
                       std::array<double, 6>{
                           1.0 / 3840, 237.0 / 3840, 1682.0 / 3840, 1682.0 / 3840, 237.0 / 3840, 1.0 / 3840 }),
          sample_error(IBKernelEvaluatorIB3{}, 1.0, std::array<double, 3>{ 1.0 / 6, 2.0 / 3, 1.0 / 6 }),
          sample_error(IBKernelEvaluatorIB4{}, 1.5, std::array<double, 4>{ a, b, b, a }),
          sample_error(
              IBKernelEvaluatorIB5{},
              1.5,
              std::array<double, 5>{
                  0.0612224005711746881, 0.438777599428825312, 0.438777599428825312, 0.0612224005711746881, 0 }),
          sample_error(
              IBKernelEvaluatorIB6{},
              3.0,
              std::array<double, 6>{ 0, -1.0 / 16 + K6 / 8, 0.25, 5.0 / 8 - K6 / 4, 0.25, -1.0 / 16 + K6 / 8 }) });
    // Exact rational values from the truncated-power definition at r = 11/4.
    error = std::max(error,
                     sample_error(IBKernelEvaluatorBSpline<7>{},
                                  2.75,
                                  std::array<double, 7>{ 729.0 / 2949120,
                                                         112546.0 / 2949120,
                                                         963327.0 / 2949120,
                                                         1434812.0 / 2949120,
                                                         422087.0 / 2949120,
                                                         15618.0 / 2949120,
                                                         1.0 / 2949120 }));
    // Exact rational values from the truncated-power definition at r = 7/2.
    error = std::max(error,
                     sample_error(IBKernelEvaluatorBSpline<8>{},
                                  3.5,
                                  std::array<double, 8>{ 1.0 / 645120,
                                                         2179.0 / 645120,
                                                         60657.0 / 645120,
                                                         259723.0 / 645120,
                                                         259723.0 / 645120,
                                                         60657.0 / 645120,
                                                         2179.0 / 645120,
                                                         1.0 / 645120 }));
    // Independently evaluated Fortran definitions, with natural odd-width
    // coordinates on either side of the nearest-center change.
    error = std::max({ error,
                       sample_error(IBKernelEvaluatorIB5{},
                                    2.25,
                                    std::array<double, 5>{ 0.000539644595320609716,
                                                           0.128737522475479593,
                                                           0.514244366143986938,
                                                           0.333140121904304905,
                                                           0.0233383448809079538 }),
                       sample_error(IBKernelEvaluatorIB5{},
                                    1.75,
                                    std::array<double, 5>{ 0.0233383448809079538,
                                                           0.333140121904304905,
                                                           0.514244366143986938,
                                                           0.128737522475479593,
                                                           0.000539644595320609716 }),
                       sample_error(IBKernelEvaluatorIB6{},
                                    2.25,
                                    std::array<double, 6>{ 0.00965617417165844278,
                                                           0.174648694040214713,
                                                           0.431221688477088836,
                                                           0.325168575099164853,
                                                           0.0591221373512527211,
                                                           0.000182730860620434541 }),
                       sample_error(IBKernelEvaluatorIB6{},
                                    2.75,
                                    std::array<double, 6>{ 0.000182730860620434541,
                                                           0.0591221373512527211,
                                                           0.325168575099164853,
                                                           0.431221688477088836,
                                                           0.174648694040214713,
                                                           0.00965617417165844278 }) });
    const double moments = std::max({ moment_error(IBKernelEvaluatorBSpline<2>{}),
                                      moment_error(IBKernelEvaluatorBSpline<3>{}),
                                      moment_error(IBKernelEvaluatorBSpline<4>{}),
                                      moment_error(IBKernelEvaluatorBSpline<5>{}),
                                      moment_error(IBKernelEvaluatorBSpline<6>{}),
                                      moment_error(IBKernelEvaluatorBSpline<7>{}),
                                      moment_error(IBKernelEvaluatorBSpline<8>{}),
                                      moment_error(IBKernelEvaluatorIB3{}),
                                      moment_error(IBKernelEvaluatorIB4{}),
                                      moment_error(IBKernelEvaluatorIB5{}),
                                      moment_error(IBKernelEvaluatorIB6{}) });

    // Exercise every component axis of a 3D tensor product in this 2D executable.
    // These are evaluator checks, not a 3D hierarchy or matrix test.
    const IBKernelTensorProductEvaluator product{ IBKernelEvaluatorIB4{}, IBKernelEvaluatorIB3{} };
    const IBKernelTensorProductEvaluator isotropic3{ IBKernelEvaluatorBSpline<3>{} };
    const IBKernelTensorProductEvaluator isotropic5{ IBKernelEvaluatorBSpline<5>{} };
    const auto weights27 = isotropic3.evaluate<0>(std::array<double, 3>{ 1.0, 1.0, 1.0 });
    const auto weights125 = isotropic5.evaluate<2>(std::array<double, 3>{ 1.5, 1.5, 1.5 });
    static_assert(weights27.size() == 27 && weights125.size() == 125, "Natural 3D stencil sizes");
    const auto w0 = product.evaluate<0>(std::array<double, 3>{ 1.5, 1.0, 1.0 });
    const auto w1 = product.evaluate<1>(std::array<double, 3>{ 1.0, 1.5, 1.0 });
    const auto w2 = product.evaluate<2>(std::array<double, 3>{ 1.0, 1.0, 1.5 });
    const auto factors = product.evaluateFactors<1>(std::array<double, 3>{ 1.0, 1.5, 1.0 });
    static_assert(std::tuple_size<std::remove_reference_t<decltype(std::get<0>(factors))>>::value == 3,
                  "Natural factor width");
    static_assert(w0.size() == 36 && w1.size() == 36 && w2.size() == 36, "No tensor-product padding");
    const std::array<double, 4> normal = { a, b, b, a };
    const std::array<double, 3> tangent = { 1.0 / 6, 2.0 / 3, 1.0 / 6 };
    double tensor_error = 0.0;
    tensor_error =
        std::max(std::abs(weights27[13] - 0.75 * 0.75 * 0.75), std::abs(weights125[0] - 1.0 / (24.0 * 24.0 * 24.0)));
    for (int i = 0; i < 3; ++i)
        tensor_error = std::max({ tensor_error,
                                  std::abs(std::get<0>(factors)[i] - tangent[i]),
                                  std::abs(std::get<2>(factors)[i] - tangent[i]) });
    for (int i = 0; i < 4; ++i) tensor_error = std::max(tensor_error, std::abs(std::get<1>(factors)[i] - normal[i]));
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
            {
                const double expected = normal[i] * tangent[j] * tangent[k];
                tensor_error = std::max({ tensor_error,
                                          std::abs(w0[i + 4 * (j + 3 * k)] - expected),
                                          std::abs(w1[j + 3 * (i + 4 * k)] - expected),
                                          std::abs(w2[j + 3 * (k + 3 * i)] - expected) });
            }
    pout << "kernel_sample_max_error = " << error << '\n';
    pout << "kernel_moment_max_error = " << moments << '\n';
    pout << "tensor_product_3d_max_error = " << tensor_error << '\n';
    return error > 1.0e-12 || moments > 1.0e-12 || tensor_error > 1.0e-12;
}

void
register_probe_kernels()
{
    const IBKernel probe_kernel("PROBE");
    const auto register_pair = [&](const IBKernel& other, auto evaluator)
    {
        using E = decltype(evaluator);
        IBOperatorRegistry::register_interpolation_matrix_sc({ probe_kernel, other },
                                                             IBKernelTensorProductEvaluator{ ProbeEvaluator{}, E{} });
        IBOperatorRegistry::register_interpolation_matrix_sc({ other, probe_kernel },
                                                             IBKernelTensorProductEvaluator{ E{}, ProbeEvaluator{} });
    };
    register_pair(IBKernel::BSPLINE_1, IBKernelEvaluatorBSpline<1>{});
    register_pair(IBKernel::BSPLINE_2, IBKernelEvaluatorBSpline<2>{});
    register_pair(IBKernel::IB_4, IBKernelEvaluatorIB4{});
    IBOperatorRegistry::register_interpolation_matrix_sc(
        probe_kernel, IBKernelTensorProductEvaluator{ ProbeEvaluator{}, ProbeEvaluator{} });
}

void
generate_probes(const unsigned int& structure,
                const int& level,
                int& count,
                std::vector<IBTK::Point>& positions,
                void* ctx)
{
    if (ctx && *static_cast<bool*>(ctx))
    {
        count = 1;
        positions.resize(1);
        for (int d = 0; d < NDIM; ++d) positions[0](d) = 0.25 / 16.0;
        return;
    }
    TBOX_ASSERT(structure == 0 && level == 0);
    count = probe.size() * probe.size();
    positions.resize(count);
    for (unsigned int j = 0; j < probe.size(); ++j)
        for (unsigned int i = 0; i < probe.size(); ++i)
        {
            positions[j * probe.size() + i](0) = probe[i] / 16.0;
            positions[j * probe.size() + i](1) = probe[j] / 16.0;
        }
}

double
expected_weight(int width, int offset, double distance, bool bspline = false)
{
    if (bspline && width > 2)
    {
        // Independent truncated-power definition of the centered cardinal
        // B-spline, rather than the production knot-interval recurrence.
        double result = 0.0, binomial = 1.0, factorial = 1.0;
        for (int i = 2; i < width; ++i) factorial *= i;
        for (int k = 0; k <= width; ++k)
        {
            const double x = std::max(0.0, distance + 0.5 * width - k);
            result += (k % 2 ? -1.0 : 1.0) * binomial * std::pow(x, width - 1) / factorial;
            binomial *= static_cast<double>(width - k) / (k + 1);
        }
        return result;
    }
    if (width == 1) return 1.0;
    if (width == 2) return std::max(0.0, 1.0 - std::abs(distance));
    if (width == 3)
    {
        const double r_lower = distance + offset;
        return offset == 0 ? 0.125 + r_lower / 32.0 : (offset == 1 ? 0.375 - r_lower / 32.0 : 0.5);
    }
    // The symmetric radial definition of the original even-width IB_4 kernel.
    const double r = std::abs(distance);
    if (r <= 1.0) return (3.0 - 2.0 * r + std::sqrt(1.0 + 4.0 * r - 4.0 * r * r)) / 8.0;
    return (5.0 - 2.0 * r - std::sqrt(-7.0 + 12.0 * r - 4.0 * r * r)) / 8.0;
}

bool
check_matrix(Mat matrix,
             Vec positions,
             Pointer<SideData<NDIM, int>> dofs,
             int component_width,
             int transverse_width,
             bool bspline = false)
{
    PetscErrorCode ierr;
    PetscInt begin, end;
    ierr = MatGetOwnershipRange(matrix, &begin, &end);
    IBTK_CHKERRQ(ierr);
    if (begin != 0 || end != static_cast<PetscInt>(NDIM * probe.size() * probe.size())) return false;
    const PetscScalar* coordinates;
    ierr = VecGetArrayRead(positions, &coordinates);
    IBTK_CHKERRQ(ierr);
    bool valid = true;
    for (PetscInt row = begin; row < end; ++row)
    {
        const int axis = row % NDIM;
        std::array<int, NDIM> lower, width, sample;
        for (int d = 0; d < NDIM; ++d)
        {
            const double q = 16.0 * PetscRealPart(coordinates[(row - begin) / NDIM * NDIM + d]);
            sample[d] = 0;
            while (sample[d] < static_cast<int>(probe.size()) && q != probe[sample[d]]) ++sample[d];
            if (sample[d] == static_cast<int>(probe.size())) TBOX_ERROR("Unexpected interpolation probe\n");
            width[d] = d == axis ? component_width : transverse_width;
            if (width[d] % 2)
                lower[d] = (d == axis ? side_nearest : cell_nearest)[sample[d]] - width[d] / 2;
            else
                lower[d] = (d == axis ? side_even_lower : cell_even_lower)[sample[d]] - (width[d] / 2 - 1);
        }
        std::map<PetscInt, double> expected;
        for (int j = 0; j < width[1]; ++j)
            for (int i = 0; i < width[0]; ++i)
            {
                SAMRAI::hier::Index<NDIM> index;
                index(0) = lower[0] + i;
                index(1) = lower[1] + j;
                const SideIndex<NDIM> side(index, axis, SideIndex<NDIM>::Lower);
                const double x = probe[sample[0]] - (index(0) + (axis == 0 ? 0.0 : 0.5));
                const double y = probe[sample[1]] - (index(1) + (axis == 1 ? 0.0 : 0.5));
                const int column = (*dofs)(side);
                valid = valid && column >= 0;
                expected[column] = expected_weight(width[0], i, x, bspline) * expected_weight(width[1], j, y, bspline);
            }
        PetscInt count;
        const PetscInt* columns;
        const PetscScalar* values;
        ierr = MatGetRow(matrix, row, &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
        valid = valid && count == static_cast<PetscInt>(expected.size());
        for (PetscInt k = 0; k < count; ++k)
        {
            const auto found = expected.find(columns[k]);
            valid = valid && found != expected.end() && std::isfinite(PetscRealPart(values[k])) &&
                    std::abs(PetscRealPart(values[k]) - found->second) < 1.0e-12;
        }
        ierr = MatRestoreRow(matrix, row, &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecRestoreArrayRead(positions, &coordinates);
    IBTK_CHKERRQ(ierr);
    return valid;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    // Keep optional visualization warnings out of the compared test output.
    Logger::getInstance()->setWarning(false);
    const std::string input_file = argc > 1 ? argv[1] : "";
    if (input_file.find("registration.") != std::string::npos)
    {
        Pointer<Logger::Appender> appender = new TestAppender();
        Logger::getInstance()->setAbortAppender(appender);
        PIO::logOnlyNodeZero("output");
        if (input_file.find("duplicate") != std::string::npos ||
            input_file.find("registration.unknown") != std::string::npos)
        {
            std::ifstream input(input_file);
            std::string kernel_name, transverse_name;
            input >> kernel_name;
            const IBKernelTensorProduct kernel =
                input >> transverse_name ? IBKernelTensorProduct{ IBKernel(kernel_name), IBKernel(transverse_name) } :
                                           IBKernelTensorProduct{ IBKernel(kernel_name) };
            IBOperatorRegistry::register_interpolation_matrix_sc(
                kernel, IBKernelTensorProductEvaluator{ IBKernelEvaluatorIB4{} });
            return 0;
        }
    }
    PetscErrorCode ierr;
    int failures = 0;
    {
        Pointer<AppInitializer> app = new AppInitializer(argc, argv, "interpolation.log");
        const std::string setup_test = app->getInputDatabase()->getStringWithDefault("setup_test", "");
        bool setup_probe = !setup_test.empty();
        const bool late_fixed = setup_test == "late_fixed";
        if (setup_probe)
        {
            Pointer<Logger::Appender> appender = new TestAppender();
            Logger::getInstance()->setAbortAppender(appender);
        }
        constexpr int max_bspline_order = IBTK_MAX_BSPLINE_ORDER;
        const bool unsupported = input_file.find("registration.unsupported") != std::string::npos;
        if (!unsupported && !setup_probe) failures += check_kernels();
        Pointer<IBMethod> method = new IBMethod("IBMethod", app->getComponentDatabase("IBMethod"));
        method->setUseFixedLEOperators(!late_fixed);
        Pointer<IBStandardForceGen> force = new IBStandardForceGen();
        method->registerIBLagrangianForceFunction(force);
        Pointer<CartesianGridGeometry<NDIM>> geometry =
            new CartesianGridGeometry<NDIM>("CartesianGeometry", app->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", geometry);
        Pointer<StandardTagAndInitialize<NDIM>> tagger = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", method, app->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
            "GriddingAlgorithm", app->getComponentDatabase("GriddingAlgorithm"), tagger, boxes, balancer);
        Pointer<IBRedundantInitializer> initializer =
            new IBRedundantInitializer("IBRedundantInitializer", app->getComponentDatabase("IBRedundantInitializer"));
        initializer->setStructureNamesOnLevel(0, { "probes" });
        initializer->registerInitStructureFunction(generate_probes, &setup_probe);
        method->registerLInitStrategy(initializer);
        gridding->makeCoarsestLevel(hierarchy, 0.0);
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        // Interior probes check centering; setup probes exercise periodic ghost support.
        if (IBTK_MPI::getNodes() != 1 || level->getNumberOfPatches() != 1)
            TBOX_ERROR("Interpolation fixture requires one patch on one rank\n");
        VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = variables->getContext("interpolation");
        Pointer<SideVariable<NDIM, double>> velocity = new SideVariable<NDIM, double>("velocity");
        Pointer<SideVariable<NDIM, int>> indices = new SideVariable<NDIM, int>("indices");
        const int u = variables->registerVariableAndContext(velocity, context, method->getMinimumGhostCellWidth());
        // Match the implicit integrator's DOF allocation, including IBMethod's input minimum.
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

        if (setup_probe)
        {
            method->preprocessIntegrateData(0.0, 0.125, 1);
            if (late_fixed) method->setUseFixedLEOperators(true);
            method->updateFixedLEOperators();
            if (setup_test == "wide_kernel")
            {
                // Applications may still register this evaluator with a smaller built-in catalog.
                if (max_bspline_order < 8)
                    IBOperatorRegistry::register_interpolation_matrix_sc(
                        IBKernel("BSPLINE_8"), IBKernelTensorProductEvaluator{ IBKernelEvaluatorBSpline<8>{} });
                Mat matrix = nullptr;
                method->constructInterpOp(matrix, IBKernel("BSPLINE_8"), counts, dof, 0.125);
                Vec ones = nullptr, residual = nullptr;
                ierr = MatCreateVecs(matrix, &ones, &residual);
                IBTK_CHKERRQ(ierr);
                ierr = VecSet(ones, 1.0);
                IBTK_CHKERRQ(ierr);
                ierr = MatMult(matrix, ones, residual);
                IBTK_CHKERRQ(ierr);
                ierr = VecShift(residual, -1.0);
                IBTK_CHKERRQ(ierr);
                PetscReal error = 0.0;
                ierr = VecNorm(residual, NORM_INFINITY, &error);
                IBTK_CHKERRQ(ierr);
                failures += !std::isfinite(error) || error > 1.0e-12;
                pout << "dof_ghost_width = " << method->getMinimumGhostCellWidth()(0) << '\n'
                     << "constant_interpolation_error = " << error << '\n';
                ierr = VecDestroy(&ones);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&residual);
                IBTK_CHKERRQ(ierr);
                ierr = MatDestroy(&matrix);
                IBTK_CHKERRQ(ierr);
            }
            method->postprocessIntegrateData(0.0, 0.125, 1);
            method->postprocessData();
            level->deallocatePatchData(u);
            level->deallocatePatchData(dof);
            variables->removePatchDataIndex(u);
            variables->removePatchDataIndex(dof);
            pout << "test_failures = " << failures << std::endl;
            return failures;
        }

        if (unsupported)
        {
            method->preprocessIntegrateData(0.0, 0.125, 1);
            method->updateFixedLEOperators();
            Mat matrix = nullptr;
            method->constructInterpOp(
                matrix, IBKernel(app->getInputDatabase()->getString("matrix_kernel")), counts, dof, 0.125);
            ierr = MatDestroy(&matrix);
            IBTK_CHKERRQ(ierr);
            method->postprocessIntegrateData(0.0, 0.125, 1);
            method->postprocessData();
        }

        const std::array<IBKernel, 4> kernel = {
            IBKernel::BSPLINE_1, IBKernel::BSPLINE_2, IBKernel("PROBE"), IBKernel::IB_4
        };
        bool placement = true, scalar_equivalence = true, lifecycle = true;
        // Repeat setup/use/cleanup on the same IBMethod to exercise scratch state
        // invalidation, with fixed-operator updates at both new and half times.
        for (int step = 0; step < 2; ++step)
        {
            const double current_time = step * 0.125, new_time = current_time + 0.125;
            method->preprocessIntegrateData(current_time, new_time, 1);
            method->updateFixedLEOperators();
            Vec X = method->getLDataManager()->getLData("X", 0)->getVec();
            // Direct construction precedes registration and uses distinct functions
            // with the same widths. Both paths must retain the actual weights.
            Mat direct = nullptr, builtin = nullptr, registered = nullptr;
            PETScMatUtilities::constructPatchLevelSCInterpOp(
                direct,
                IBKernelTensorProductEvaluator{ ProbeEvaluator{}, IBKernelEvaluatorBSpline<2>{} },
                X,
                counts,
                dof,
                level);
            PETScMatUtilities::constructPatchLevelSCInterpOp(
                builtin,
                IBKernelTensorProductEvaluator{ IBKernelEvaluatorBSpline<3>{}, IBKernelEvaluatorBSpline<2>{} },
                X,
                counts,
                dof,
                level);
            placement = check_matrix(direct, X, dofs, 3, 2) && placement;
            placement = check_matrix(builtin, X, dofs, 3, 2, true) && placement;
            PetscBool equal;
            ierr = MatEqual(direct, builtin, &equal);
            IBTK_CHKERRQ(ierr);
            if (equal) ++failures;
            if (step == 0) register_probe_kernels();
            method->constructInterpOp(registered, { IBKernel("PROBE"), IBKernel::BSPLINE_2 }, counts, dof, new_time);
            ierr = MatEqual(direct, registered, &equal);
            IBTK_CHKERRQ(ierr);
            if (!equal) ++failures;
            if (max_bspline_order >= 3)
            {
                method->constructInterpOp(
                    registered, { IBKernel::BSPLINE_3, IBKernel::BSPLINE_2 }, counts, dof, new_time);
                ierr = MatEqual(builtin, registered, &equal);
                IBTK_CHKERRQ(ierr);
                if (!equal) ++failures;
            }

            // Direct evaluation is independent of the configured registrations.
            // Supply this combination when it is outside the automatic catalog.
            PETScMatUtilities::constructPatchLevelSCInterpOp(
                direct,
                IBKernelTensorProductEvaluator{ IBKernelEvaluatorBSpline<7>{}, IBKernelEvaluatorBSpline<2>{} },
                X,
                counts,
                dof,
                level);
            placement = check_matrix(direct, X, dofs, 7, 2, true) && placement;
            if (step == 0 && max_bspline_order < 7)
                IBOperatorRegistry::register_interpolation_matrix_sc(
                    { IBKernel("BSPLINE_7"), IBKernel::BSPLINE_2 },
                    IBKernelTensorProductEvaluator{ IBKernelEvaluatorBSpline<7>{}, IBKernelEvaluatorBSpline<2>{} });
            method->constructInterpOp(
                registered, { IBKernel("BSPLINE_7"), IBKernel::BSPLINE_2 }, counts, dof, new_time);
            ierr = MatEqual(direct, registered, &equal);
            IBTK_CHKERRQ(ierr);
            if (!equal) ++failures;

            // The B-spline limit does not restrict the other supplied kernels.
            PETScMatUtilities::constructPatchLevelSCInterpOp(
                direct, IBKernelTensorProductEvaluator{ IBKernelEvaluatorIB6{} }, X, counts, dof, level);
            method->constructInterpOp(registered, IBKernel::IB_6, counts, dof, new_time);
            ierr = MatEqual(direct, registered, &equal);
            IBTK_CHKERRQ(ierr);
            if (!equal) ++failures;
            ierr = MatDestroy(&direct);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&builtin);
            IBTK_CHKERRQ(ierr);
            ierr = MatDestroy(&registered);
            IBTK_CHKERRQ(ierr);
            for (int cw = 1; cw <= 4; ++cw)
                for (int tw = 1; tw <= 4; ++tw)
                {
                    if (max_bspline_order < 2 && (cw == 2 || tw == 2)) continue;
                    Mat matrix = nullptr;
                    IBImplicitStrategy& strategy = *method;
                    strategy.constructInterpOp(matrix, { kernel[cw - 1], kernel[tw - 1] }, counts, dof, new_time);
                    placement = check_matrix(matrix, X, dofs, cw, tw) && placement;
                    if (cw == tw)
                    {
                        Mat scalar = nullptr;
                        strategy.constructInterpOp(scalar, kernel[cw - 1], counts, dof, new_time);
                        PetscBool equal;
                        ierr = MatEqual(matrix, scalar, &equal);
                        IBTK_CHKERRQ(ierr);
                        scalar_equivalence = scalar_equivalence && equal;
                        IBOperatorRegistry::construct_interpolation_matrix_sc(
                            scalar, kernel[cw - 1], X, counts, dof, level);
                        ierr = MatEqual(matrix, scalar, &equal);
                        IBTK_CHKERRQ(ierr);
                        scalar_equivalence = scalar_equivalence && equal;
                        ierr = MatDestroy(&scalar);
                        IBTK_CHKERRQ(ierr);
                    }
                    ierr = MatDestroy(&matrix);
                    IBTK_CHKERRQ(ierr);
                }
            // Keep the highest-order sample inside this fixed patch geometry.
            const int highest_order = std::min(max_bspline_order, 8);
            const int transverse_order = std::min(max_bspline_order, 2);
            for (bool normal : { true, false })
            {
                const int cw = normal ? highest_order : transverse_order;
                const int tw = normal ? transverse_order : highest_order;
                Mat matrix = nullptr;
                method->constructInterpOp(
                    matrix,
                    { IBKernel("BSPLINE_" + std::to_string(cw)), IBKernel("BSPLINE_" + std::to_string(tw)) },
                    counts,
                    dof,
                    new_time);
                placement = check_matrix(matrix, X, dofs, cw, tw, true) && placement;
                ierr = MatDestroy(&matrix);
                IBTK_CHKERRQ(ierr);
            }
            // Natural named odd widths in both orientations, including mixed
            // parity, checked against independent scalar B-spline values.
            for (int cw : { 2, 3, 5 })
                for (int tw : { 2, 3, 5 })
                {
                    if (cw > max_bspline_order || tw > max_bspline_order) continue;
                    Mat matrix = nullptr;
                    method->constructInterpOp(
                        matrix,
                        { IBKernel("BSPLINE_" + std::to_string(cw)), IBKernel("BSPLINE_" + std::to_string(tw)) },
                        counts,
                        dof,
                        new_time);
                    placement = check_matrix(matrix, X, dofs, cw, tw, true) && placement;
                    ierr = MatDestroy(&matrix);
                    IBTK_CHKERRQ(ierr);
                }
            Mat half = nullptr;
            method->constructInterpOp(half, kernel[0], counts, dof, current_time + 0.0625);
            lifecycle = check_matrix(half, X, dofs, 1, 1) && lifecycle;
            ierr = MatDestroy(&half);
            IBTK_CHKERRQ(ierr);
            method->postprocessIntegrateData(current_time, new_time, 1);
            method->postprocessData();
        }
        failures += !placement + !scalar_equivalence + !lifecycle;
        pout << "stencil_columns_and_weights_valid = " << (placement ? "true" : "false") << '\n';
        pout << "isotropic_equivalence_valid = " << (scalar_equivalence ? "true" : "false") << '\n';
        pout << "standalone_lifecycle_valid = " << (lifecycle ? "true" : "false") << '\n';
        level->deallocatePatchData(u);
        level->deallocatePatchData(dof);
        variables->removePatchDataIndex(u);
        variables->removePatchDataIndex(dof);
        pout << "test_failures = " << failures << std::endl;
    }
    return failures;
}
