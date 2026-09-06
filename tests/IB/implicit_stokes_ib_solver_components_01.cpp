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

// Direct interpolation and implicit-operator contracts with a live IBMethod.
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/StaggeredStokesIBJacobianFACPreconditioner.h>
#include <ibamr/StaggeredStokesIBJacobianOperator.h>
#include <ibamr/StaggeredStokesIBLevelRelaxationFACOperator.h>
#include <ibamr/StaggeredStokesIBOperator.h>
#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/ibamr_enums.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CCPoissonPETScLevelSolver.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/HierarchyMathOps.h>
#include <ibtk/IBKernelEvaluators.h>
#include <ibtk/IBKernelTensorProductEvaluator.h>
#include <ibtk/IBOperatorRegistry.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/PETScKrylovLinearSolver.h>
#include <ibtk/PETScMFFDJacobianOperator.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>
#include <ibtk/SCPoissonPETScLevelSolver.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <tbox/Logger.h>
#include <tbox/MemoryDatabase.h>

#if defined(PETSC_USE_LOG)
#include <petsclog.h>
#endif

#include <BergerRigoutsos.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <HierarchyCellDataOpsReal.h>
#include <HierarchySideDataOpsReal.h>
#include <LoadBalancer.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PatchHierarchy.h>
#include <ProcessorMapping.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SAMRAI_config.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
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
run_interpolation(Pointer<AppInitializer> app, const std::string& input_file)
{
    PetscErrorCode ierr;
    int failures = 0;
    {
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

namespace
{
using HierarchyVector = SAMRAIVectorReal<NDIM, double>;

class CountedIBOperator : public StaggeredStokesIBOperator
{
public:
    CountedIBOperator() : StaggeredStokesIBOperator("nonlinear")
    {
    }
    void apply(HierarchyVector& x, HierarchyVector& y) override
    {
        ++evaluations;
        StaggeredStokesIBOperator::apply(x, y);
    }
    int evaluations = 0;
};

// Record dispatch while retaining the actual Stokes boundary operations.
class BoundaryCheckedStokesOperator : public StaggeredStokesOperator
{
public:
    BoundaryCheckedStokesOperator() : StaggeredStokesOperator("operator_test::stokes")
    {
    }
    void modifyRhsForBcs(HierarchyVector& y) override
    {
        ++rhs_calls;
        StaggeredStokesOperator::modifyRhsForBcs(y);
    }
    void imposeSolBcs(HierarchyVector& x) override
    {
        ++sol_calls;
        StaggeredStokesOperator::imposeSolBcs(x);
    }
    int rhs_calls = 0, sol_calls = 0;
};

void
generate_operator_structure(const unsigned int&, const int&, int& n, std::vector<IBTK::Point>& X, void* ctx)
{
    n = 16;
    X.resize(n);
    for (int k = 0; k < n; ++k)
    {
        const bool physical_boundary = ctx && *static_cast<bool*>(ctx);
        X[k](0) = (physical_boundary ? 0.1 : 0.5) + (physical_boundary ? 0.06 : 0.16) * std::cos(2.0 * M_PI * k / n);
        X[k](1) = 0.5 + 0.12 * std::sin(2.0 * M_PI * k / n);
    }
}

void
generate_operator_springs(
    const unsigned int&,
    const int&,
    std::multimap<int, IBRedundantInitializer::Edge>& edges,
    std::map<IBRedundantInitializer::Edge, IBRedundantInitializer::SpringSpec, IBRedundantInitializer::EdgeComp>& specs,
    void*)
{
    for (int k = 0; k < 16; ++k)
    {
        IBRedundantInitializer::Edge edge = { k, (k + 1) % 16 };
        if (edge.first > edge.second) std::swap(edge.first, edge.second);
        edges.emplace(edge.first, edge);
        IBRedundantInitializer::SpringSpec spring;
        spring.force_fcn_idx = 0;
        // Nonzero rest length makes the force Jacobian depend on the base state.
        spring.parameters = { 8.0, 0.02 };
        specs.emplace(edge, spring);
    }
}

void
set_operator_velocity(int idx, Pointer<PatchLevel<NDIM>> level, double amplitude, double phase)
{
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> data = patch->getPatchData(idx);
        Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
        for (int axis = 0; axis < NDIM; ++axis)
            for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(data->getGhostBox(), axis)); b; b++)
            {
                const int other = 1 - axis;
                const double q = geometry->getXLower()[other] +
                                 (b()(other) - patch->getBox().lower()(other) + 0.5) * geometry->getDx()[other];
                (*data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower)) =
                    amplitude * (axis == 0 ? 1.0 : -0.7) * std::sin(2.0 * M_PI * q + phase);
            }
    }
}

int
run_operators(Pointer<AppInitializer> app)
{
    PetscErrorCode ierr;
    bool physical_boundary = app->getInputDatabase()->getBoolWithDefault("physical_boundary", false);
    const TimeSteppingType type =
        IBAMR::string_to_enum<TimeSteppingType>(app->getInputDatabase()->getString("time_stepping"));
    const double current = 0.25, dt = 0.125, next = current + dt;
    const bool midpoint = type == MIDPOINT_RULE;
    const double force_scale = type == TRAPEZOIDAL_RULE ? 0.5 : 1.0;
    const double position_scale = type == BACKWARD_EULER ? 1.0 : (midpoint ? 0.25 : 0.5);
    const double force_time = midpoint ? current + dt / 2 : next;
    // Check the time/state choices through actual positions and independently
    // composed residual/Jacobian actions below, not the implementation's table.
    bool time_valid = true;
    Pointer<IBMethod> method = new IBMethod("IBMethod", app->getComponentDatabase("IBMethod"));
    method->setUseFixedLEOperators(true);
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
    initializer->setStructureNamesOnLevel(0, { "curve" });
    initializer->registerInitStructureFunction(generate_operator_structure, &physical_boundary);
    initializer->registerInitSpringDataFunction(generate_operator_springs);
    method->registerLInitStrategy(initializer);
    gridding->makeCoarsestLevel(hierarchy, current);
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    if (IBTK_MPI::getNodes() != 1 || level->getNumberOfPatches() != 1)
        TBOX_ERROR("Operator fixture requires one patch on one rank\n");

    auto* variables = VariableDatabase<NDIM>::getDatabase();
    auto context = variables->getContext("operators");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    Pointer<SideVariable<NDIM, int>> u_dof_var = new SideVariable<NDIM, int>("u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_var = new CellVariable<NDIM, int>("p_dof");
    const auto ghosts = method->getMinimumGhostCellWidth();
    std::vector<int> allocated;
    auto register_data =
        [&](Pointer<SAMRAI::hier::Variable<NDIM>> variable, const std::string& name, IntVector<NDIM> width)
    {
        const int idx = variables->registerVariableAndContext(variable, variables->getContext(name), width);
        level->allocatePatchData(idx, current);
        allocated.push_back(idx);
        return idx;
    };
    const int u_current = register_data(u_var, "current", ghosts);
    const int scratch = register_data(u_var, "scratch", ghosts);
    const int f_scratch = register_data(u_var, "force", ghosts);
    const int u_dof = register_data(u_dof_var, "dofs", ghosts);
    const int p_dof = register_data(p_dof_var, "dofs", IntVector<NDIM>(0));
    set_operator_velocity(u_current, level, 0.03, 0.2);
    std::vector<Pointer<CoarsenSchedule<NDIM>>> synch(1);
    std::vector<Pointer<RefineSchedule<NDIM>>> fill(1), prolong(1);
    std::vector<Pointer<LocationIndexRobinBcCoefs<NDIM>>> physical_coefs(NDIM);
    std::vector<RobinBcCoefStrategy<NDIM>*> bc_coefs(NDIM, nullptr);
    Pointer<CartSideRobinPhysBdryOp> physical_bc;
    if (physical_boundary)
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            physical_coefs[axis] = new LocationIndexRobinBcCoefs<NDIM>("velocity_bc", nullptr);
            for (int face = 0; face < 2 * NDIM; ++face)
                physical_coefs[axis]->setBoundaryValue(face, axis == 1 && face < 2 ? 0.3 : 0.0);
            bc_coefs[axis] = physical_coefs[axis];
        }
        physical_bc = new CartSideRobinPhysBdryOp(u_current, bc_coefs, false);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            physical_bc->setPhysicalBoundaryConditions(*level->getPatch(p()), current, ghosts);
        method->beginDataRedistribution(hierarchy, gridding);
        method->endDataRedistribution(hierarchy, gridding);
    }
    method->initializePatchHierarchy(hierarchy, gridding, u_current, synch, fill, 0, current, true);
    method->freeLInitStrategy();
    initializer.setNull();
    method->preprocessIntegrateData(current, next, 1);
    method->updateFixedLEOperators();
    method->interpolateVelocity(u_current, synch, fill, current);
    if (physical_boundary)
    {
        physical_bc->setPatchDataIndex(scratch);
        RefineAlgorithm<NDIM> ghost_fill;
        ghost_fill.registerRefine(scratch, scratch, scratch, nullptr);
        fill[0] = ghost_fill.createSchedule(level, physical_bc.getPointer());
    }

    Pointer<HierarchySideDataOpsReal<NDIM, double>> side_ops =
        new HierarchySideDataOpsReal<NDIM, double>(hierarchy, 0, 0);
    Pointer<HierarchyCellDataOpsReal<NDIM, double>> cell_ops =
        new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, 0, 0);
    HierarchyMathOps math_ops("operator_test::math", hierarchy);
    const int u = variables->registerVariableAndContext(u_var, context, ghosts);
    const int p = variables->registerVariableAndContext(p_var, context, IntVector<NDIM>(1));
    Pointer<HierarchyVector> base = new HierarchyVector("base", hierarchy, 0, 0);
    base->addComponent(u_var, u, math_ops.getSideWeightPatchDescriptorIndex(), side_ops);
    base->addComponent(p_var, p, math_ops.getCellWeightPatchDescriptorIndex(), cell_ops);
    base->allocateVectorData();
    std::vector<Pointer<HierarchyVector>> vectors = { base };
    auto clone = [&](const std::string& name)
    {
        auto vector = base->cloneVector(name);
        vector->allocateVectorData();
        vector->setToScalar(0.0);
        vectors.push_back(vector);
        return vector;
    };
    auto direction = clone("direction"), residual = clone("residual"), expected = clone("expected"),
         action = clone("action"), finite_difference = clone("finite_difference"), plus = clone("plus"),
         minus = clone("minus"), work = clone("work"), difference = clone("difference"),
         first_action = clone("first_action");
    set_operator_velocity(direction->getComponentDescriptorIndex(0), level, 0.2, 0.8);
    cell_ops->setToScalar(direction->getComponentDescriptorIndex(1), -0.25);

    if (physical_boundary)
    {
        // Exercise the live strategy near the wall, without assembling an IB matrix.
        Pointer<StaggeredStokesOperator> stokes = new StaggeredStokesOperator("boundary_stokes");
        PoissonSpecifications coefs("boundary_coefs");
        coefs.setCConstant(2.0);
        coefs.setDConstant(-0.01);
        stokes->setVelocityPoissonSpecifications(coefs);
        stokes->setPhysicalBcCoefs(bc_coefs, nullptr);
        StaggeredStokesIBOperator::Context ctx;
        ctx.ib_implicit_ops = method;
        ctx.stokes_op = stokes;
        ctx.hier_velocity_data_ops = side_ops;
        ctx.u_phys_bdry_op = physical_bc;
        ctx.u_synch_scheds = synch;
        ctx.u_ghost_fill_scheds = fill;
        ctx.f_prolongation_scheds = prolong;
        ctx.u_idx = scratch;
        ctx.f_idx = f_scratch;
        ctx.u_current_idx = u_current;
        ctx.time_stepping_type = BACKWARD_EULER;
        CountedIBOperator nonlinear;
        StaggeredStokesIBJacobianOperator jacobian("boundary_jacobian");
        nonlinear.setOperatorContext(ctx);
        jacobian.setOperatorContext(ctx);
        for (GeneralOperator* op :
             { static_cast<GeneralOperator*>(&nonlinear), static_cast<GeneralOperator*>(&jacobian) })
        {
            op->setTimeInterval(current, next);
            op->setSolutionTime(next);
        }
        Vec physical_position = nullptr;
        Vec X0 = method->getLDataManager()->getLData("X", 0)->getVec();
        ierr = VecDuplicate(X0, &physical_position);
        IBTK_CHKERRQ(ierr);
        double fd_error = 0.0, boundary_position_change = 0.0, coupling_norm = 0.0;
        int failures = 0;
        for (int cycle = 0; cycle < 2; ++cycle)
        {
            method->setUpdatedPosition(X0);
            nonlinear.initializeOperatorState(*base, *residual);
            jacobian.initializeOperatorState(*base, *residual);
            for (int state = 0; state < 2; ++state)
            {
                base->setToScalar(0.0);
                set_operator_velocity(u, level, state == 0 ? 0.12 : 0.24, state == 0 ? 0.1 : 0.6);
                nonlinear.apply(*base, *residual);
                std::vector<Pointer<LData>>* positions = nullptr;
                bool* needs_fill = nullptr;
                method->getPositionData(&positions, &needs_fill, TimePoint::NEW_TIME);
                ierr = VecCopy((*positions)[0]->getVec(), physical_position);
                IBTK_CHKERRQ(ierr);
                // A zero-boundary control proves the marker interpolation sees the wall data.
                for (int face = 0; face < 2; ++face) physical_coefs[1]->setBoundaryValue(face, 0.0);
                nonlinear.apply(*base, *expected);
                ierr = VecAXPY(physical_position, -1.0, (*positions)[0]->getVec());
                IBTK_CHKERRQ(ierr);
                PetscReal position_change = 0.0;
                ierr = VecNorm(physical_position, NORM_INFINITY, &position_change);
                IBTK_CHKERRQ(ierr);
                boundary_position_change = std::max(boundary_position_change, position_change);
                for (int face = 0; face < 2; ++face) physical_coefs[1]->setBoundaryValue(face, 0.3);
                nonlinear.apply(*base, *residual);
                jacobian.formJacobian(*base);
                jacobian.apply(*direction, *action);
                stokes->apply(*direction, *expected);
                difference->subtract(action, expected);
                const double ib_norm = difference->maxNorm();
                coupling_norm = std::max(coupling_norm, ib_norm);
                const double h = 1.0e-5;
                work->linearSum(1.0, base, h, direction);
                nonlinear.apply(*work, *plus);
                work->linearSum(1.0, base, -h, direction);
                nonlinear.apply(*work, *minus);
                finite_difference->linearSum(0.5 / h, plus, -0.5 / h, minus);
                difference->subtract(action, finite_difference);
                const double error = difference->maxNorm() / std::max(1.0, finite_difference->maxNorm());
                fd_error = std::max(fd_error, error);
                failures += !std::isfinite(error) || error > 1.0e-6 || !std::isfinite(position_change) ||
                            position_change <= 1.0e-4 || !std::isfinite(ib_norm) || ib_norm <= 1.0e-4;
            }
            jacobian.deallocateOperatorState();
            nonlinear.deallocateOperatorState();
        }
        pout << "physical_base_fd_error = " << std::scientific << std::setprecision(0) << fd_error << std::defaultfloat
             << std::setprecision(6) << '\n'
             << "boundary_position_change = " << boundary_position_change << '\n'
             << "nonzero_ib_action = " << coupling_norm << '\n';
        ierr = VecDestroy(&physical_position);
        IBTK_CHKERRQ(ierr);
        method->postprocessIntegrateData(current, next, 1);
        for (auto& vector : vectors) free_vector_components(*vector);
        for (int idx : allocated)
        {
            level->deallocatePatchData(idx);
            variables->removePatchDataIndex(idx);
        }
        pout << "test_failures = " << failures << std::endl;
        return failures;
    }

    std::vector<int> counts;
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(counts, u_dof, p_dof, level);
    Mat J = nullptr, A = nullptr;
    method->constructInterpOp(J, IBKernel::IB_4, counts, u_dof, force_time);
    method->constructLagrangianForceJacobian(A, MATAIJ, force_time);
    Vec eulerian = nullptr, spread = nullptr, interpolated = nullptr, position = nullptr, expected_position = nullptr;
    ierr = MatCreateVecs(J, &eulerian, &interpolated);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(eulerian, &spread);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(interpolated, &position);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(position, &expected_position);
    IBTK_CHKERRQ(ierr);
    Vec X0 = method->getLDataManager()->getLData("X", 0)->getVec();
    Pointer<LData> expected_force = new LData("expected_force", 16, NDIM);
    const double cell_volume = geometry->getDx()[0] * geometry->getDx()[1];
    Pointer<BoundaryCheckedStokesOperator> stokes = new BoundaryCheckedStokesOperator();
    PoissonSpecifications coefs("coefs");
    coefs.setCConstant(2.0);
    coefs.setDConstant(-0.01);
    stokes->setVelocityPoissonSpecifications(coefs);
    stokes->setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr), nullptr);
    StaggeredStokesIBOperator::Context ctx;
    ctx.ib_implicit_ops = method;
    ctx.stokes_op = stokes;
    ctx.hier_velocity_data_ops = side_ops;
    ctx.u_synch_scheds = synch;
    ctx.u_ghost_fill_scheds = fill;
    ctx.f_prolongation_scheds = prolong;
    ctx.patch_level = level;
    ctx.u_idx = scratch;
    ctx.f_idx = f_scratch;
    ctx.u_current_idx = u_current;
    ctx.u_dof_index_idx = u_dof;
    ctx.p_dof_index_idx = p_dof;
    ctx.time_stepping_type = type;
    CountedIBOperator nonlinear;
    StaggeredStokesIBJacobianOperator jacobian("jacobian");
    PETScMFFDJacobianOperator mffd("mffd");
    nonlinear.setOperatorContext(ctx);
    jacobian.setOperatorContext(ctx);
    mffd.setOperator(Pointer<GeneralOperator>(&nonlinear, false));
    for (GeneralOperator* op : { static_cast<GeneralOperator*>(&nonlinear),
                                 static_cast<GeneralOperator*>(&jacobian),
                                 static_cast<GeneralOperator*>(&mffd) })
    {
        op->setTimeInterval(current, next);
        op->setSolutionTime(force_time);
    }
    auto copy_to_petsc = [&](Vec v, Pointer<HierarchyVector> x)
    {
        StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            v, x->getComponentDescriptorIndex(0), u_dof, x->getComponentDescriptorIndex(1), p_dof, level);
    };
    auto close = [&](Pointer<HierarchyVector> lhs, Pointer<HierarchyVector> rhs, double tol)
    {
        difference->subtract(lhs, rhs);
        const double error = difference->maxNorm(), norm = rhs->maxNorm();
        return std::isfinite(error) && std::isfinite(norm) && error <= tol * std::max(1.0, norm);
    };
    bool residual_valid = true, derivative_valid = true, assembled_valid = true, base_valid = true,
         lifecycle_valid = true, boundary_valid = true, nontrivial = true;
    std::array<double, 5> errors = {};
    auto record_error = [&](int slot, Pointer<HierarchyVector> lhs, Pointer<HierarchyVector> rhs)
    {
        difference->subtract(lhs, rhs);
        errors[slot] = std::max(errors[slot], difference->maxNorm() / std::max(1.0, rhs->maxNorm()));
    };
    // Reinitialize all operator storage, and change the base twice per lifetime.
    for (int cycle = 0; cycle < 2; ++cycle)
    {
        ierr = VecCopy(X0, position);
        IBTK_CHKERRQ(ierr);
        // The second lifetime uses different fixed current/endpoint coupling
        // positions. Trapezoidal stepping must retain stored Lagrangian Ucurrent.
        if (cycle == 1)
        {
            ierr = VecShift(position, 0.01);
            IBTK_CHKERRQ(ierr);
        }
        method->setUpdatedPosition(position);
        nonlinear.initializeOperatorState(*base, *residual);
        jacobian.initializeOperatorState(*base, *residual);
        mffd.initializeOperatorState(*base, *residual);
        method->constructInterpOp(J, IBKernel::IB_4, counts, u_dof, force_time);
        for (int state = 0; state < 2; ++state)
        {
            base->setToScalar(0.0);
            set_operator_velocity(u, level, state == 0 ? 0.12 : 0.24, state == 0 ? 0.1 : 0.6);
            nonlinear.apply(*base, *residual);
            // Independent position/force composition using the live force law
            // and a directly assembled interpolation matrix, with no FAC helper.
            side_ops->linearSum(scratch, position_scale, u, midpoint ? position_scale : 0.0, u_current);
            side_ops->copyData(work->getComponentDescriptorIndex(0), scratch);
            cell_ops->setToScalar(work->getComponentDescriptorIndex(1), 0.0);
            copy_to_petsc(eulerian, work);
            ierr = MatMult(J, eulerian, interpolated);
            IBTK_CHKERRQ(ierr);
            ierr = VecWAXPY(position, dt, interpolated, X0);
            IBTK_CHKERRQ(ierr);
            if (type == TRAPEZOIDAL_RULE)
            {
                std::vector<Pointer<LData>>* U_current_data;
                method->getVelocityData(&U_current_data, TimePoint::CURRENT_TIME);
                ierr = VecAXPY(position, dt / 2, (*U_current_data)[0]->getVec());
                IBTK_CHKERRQ(ierr);
            }
            ierr = VecCopy(position, expected_position);
            IBTK_CHKERRQ(ierr);
            std::vector<Pointer<LData>>*X_data, *U_data;
            bool* X_ghost;
            const TimePoint force_point = midpoint ? TimePoint::HALF_TIME : TimePoint::NEW_TIME;
            method->getPositionData(&X_data, &X_ghost, force_point);
            method->getVelocityData(&U_data, force_point);
            ierr = VecAXPY(position, -1.0, (*X_data)[0]->getVec());
            IBTK_CHKERRQ(ierr);
            PetscReal position_error;
            ierr = VecNorm(position, NORM_INFINITY, &position_error);
            IBTK_CHKERRQ(ierr);
            time_valid = time_valid && std::isfinite(position_error) && position_error < 1.0e-12;
            stokes->apply(*base, *expected);
            ierr = VecSet(expected_force->getVec(), 0.0);
            IBTK_CHKERRQ(ierr);
            force->computeLagrangianForce(
                expected_force, (*X_data)[0], (*U_data)[0], hierarchy, 0, force_time, method->getLDataManager());
            ierr = MatMultTranspose(J, expected_force->getVec(), spread);
            IBTK_CHKERRQ(ierr);
            copy_to_petsc(eulerian, expected);
            ierr = VecAXPY(eulerian, -force_scale / cell_volume, spread);
            IBTK_CHKERRQ(ierr);
            StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(eulerian,
                                                                    expected->getComponentDescriptorIndex(0),
                                                                    u_dof,
                                                                    expected->getComponentDescriptorIndex(1),
                                                                    p_dof,
                                                                    level,
                                                                    nullptr,
                                                                    nullptr);
            residual_valid = close(residual, expected, 1.0e-11) && residual_valid;
            record_error(0, residual, expected);
            nontrivial = nontrivial && residual->maxNorm() > 1.0e-6;
            // Assemble -gamma*dt*alpha J^T A J at the actual force position.
            ierr = MatZeroEntries(A);
            IBTK_CHKERRQ(ierr);
            force->computeLagrangianForceJacobian(A,
                                                  MAT_FINAL_ASSEMBLY,
                                                  1.0,
                                                  (*X_data)[0],
                                                  0.0,
                                                  nullptr,
                                                  hierarchy,
                                                  0,
                                                  force_time,
                                                  method->getLDataManager());
            Mat coupling = nullptr;
            ierr = MatPtAP(A, J, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &coupling);
            IBTK_CHKERRQ(ierr);
            ierr = MatScale(coupling, -force_scale * dt * position_scale / cell_volume);
            IBTK_CHKERRQ(ierr);
            // Leave stale nonlinear state from a different input deliberately:
            // formJacobian must reconstruct its own supplied base.
            nonlinear.apply(*direction, *work);
            jacobian.formJacobian(*base);
            if (type != BACKWARD_EULER)
            {
                method->getPositionData(&X_data, &X_ghost, force_point);
                ierr = VecWAXPY(position, -1.0, expected_position, (*X_data)[0]->getVec());
                IBTK_CHKERRQ(ierr);
                ierr = VecNorm(position, NORM_INFINITY, &position_error);
                IBTK_CHKERRQ(ierr);
                time_valid = time_valid && std::isfinite(position_error) && position_error < 1.0e-12;
            }
            base_valid = close(jacobian.getBaseVector(), base, 0.0) && base_valid;
            jacobian.apply(*direction, *action);
            stokes->apply(*direction, *expected);
            nontrivial = nontrivial && !close(action, expected, 1.0e-5);
            if (state == 0)
                first_action->copyVector(action);
            else
                nontrivial = nontrivial && !close(action, first_action, 1.0e-7);
            jacobian.setIBCouplingJacobian(coupling);
            // Release the creator's reference, then reinstall the borrowed handle.
            // The operator must retain its own reference throughout replacement.
            Mat coupling_alias = coupling;
            ierr = MatDestroy(&coupling);
            IBTK_CHKERRQ(ierr);
            jacobian.setIBCouplingJacobian(coupling_alias);
            jacobian.apply(*direction, *expected);
            assembled_valid = close(action, expected, 1.0e-9) && assembled_valid;
            record_error(1, action, expected);
            Mat no_coupling = nullptr;
            jacobian.setIBCouplingJacobian(no_coupling);

            const double h = 1.0e-5;
            plus->linearSum(1.0, base, h, direction);
            minus->linearSum(1.0, base, -h, direction);
            nonlinear.apply(*plus, *expected);
            nonlinear.apply(*minus, *work);
            finite_difference->linearSum(0.5 / h, expected, -0.5 / h, work);
            derivative_valid = close(action, finite_difference, 1.0e-6) && derivative_valid;
            record_error(2, action, finite_difference);
            const int evaluations_before_form = nonlinear.evaluations;
            mffd.formJacobian(*base);
            base_valid = base_valid && nonlinear.evaluations == evaluations_before_form + 1;
            base_valid = close(mffd.getBaseVector(), base, 0.0) && base_valid;
            mffd.apply(*direction, *expected);
            derivative_valid = close(expected, finite_difference, 2.0e-5) && derivative_valid;
            record_error(3, expected, finite_difference);

            nonlinear.applyAdd(*base, *direction, *expected);
            work->add(residual, direction);
            residual_valid = close(expected, work, 1.0e-11) && residual_valid;
            record_error(4, expected, work);
            jacobian.applyAdd(*direction, *base, *expected);
            work->add(action, base);
            derivative_valid = close(expected, work, 1.0e-9) && derivative_valid;
        }
        for (GeneralOperator* op :
             { static_cast<GeneralOperator*>(&nonlinear), static_cast<GeneralOperator*>(&jacobian) })
            for (bool homogeneous : { false, true })
            {
                op->setHomogeneousBc(homogeneous);
                op->setSolutionTime(current + 0.03125);
                stokes->setTimeInterval(-2.0, -1.0);
                stokes->setSolutionTime(-1.0);
                stokes->setHomogeneousBc(!homogeneous);
                op->modifyRhsForBcs(*residual);
                boundary_valid = boundary_valid && stokes->getTimeInterval() == std::make_pair(current, next) &&
                                 stokes->getSolutionTime() == current + 0.03125 &&
                                 stokes->getHomogeneousBc() == homogeneous;
                stokes->setTimeInterval(-2.0, -1.0);
                stokes->setSolutionTime(-1.0);
                stokes->setHomogeneousBc(!homogeneous);
                op->imposeSolBcs(*base);
                boundary_valid = boundary_valid && stokes->getTimeInterval() == std::make_pair(current, next) &&
                                 stokes->getSolutionTime() == current + 0.03125 &&
                                 stokes->getHomogeneousBc() == homogeneous;
                op->setSolutionTime(force_time);
            }
        const std::array<int, 4> base_indices = { jacobian.getBaseVector()->getComponentDescriptorIndex(0),
                                                  jacobian.getBaseVector()->getComponentDescriptorIndex(1),
                                                  mffd.getBaseVector()->getComponentDescriptorIndex(0),
                                                  mffd.getBaseVector()->getComponentDescriptorIndex(1) };
        mffd.deallocateOperatorState();
        jacobian.deallocateOperatorState();
        nonlinear.deallocateOperatorState();
        lifecycle_valid = lifecycle_valid && !mffd.getIsInitialized() && !jacobian.getIsInitialized() &&
                          !nonlinear.getIsInitialized() && !jacobian.getBaseVector() && !mffd.getBaseVector();
        for (const int idx : base_indices)
        {
            Pointer<SAMRAI::hier::Variable<NDIM>> variable;
            lifecycle_valid = !variables->mapIndexToVariable(idx, variable) && lifecycle_valid;
        }
        for (const int idx : { u, p, u_current, scratch, f_scratch })
        {
            Pointer<SAMRAI::hier::Variable<NDIM>> variable;
            lifecycle_valid = variables->mapIndexToVariable(idx, variable) && level->getPatch(0)->getPatchData(idx) &&
                              lifecycle_valid;
        }
    }
    // A zero supplied contribution distinguishes its action from the nonzero
    // strategy Jacobian without inspecting which implementation is selected.
    Mat zero_coupling = nullptr;
    ierr = MatPtAP(A, J, MAT_INITIAL_MATRIX, 1.0, &zero_coupling);
    IBTK_CHKERRQ(ierr);
    ierr = MatZeroEntries(zero_coupling);
    IBTK_CHKERRQ(ierr);
    jacobian.initializeOperatorState(*base, *residual);
    jacobian.setIBCouplingJacobian(zero_coupling);
    jacobian.apply(*direction, *first_action);
    stokes->apply(*direction, *work);
    assembled_valid = close(first_action, work, 1.0e-12) && assembled_valid;

    PETScKrylovLinearSolver outer("initialization_sequence", nullptr, "initialization_sequence_");
    outer.setOperator(Pointer<LinearOperator>(&jacobian, false));
    outer.setTimeInterval(current, next);
    outer.setSolutionTime(force_time);
    outer.initializeSolverState(*base, *residual);
    lifecycle_valid = !jacobian.getBaseVector() && lifecycle_valid;
    jacobian.formJacobian(*base);
    jacobian.apply(*direction, *action);
    Mat no_coupling = nullptr;
    jacobian.setIBCouplingJacobian(no_coupling);
    jacobian.apply(*direction, *expected);
    difference->subtract(action, expected);
    const double strategy_error = difference->maxNorm();
    difference->subtract(action, first_action);
    const double selection_difference = difference->maxNorm();
    assembled_valid = close(action, expected, 1.0e-12) && !close(action, first_action, 1.0e-5) && assembled_valid;

    // Install the supplied matrix after the outer solver has initialized its
    // operator. This must restore the zero-coupling (Stokes-only) action.
    jacobian.setIBCouplingJacobian(zero_coupling);
    jacobian.apply(*direction, *expected);
    difference->subtract(expected, first_action);
    const double supplied_error = difference->maxNorm();
    assembled_valid = close(expected, first_action, 1.0e-12) && assembled_valid;
    outer.deallocateSolverState();
    lifecycle_valid = !jacobian.getIsInitialized() && !jacobian.getBaseVector() && lifecycle_valid;
    ierr = MatDestroy(&zero_coupling);
    IBTK_CHKERRQ(ierr);
    pout << "outer_initialization_strategy_error = " << strategy_error << '\n'
         << "outer_initialization_action_change = " << selection_difference << '\n'
         << "post_initialization_supplied_error = " << supplied_error << '\n';
    boundary_valid = boundary_valid && stokes->rhs_calls == 8 && stokes->sol_calls == 8;
    // Roundoff-sensitive output is compact; the checks above retain full precision.
    pout << "comparison_errors (residual, assembled, centered_fd, mffd, apply_add) = " << std::scientific
         << std::setprecision(0) << errors[0] << ' ' << errors[1] << ' ' << errors[2] << ' ' << errors[3] << ' '
         << errors[4] << std::defaultfloat << std::setprecision(6) << '\n';
    int failures = !residual_valid + !derivative_valid + !assembled_valid;
    for (const auto& check :
         std::vector<std::pair<std::string, bool>>{ { "time_state_scaling_valid", time_valid },
                                                    { "nontrivial_coupling_valid", nontrivial },
                                                    { "updated_base_state_valid", base_valid },
                                                    { "boundary_forwarding_valid", boundary_valid },
                                                    { "operator_lifecycle_valid", lifecycle_valid } })
    {
        pout << check.first << " = " << (check.second ? "true" : "false") << '\n';
        failures += !check.second;
    }
    method->postprocessIntegrateData(current, next, 1);
    for (Vec* v : { &eulerian, &spread, &interpolated, &position, &expected_position })
    {
        ierr = VecDestroy(v);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatDestroy(&A);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&J);
    IBTK_CHKERRQ(ierr);
    for (auto& vector : vectors) free_vector_components(*vector);
    for (int idx : allocated)
    {
        level->deallocatePatchData(idx);
        variables->removePatchDataIndex(idx);
    }
    pout << "test_failures = " << failures << std::endl;
    return failures;
}
} // namespace

namespace
{
// Expose existing protected state for lifecycle assertions, without changing
// the production API or replacing any solver operation.
template <class Solver>
class LevelSolverProbe : public Solver
{
public:
    LevelSolverProbe(const std::string& name, Pointer<Database> db) : Solver(name, db, "")
    {
    }
    Mat matrixBeforeKSP() const
    {
        return d_matrix_before_ksp;
    }
    PetscInt referencesBeforeKSP() const
    {
        return d_references_before_ksp;
    }
    bool shellStorageEmpty() const
    {
        return this->d_sub_x.empty() && this->d_sub_y.empty() && this->d_sub_ksp.empty() &&
               this->d_restriction.empty() && this->d_prolongation.empty();
    }
    std::vector<Vec> retainShellVectors()
    {
        std::vector<Vec> result = this->d_sub_x;
        result.insert(result.end(), this->d_sub_y.begin(), this->d_sub_y.end());
        for (Vec v : result)
        {
            PetscErrorCode ierr = PetscObjectReference(reinterpret_cast<PetscObject>(v));
            IBTK_CHKERRQ(ierr);
        }
        return result;
    }

protected:
    void initializeSolverStateSpecialized(const HierarchyVector& x, const HierarchyVector& b) override
    {
        Solver::initializeSolverStateSpecialized(x, b);
        d_matrix_before_ksp = this->d_petsc_mat;
        PetscErrorCode ierr =
            PetscObjectGetReference(reinterpret_cast<PetscObject>(d_matrix_before_ksp), &d_references_before_ksp);
        IBTK_CHKERRQ(ierr);
    }

private:
    Mat d_matrix_before_ksp = nullptr;
    PetscInt d_references_before_ksp = 0;
};

struct LevelFixture
{
    Pointer<PatchHierarchy<NDIM>> hierarchy;
    Pointer<PatchLevel<NDIM>> level;
    Pointer<HierarchyVector> x, b;
    std::vector<int> indices;
    std::vector<PetscInt> velocity_ids;
    int full_size = 0;

    LevelFixture(Pointer<Database> geometry_db, int ln = 0, bool full = true)
    {
        if (IBTK_MPI::getNodes() != 1) TBOX_ERROR("Level fixture requires one rank\n");
        Pointer<CartesianGridGeometry<NDIM>> geometry = new CartesianGridGeometry<NDIM>("level_geometry", geometry_db);
        hierarchy = new PatchHierarchy<NDIM>("level_hierarchy", geometry);
        BoxArray<NDIM> boxes(1);
        boxes[0] = Box<NDIM>(SAMRAI::hier::Index<NDIM>(0), SAMRAI::hier::Index<NDIM>(15));
        ProcessorMapping mapping(1);
        mapping.setProcessorAssignment(0, 0);
        hierarchy->makeNewPatchLevel(0, IntVector<NDIM>(1), boxes, mapping);
        if (ln == 1)
        {
            boxes[0] = full ? Box<NDIM>(SAMRAI::hier::Index<NDIM>(0), SAMRAI::hier::Index<NDIM>(31)) :
                              Box<NDIM>(SAMRAI::hier::Index<NDIM>(8), SAMRAI::hier::Index<NDIM>(23));
            hierarchy->makeNewPatchLevel(1, IntVector<NDIM>(2), boxes, mapping);
        }
        level = hierarchy->getPatchLevel(ln);
        auto* db = VariableDatabase<NDIM>::getDatabase();
        auto context = db->getContext("level_fixture");
        Pointer<SideVariable<NDIM, double>> u = new SideVariable<NDIM, double>("level_u");
        Pointer<CellVariable<NDIM, double>> p = new CellVariable<NDIM, double>("level_p");
        if (db->checkVariableExists("level_u")) u = db->getVariable("level_u");
        if (db->checkVariableExists("level_p")) p = db->getVariable("level_p");
        const int ui = db->registerVariableAndContext(u, context, IntVector<NDIM>(1));
        const int pi = db->registerVariableAndContext(p, context, IntVector<NDIM>(1));
        x = new HierarchyVector("level_x", hierarchy, ln, ln);
        x->addComponent(u, ui, -1, new HierarchySideDataOpsReal<NDIM, double>(hierarchy, ln, ln));
        x->addComponent(p, pi, -1, new HierarchyCellDataOpsReal<NDIM, double>(hierarchy, ln, ln));
        x->allocateVectorData();
        x->setToScalar(0.0);
        b = x->cloneVector("level_b");
        b->allocateVectorData();
        b->setToScalar(0.0);
        Pointer<SideVariable<NDIM, int>> ud = new SideVariable<NDIM, int>("level_ud");
        Pointer<CellVariable<NDIM, int>> pd = new CellVariable<NDIM, int>("level_pd");
        if (db->checkVariableExists("level_ud")) ud = db->getVariable("level_ud");
        if (db->checkVariableExists("level_pd")) pd = db->getVariable("level_pd");
        const int udi = db->registerVariableAndContext(ud, context, IntVector<NDIM>(1));
        const int pdi = db->registerVariableAndContext(pd, context, IntVector<NDIM>(1));
        indices = { udi, pdi };
        for (int idx : indices) level->allocatePatchData(idx);
        std::vector<int> counts;
        StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(counts, udi, pdi, level);
        full_size = counts[0];
        std::set<int> velocity;
        Pointer<SideData<NDIM, int>> data = level->getPatch(0)->getPatchData(udi);
        for (int axis = 0; axis < NDIM; ++axis)
            for (Box<NDIM>::Iterator i(SideGeometry<NDIM>::toSideBox(level->getPatch(0)->getBox(), axis)); i; i++)
                velocity.insert((*data)(SideIndex<NDIM>(i(), axis, SideIndex<NDIM>::Lower)));
        velocity_ids.assign(velocity.begin(), velocity.end());
    }
    ~LevelFixture()
    {
        free_vector_components(*b);
        free_vector_components(*x);
        for (int idx : indices)
        {
            level->deallocatePatchData(idx);
            VariableDatabase<NDIM>::getDatabase()->removePatchDataIndex(idx);
        }
    }
};

Pointer<Database>
level_solver_database(const std::string& pc = "none", int overlap = 0)
{
    Pointer<Database> db = new MemoryDatabase("level_solver");
    db->putString("ksp_type", "gmres");
    db->putString("options_prefix", "level_");
    db->putString("pc_type", pc);
    db->putString("shell_pc_type", "additive");
    db->putBool("initial_guess_nonzero", false);
    db->putDouble("rel_residual_tol", 1.0e-12);
    db->putInteger("max_iterations", 100);
    const int size[NDIM] = { 4, 8 }, width[NDIM] = { overlap, overlap };
    db->putIntegerArray("subdomain_box_size", size, NDIM);
    db->putIntegerArray("subdomain_overlap_size", width, NDIM);
    return db;
}

Mat
level_test_matrix(PetscInt n, double shift)
{
    Mat matrix;
    PetscErrorCode ierr = MatCreateAIJ(PETSC_COMM_WORLD, n, n, n, n, 2, nullptr, 0, nullptr, &matrix);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = 0; i < n; ++i)
    {
        const PetscInt cols[2] = { i, (i + 3) % n };
        const PetscScalar vals[2] = { shift + 0.001 * i, 0.125 };
        ierr = MatSetValues(matrix, 1, &i, 2, cols, vals, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    return matrix;
}

// Maximum errors across the matrix replacements and solver lifetimes in one case.
double max_matrix_error = 0.0, max_level_solve_error = 0.0, max_mapping_error = 0.0;

bool
matrices_equal(Mat a, Mat b)
{
    Mat diff;
    PetscErrorCode ierr = MatDuplicate(a, MAT_COPY_VALUES, &diff);
    IBTK_CHKERRQ(ierr);
    ierr = MatAXPY(diff, -1.0, b, DIFFERENT_NONZERO_PATTERN);
    IBTK_CHKERRQ(ierr);
    PetscReal norm;
    ierr = MatNorm(diff, NORM_INFINITY, &norm);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&diff);
    IBTK_CHKERRQ(ierr);
    max_matrix_error = std::max(max_matrix_error, norm);
    return std::isfinite(norm) && norm < 1.0e-12;
}

// Exercise the installed solver KSP and its actual configured matrix/PC.
bool
check_level_solve(PETScLevelSolver& solver)
{
    Mat matrix;
    PetscErrorCode ierr = KSPGetOperators(solver.getPETScKSP(), &matrix, nullptr);
    IBTK_CHKERRQ(ierr);
    Vec exact, rhs, solution;
    ierr = MatCreateVecs(matrix, &exact, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(exact, &solution);
    IBTK_CHKERRQ(ierr);
    PetscInt n;
    ierr = VecGetSize(exact, &n);
    IBTK_CHKERRQ(ierr);
    PetscScalar* values;
    ierr = VecGetArray(exact, &values);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = 0; i < n; ++i) values[i] = std::sin(0.13 * i) + 0.5;
    ierr = VecRestoreArray(exact, &values);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(matrix, exact, rhs);
    IBTK_CHKERRQ(ierr);
    ierr = KSPSolve(solver.getPETScKSP(), rhs, solution);
    IBTK_CHKERRQ(ierr);
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason(solver.getPETScKSP(), &reason);
    IBTK_CHKERRQ(ierr);
    ierr = VecAXPY(solution, -1.0, exact);
    IBTK_CHKERRQ(ierr);
    PetscReal error;
    ierr = VecNorm(solution, NORM_INFINITY, &error);
    IBTK_CHKERRQ(ierr);
    for (Vec* v : { &exact, &rhs, &solution })
    {
        ierr = VecDestroy(v);
        IBTK_CHKERRQ(ierr);
    }
    max_level_solve_error = std::max(max_level_solve_error, error);
    return reason > 0 && std::isfinite(error) && error < 1.0e-9;
}

bool
check_stokes_vector_mapping(StaggeredStokesPETScLevelSolver& solver, LevelFixture& fixture)
{
    const int u = fixture.x->getComponentDescriptorIndex(0), p = fixture.x->getComponentDescriptorIndex(1);
    set_operator_velocity(u, fixture.level, 0.3, 0.4);
    Pointer<CellData<NDIM, double>> pressure = fixture.level->getPatch(0)->getPatchData(p);
    pressure->fill(0.25);
    Mat matrix;
    PetscErrorCode ierr = KSPGetOperators(solver.getPETScKSP(), &matrix, nullptr);
    IBTK_CHKERRQ(ierr);
    Vec exact, rhs, result;
    ierr = MatCreateVecs(matrix, &exact, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(exact, &result);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        exact, u, fixture.indices[0], p, fixture.indices[1], fixture.level);
    ierr = MatMult(matrix, exact, rhs);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(rhs,
                                                            fixture.b->getComponentDescriptorIndex(0),
                                                            fixture.indices[0],
                                                            fixture.b->getComponentDescriptorIndex(1),
                                                            fixture.indices[1],
                                                            fixture.level,
                                                            nullptr,
                                                            nullptr);
    fixture.x->setToScalar(0.0);
    const bool converged = solver.solveSystem(*fixture.x, *fixture.b);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        result, u, fixture.indices[0], p, fixture.indices[1], fixture.level);
    ierr = VecAXPY(result, -1.0, exact);
    IBTK_CHKERRQ(ierr);
    PetscReal error;
    ierr = VecNorm(result, NORM_INFINITY, &error);
    IBTK_CHKERRQ(ierr);
    for (Vec* v : { &exact, &rhs, &result })
    {
        ierr = VecDestroy(v);
        IBTK_CHKERRQ(ierr);
    }
    max_mapping_error = std::max(max_mapping_error, error);
    return converged && std::isfinite(error) && error < 1.0e-9;
}

PetscInt
matrix_references(Mat matrix)
{
    PetscInt references;
    PetscErrorCode ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(matrix), &references);
    IBTK_CHKERRQ(ierr);
    return references;
}

bool
check_matrix_reference_lifetime(LevelFixture& fixture)
{
    Mat a = level_test_matrix(fixture.full_size, 4.0), b = level_test_matrix(fixture.full_size, 5.0);
    bool valid = true;
    {
        StaggeredStokesPETScLevelSolver solver("uninitialized_lifetime", level_solver_database(), "");
        // Each setter owns a reference, even when both inputs are the same Mat.
        solver.setOperatorMat(a);
        solver.setAugmentedOperatorMat(a);
        valid = matrix_references(a) == 3 && valid;
        solver.setOperatorMat(a);
        solver.setAugmentedOperatorMat(a);
        valid = matrix_references(a) == 3 && valid;
        solver.setOperatorMat(b);
        valid = matrix_references(a) == 2 && matrix_references(b) == 2 && valid;
        solver.setAugmentedOperatorMat(b);
        valid = matrix_references(a) == 1 && matrix_references(b) == 3 && valid;
        solver.setOperatorMat(nullptr);
        solver.setOperatorMat(nullptr);
        valid = matrix_references(b) == 2 && valid;
        solver.setOperatorMat(a);
        valid = matrix_references(a) == 2 && valid;
        // Destruction must release inputs even without initialization.
    }
    valid = matrix_references(a) == 1 && matrix_references(b) == 1 && valid;
    {
        StaggeredStokesPETScLevelSolver solver("initialized_lifetime", level_solver_database(), "");
        solver.setTimeInterval(0.0, 1.0);
        solver.setSolutionTime(1.0);
        solver.setOperatorMat(a);
        solver.setAugmentedOperatorMat(b);
        solver.initializeSolverState(*fixture.x, *fixture.b);
        valid = check_level_solve(solver) && valid;
        // Destruction also tears down KSP state before releasing both inputs.
    }
    valid = matrix_references(a) == 1 && matrix_references(b) == 1 && valid;
    PetscErrorCode ierr = MatDestroy(&a);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&b);
    IBTK_CHKERRQ(ierr);
    return valid;
}

int
run_level_operator(Pointer<AppInitializer> app, bool augmentation)
{
    LevelFixture fixture(app->getComponentDatabase("CartesianGeometry"));
    LevelSolverProbe<StaggeredStokesPETScLevelSolver> solver("supplied_stokes", level_solver_database());
    solver.setTimeInterval(0.0, 1.0);
    solver.setSolutionTime(1.0);
    bool identity = true, creator_valid = true, values_valid = true, solves = true;
    bool references_valid = check_matrix_reference_lifetime(fixture);
    for (int cycle = 0; cycle < (augmentation ? 4 : 2); ++cycle)
    {
        const bool full_augmentation = cycle % 2 == 0;
        Mat creator = level_test_matrix(fixture.full_size, 4.0 + cycle), original;
        PetscErrorCode ierr = MatDuplicate(creator, MAT_COPY_VALUES, &original);
        IBTK_CHKERRQ(ierr);
        solver.setOperatorMat(creator);
        references_valid = matrix_references(creator) == 2 && references_valid;
        Mat augmented = nullptr, expected = nullptr, augmented_original = nullptr;
        if (augmentation)
        {
            const PetscInt n = full_augmentation ? fixture.full_size : fixture.velocity_ids.size();
            augmented = level_test_matrix(n, 1.0);
            ierr = MatDuplicate(augmented, MAT_COPY_VALUES, &augmented_original);
            IBTK_CHKERRQ(ierr);
            ierr = MatDuplicate(creator, MAT_COPY_VALUES, &expected);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetOption(expected, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
            IBTK_CHKERRQ(ierr);
            // Independently embed known entries, using actual coupled DOF data.
            for (PetscInt row = 0; row < n; ++row)
            {
                const PetscInt next = (row + 3) % n;
                const PetscInt full_row = full_augmentation ? row : fixture.velocity_ids[row];
                const PetscInt cols[2] = { full_row, full_augmentation ? next : fixture.velocity_ids[next] };
                const PetscScalar vals[2] = { 1.0 + 0.001 * row, 0.125 };
                ierr = MatSetValues(expected, 1, &full_row, 2, cols, vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            ierr = MatAssemblyBegin(expected, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(expected, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
        }
        solver.setAugmentedOperatorMat(augmented);
        if (augmentation) references_valid = matrix_references(augmented) == 2 && references_valid;
        // Reject missing retention before testing caller release, so that a
        // regression reports failure instead of dereferencing a dangling Mat.
        if (!references_valid) TBOX_ERROR("Installed matrix references were not retained correctly.\n");
        const Mat operator_alias = creator, augmentation_alias = augmented;
        ierr = MatDestroy(&creator);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&augmented);
        IBTK_CHKERRQ(ierr);
        for (int lifetime = 0; lifetime < 2; ++lifetime)
        {
            // No caller-owned reference remains, even before initialization.
            references_valid = matrix_references(operator_alias) == 1 && references_valid;
            if (augmentation) references_valid = matrix_references(augmentation_alias) == 1 && references_valid;
            solver.initializeSolverState(*fixture.x, *fixture.b);
            Mat installed;
            ierr = KSPGetOperators(solver.getPETScKSP(), &installed, nullptr);
            IBTK_CHKERRQ(ierr);
            identity = identity && installed == solver.matrixBeforeKSP() &&
                       (augmentation ? installed != operator_alias : installed == operator_alias) &&
                       solver.referencesBeforeKSP() == 1;
            values_valid = matrices_equal(installed, augmentation ? expected : original) && values_valid;
            creator_valid = matrices_equal(operator_alias, original) && creator_valid;
            if (augmentation) creator_valid = matrices_equal(augmentation_alias, augmented_original) && creator_valid;
            solves = check_level_solve(solver) && solves;
            solves = check_stokes_vector_mapping(solver, fixture) && solves;
            solver.deallocateSolverState();
            // Same-handle replacement must also work with only the solver's
            // reference remaining. Aliases are used read-only throughout.
            solver.setOperatorMat(operator_alias);
            solver.setAugmentedOperatorMat(augmentation_alias);
            references_valid = matrix_references(operator_alias) == 1 && references_valid;
            if (augmentation) references_valid = matrix_references(augmentation_alias) == 1 && references_valid;
        }
        // Retain observer references to check that clearing releases exactly
        // the solver's references, without destroying another owner's data.
        creator = operator_alias;
        augmented = augmentation_alias;
        ierr = PetscObjectReference(reinterpret_cast<PetscObject>(creator));
        IBTK_CHKERRQ(ierr);
        ierr = PetscObjectReference(reinterpret_cast<PetscObject>(augmented));
        IBTK_CHKERRQ(ierr);
        solver.setOperatorMat(nullptr);
        solver.setAugmentedOperatorMat(nullptr);
        references_valid = matrix_references(creator) == 1 && references_valid;
        creator_valid = matrices_equal(creator, original) && creator_valid;
        if (augmentation)
        {
            references_valid = matrix_references(augmented) == 1 && references_valid;
            creator_valid = matrices_equal(augmented, augmented_original) && creator_valid;
        }
        for (Mat* m : { &creator, &original, &augmented, &expected, &augmented_original })
        {
            ierr = MatDestroy(m);
            IBTK_CHKERRQ(ierr);
        }
    }
    int failures = !identity + !creator_valid + !values_valid + !solves + !references_valid;
    pout << "matrix_identity_valid = " << (identity ? "true" : "false") << '\n'
         << "creator_lifetime_valid = " << (creator_valid ? "true" : "false") << '\n'
         << "retained_references_valid = " << (references_valid ? "true" : "false") << '\n'
         << "matrix_error = " << max_matrix_error << '\n'
         << "level_solve_error = " << max_level_solve_error << '\n'
         << "velocity_pressure_mapping_error = " << max_mapping_error << '\n'
         << "test_failures = " << failures << std::endl;
    return failures;
}

int
run_initialized_matrix_setter(Pointer<AppInitializer> app, bool augmentation)
{
    // Use the repository's path-independent expected-error output.
    Pointer<Logger::Appender> abort_appender = new TestAppender();
    Logger::getInstance()->setAbortAppender(abort_appender);
    LevelFixture fixture(app->getComponentDatabase("CartesianGeometry"));
    StaggeredStokesPETScLevelSolver solver("initialized_setter", level_solver_database(), "");
    solver.setTimeInterval(0.0, 1.0);
    solver.setSolutionTime(1.0);
    Mat matrix = level_test_matrix(fixture.full_size, 4.0);
    solver.setOperatorMat(matrix);
    if (augmentation) solver.setAugmentedOperatorMat(matrix);
    solver.initializeSolverState(*fixture.x, *fixture.b);
    if (augmentation)
        solver.setAugmentedOperatorMat(matrix);
    else
        solver.setOperatorMat(matrix);
    // Even same-handle setters must reject initialized state. Returning zero
    // on unexpected continuation makes this expect_error=true case fail attest.
    solver.deallocateSolverState();
    PetscErrorCode ierr = MatDestroy(&matrix);
    IBTK_CHKERRQ(ierr);
    pout << "ERROR: initialized matrix setter unexpectedly returned.\n";
    return 0;
}

template <class Solver>
bool
check_shell_state(LevelSolverProbe<Solver>& solver,
                  Pointer<HierarchyVector> x,
                  Pointer<HierarchyVector> b,
                  PetscInt& overlap_total)
{
    bool valid = true;
    for (int cycle = 0; cycle < 2; ++cycle)
    {
        solver.initializeSolverState(*x, *b);
        std::vector<IS>*nonoverlap, *overlap;
        solver.getASMSubdomains(&nonoverlap, &overlap);
        valid = valid && overlap->size() == 8 && nonoverlap->size() == 8;
        overlap_total = 0;
        for (IS is : *overlap)
        {
            PetscInt n;
            PetscErrorCode ierr = ISGetSize(is, &n);
            IBTK_CHKERRQ(ierr);
            overlap_total += n;
        }
        valid = check_level_solve(solver) && valid;
        auto retained = solver.retainShellVectors();
        valid = valid && retained.size() == 16;
        solver.deallocateSolverState();
        valid = solver.shellStorageEmpty() && valid;
        PetscInt max_references = 0;
        for (Vec& v : retained)
        {
            PetscInt references;
            PetscErrorCode ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(v), &references);
            IBTK_CHKERRQ(ierr);
            // Only this test reference may remain after solver teardown.
            valid = references == 1 && valid;
            max_references = std::max(max_references, references);
            ierr = VecDestroy(&v);
            IBTK_CHKERRQ(ierr);
        }
        if (max_references != 1) pout << "unreleased_shell_vector_references = " << max_references << std::endl;
    }
    return valid;
}

int
run_level_state(Pointer<AppInitializer> app)
{
    bool cc_valid = true, sc_valid = true, stokes_valid = true, domain_valid = true;
    {
        LevelFixture fixture(app->getComponentDatabase("CartesianGeometry"));
        Pointer<HierarchyVector> cc_x = new HierarchyVector("cc_x", fixture.hierarchy, 0, 0);
        Pointer<HierarchyVector> sc_x = new HierarchyVector("sc_x", fixture.hierarchy, 0, 0);
        cc_x->addComponent(fixture.x->getComponentVariable(1), fixture.x->getComponentDescriptorIndex(1));
        sc_x->addComponent(fixture.x->getComponentVariable(0), fixture.x->getComponentDescriptorIndex(0));
        auto cc_b = cc_x->cloneVector("cc_b"), sc_b = sc_x->cloneVector("sc_b");
        cc_b->allocateVectorData();
        sc_b->allocateVectorData();
        PoissonSpecifications coefs("state_coefs");
        coefs.setCConstant(2.0);
        coefs.setDConstant(0.0);
        PetscInt cc_previous = 0, sc_previous = 0, stokes_previous = 0;
        for (int width : { 0, 2 })
        {
            auto db = level_solver_database("shell", width);
            // Exercise both existing shell compositions across these lifetimes.
            db->putString("shell_pc_type", width == 0 ? "multiplicative" : "additive");
            LevelSolverProbe<CCPoissonPETScLevelSolver> cc("state_cc", db);
            LevelSolverProbe<SCPoissonPETScLevelSolver> sc("state_sc", db);
            cc.setPoissonSpecifications(coefs);
            sc.setPoissonSpecifications(coefs);
            cc.setPhysicalBcCoef(nullptr);
            sc.setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr));
            cc.setTimeInterval(0.0, 1.0);
            sc.setTimeInterval(0.0, 1.0);
            PetscInt cc_total, sc_total, stokes_total;
            cc_valid = check_shell_state(cc, cc_x, cc_b, cc_total) && cc_valid;
            sc_valid = check_shell_state(sc, sc_x, sc_b, sc_total) && sc_valid;
            LevelSolverProbe<StaggeredStokesPETScLevelSolver> stokes("state_stokes", db);
            Mat creator = level_test_matrix(fixture.full_size, 4.0);
            stokes.setOperatorMat(creator);
            stokes.setTimeInterval(0.0, 1.0);
            stokes_valid = check_shell_state(stokes, fixture.x, fixture.b, stokes_total) && stokes_valid;
            stokes.setOperatorMat(nullptr);
            PetscErrorCode ierr = MatDestroy(&creator);
            IBTK_CHKERRQ(ierr);
            if (width == 2)
            {
                cc_valid = cc_valid && cc_total > cc_previous;
                sc_valid = sc_valid && sc_total > sc_previous;
                stokes_valid = stokes_valid && stokes_total > stokes_previous;
            }
            cc_previous = cc_total;
            sc_previous = sc_total;
            stokes_previous = stokes_total;
        }
        free_vector_components(*cc_b);
        free_vector_components(*sc_b);
    }
    for (int variant = 0; variant < 3; ++variant)
    {
        const bool full = variant != 2;
        LevelFixture fixture(app->getComponentDatabase("CartesianGeometry"), variant == 0 ? 0 : 1, full);
        LevelSolverProbe<StaggeredStokesPETScLevelSolver> solver("domain_stokes", level_solver_database());
        PoissonSpecifications coefs("domain_coefs");
        coefs.setCConstant(2.0);
        coefs.setDConstant(-0.01);
        solver.setVelocityPoissonSpecifications(coefs);
        solver.setPhysicalBcCoefs(std::vector<RobinBcCoefStrategy<NDIM>*>(NDIM, nullptr), nullptr);
        solver.setComponentsHaveNullSpace(false, true);
        solver.setTimeInterval(0.0, 1.0);
        solver.setSolutionTime(1.0);
        solver.initializeSolverState(*fixture.x, *fixture.b);
        Mat matrix;
        MatNullSpace nullspace;
        PetscErrorCode ierr = KSPGetOperators(solver.getPETScKSP(), &matrix, nullptr);
        IBTK_CHKERRQ(ierr);
        ierr = MatGetNullSpace(matrix, &nullspace);
        IBTK_CHKERRQ(ierr);
        domain_valid = domain_valid && (full ? nullspace != nullptr : nullspace == nullptr);
        if (full && nullspace)
        {
            PetscBool is_nullspace;
            ierr = MatNullSpaceTest(nullspace, matrix, &is_nullspace);
            IBTK_CHKERRQ(ierr);
            domain_valid = domain_valid && is_nullspace;
        }
        solver.deallocateSolverState();
    }
    int failures = !cc_valid + !sc_valid + !stokes_valid + !domain_valid;
    pout << "shell_solve_error = " << max_level_solve_error << '\n'
         << "cc_state_valid = " << (cc_valid ? "true" : "false") << '\n'
         << "sc_state_valid = " << (sc_valid ? "true" : "false") << '\n'
         << "stokes_state_valid = " << (stokes_valid ? "true" : "false") << '\n'
         << "domain_nullspace_valid = " << (domain_valid ? "true" : "false") << '\n'
         << "test_failures = " << failures << std::endl;
    return failures;
}
} // namespace

namespace
{
#if defined(PETSC_USE_LOG)
PetscInt petsc_vec_creation_count = 0;
PetscInt petsc_vec_destruction_count = 0;

#if PETSC_VERSION_GE(3, 20, 0)
PetscErrorCode
ignore_petsc_log_event(PetscLogEvent, int, PetscObject, PetscObject, PetscObject, PetscObject)
{
    return 0;
}
#endif

PetscErrorCode
count_petsc_vec_creation(PetscObject object)
{
    PetscClassId class_id;
    const PetscErrorCode ierr = PetscObjectGetClassId(object, &class_id);
    if (ierr) return ierr;
    if (class_id == VEC_CLASSID) ++petsc_vec_creation_count;
    return 0;
}

PetscErrorCode
count_petsc_vec_destruction(PetscObject object)
{
    PetscClassId class_id;
    const PetscErrorCode ierr = PetscObjectGetClassId(object, &class_id);
    if (ierr) return ierr;
    if (class_id == VEC_CLASSID) ++petsc_vec_destruction_count;
    return 0;
}
#endif

struct StructureSpec
{
    int num_curve_points = 64;
    double ds = 1.0 / 64.0;
    double x_center = 0.5;
    double y_center = 0.5;
    double x_radius = 0.2;
    double y_radius = 0.2;
    double spring_stiffness = 2.0e2;
    int finest_ln = 0;
};

void
generate_structure(const unsigned int& strct_num,
                   const int& ln,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn,
                   void* ctx)
{
    auto* spec = static_cast<StructureSpec*>(ctx);
    if (!spec)
    {
        TBOX_ERROR("generate_structure(): missing structure specification context\n");
    }

    if (ln != spec->finest_ln || strct_num != 0)
    {
        num_vertices = 0;
        vertex_posn.resize(0);
        return;
    }

    num_vertices = spec->num_curve_points;
    vertex_posn.resize(num_vertices);
    for (int k = 0; k < num_vertices; ++k)
    {
        const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_vertices);
        vertex_posn[k](0) = spec->x_center + spec->x_radius * std::cos(theta);
        vertex_posn[k](1) = spec->y_center + spec->y_radius * std::sin(theta);
    }
}

void
generate_springs(
    const unsigned int& strct_num,
    const int& ln,
    std::multimap<int, IBRedundantInitializer::Edge>& spring_map,
    std::map<IBRedundantInitializer::Edge, IBRedundantInitializer::SpringSpec, IBRedundantInitializer::EdgeComp>&
        spring_spec,
    void* ctx)
{
    auto* spec = static_cast<StructureSpec*>(ctx);
    if (!spec)
    {
        TBOX_ERROR("generate_springs(): missing structure specification context\n");
    }
    if (ln != spec->finest_ln || strct_num != 0) return;

    for (int k = 0; k < spec->num_curve_points; ++k)
    {
        IBRedundantInitializer::Edge edge = { k, (k + 1) % spec->num_curve_points };
        if (edge.first > edge.second) std::swap(edge.first, edge.second);
        spring_map.insert(std::make_pair(edge.first, edge));

        IBRedundantInitializer::SpringSpec spec_data;
        spec_data.force_fcn_idx = 0;
        spec_data.parameters.resize(2);
        spec_data.parameters[0] = spec->spring_stiffness;
        spec_data.parameters[1] = 0.0;
        spring_spec.insert(std::make_pair(edge, spec_data));
    }
}

bool
side_l2_norm_is_finite(const Pointer<HierarchySideDataOpsReal<NDIM, double>>& side_data_ops,
                       const int data_idx,
                       const int weight_idx,
                       double& l2_norm)
{
    l2_norm = side_data_ops->L2Norm(data_idx, weight_idx);
    return std::isfinite(l2_norm);
}

bool
cell_l2_norm_is_finite(const Pointer<HierarchyCellDataOpsReal<NDIM, double>>& cell_data_ops,
                       const int data_idx,
                       const int weight_idx,
                       double& l2_norm)
{
    l2_norm = cell_data_ops->L2Norm(data_idx, weight_idx);
    return std::isfinite(l2_norm);
}

void
set_divergence_free_probe_velocity(const int u_idx, Pointer<PatchHierarchy<NDIM>> patch_hierarchy)
{
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();
            const double* const x_lower = patch_geom->getXLower();
            const SAMRAI::hier::Index<NDIM>& lower = patch_box.lower();

            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
                for (Box<NDIM>::Iterator b(side_box); b; b++)
                {
                    const SideIndex<NDIM> i_s(b(), axis, SideIndex<NDIM>::Lower);
                    double x[NDIM];
                    for (int d = 0; d < NDIM; ++d)
                    {
                        const double offset = (d == axis) ? 0.0 : 0.5;
                        x[d] = x_lower[d] + (static_cast<double>(i_s(d) - lower(d)) + offset) * dx[d];
                    }
                    const double val = (axis == 0) ? std::sin(2.0 * M_PI * x[1]) : -std::sin(2.0 * M_PI * x[0]);
                    (*u_data)(i_s) = val;
                }
            }
        }
    }
}

int
run_foundation(Pointer<AppInitializer> app_initializer)
{
    int test_failures = 0;
    {
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        const double current_time = 0.0;
        const double dt = input_db->getDoubleWithDefault("DT", 0.005);
        const double rho = input_db->getDoubleWithDefault("RHO", 1.0);
        const double mu = input_db->getDoubleWithDefault("MU", 1.0);
        const double new_time = current_time + dt;
        const bool use_fixed_le_operators = input_db->getBoolWithDefault("USE_FIXED_LE_OPERATORS", true);
        StructureSpec structure_spec;
        structure_spec.ds = input_db->getDoubleWithDefault("DS", 1.0 / 64.0);
        structure_spec.x_center = input_db->getDoubleWithDefault("X_CENTER", 0.5);
        structure_spec.y_center = input_db->getDoubleWithDefault("Y_CENTER", 0.5);
        structure_spec.x_radius = input_db->getDoubleWithDefault("X_RADIUS", 0.2);
        structure_spec.y_radius = input_db->getDoubleWithDefault("Y_RADIUS", 0.2);
        structure_spec.spring_stiffness = input_db->getDoubleWithDefault("SPRING_STIFFNESS", 2.0e2);
        if (!(structure_spec.ds > 0.0))
        {
            TBOX_ERROR("DS must be positive\n");
        }
        if (!(structure_spec.x_radius > 0.0) || !(structure_spec.y_radius > 0.0))
        {
            TBOX_ERROR("X_RADIUS and Y_RADIUS must be positive\n");
        }
        const double a = structure_spec.x_radius;
        const double b = structure_spec.y_radius;
        const double h = std::pow((a - b) / (a + b), 2.0);
        const double circumference =
            M_PI * (a + b) * (1.0 + 3.0 * h / (10.0 + std::sqrt(std::max(0.0, 4.0 - 3.0 * h))));
        structure_spec.num_curve_points = std::max(3, static_cast<int>(circumference / structure_spec.ds));
        if (structure_spec.num_curve_points < 3)
        {
            TBOX_ERROR("computed num_curve_points must be >= 3\n");
        }

        Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
        ib_method_ops->setUseFixedLEOperators(true);

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               ib_method_ops,
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

        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        structure_spec.finest_ln = input_db->getIntegerWithDefault("MAX_LEVELS", 1) - 1;
        ib_initializer->setStructureNamesOnLevel(structure_spec.finest_ln, { "curve2d" });
        ib_initializer->registerInitStructureFunction(generate_structure, &structure_spec);
        ib_initializer->registerInitSpringDataFunction(generate_springs, &structure_spec);
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, current_time);
        int tag_buffer = input_db->getIntegerWithDefault("TAG_BUFFER", 1);
        int level_number = 0;
        bool done = false;
        while (!done && gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, current_time, true, tag_buffer);
            done = !patch_hierarchy->finerLevelExists(level_number);
            ++level_number;
        }

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> current_ctx = var_db->getContext("current_ctx");
        Pointer<VariableContext> scratch_ctx = var_db->getContext("scratch_ctx");
        Pointer<VariableContext> solver_ctx = var_db->getContext("solver_ctx");

        Pointer<SideVariable<NDIM, double>> f_var = new SideVariable<NDIM, double>("f_var");
        Pointer<CellVariable<NDIM, double>> g_var = new CellVariable<NDIM, double>("g_var");
        Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p_var");
        Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("p_dof_index");
        Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u_var");
        Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("u_dof_index");

        const IntVector<NDIM> ib_ghosts = ib_method_ops->getMinimumGhostCellWidth();
        const IntVector<NDIM> one_ghost = IntVector<NDIM>(1);
        const IntVector<NDIM> no_ghosts = IntVector<NDIM>(0);

        const int u_current_idx = var_db->registerVariableAndContext(u_var, current_ctx, ib_ghosts);
        const int u_sol_idx = var_db->registerVariableAndContext(u_var, solver_ctx, one_ghost);
        const int f_rhs_idx = var_db->registerVariableAndContext(f_var, solver_ctx, one_ghost);
        const int p_sol_idx = var_db->registerVariableAndContext(p_var, solver_ctx, one_ghost);
        const int g_rhs_idx = var_db->registerVariableAndContext(g_var, solver_ctx, one_ghost);
        const int u_scratch_idx = var_db->registerVariableAndContext(u_var, scratch_ctx, ib_ghosts);
        const int f_scratch_idx = var_db->registerVariableAndContext(f_var, scratch_ctx, ib_ghosts);
        const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, scratch_ctx, ib_ghosts);
        const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, scratch_ctx, no_ghosts);

        const std::vector<int> allocated_patch_data_indices = {
            u_current_idx, u_scratch_idx, f_scratch_idx, u_dof_index_idx, p_dof_index_idx
        };
        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int data_idx : allocated_patch_data_indices) level->allocatePatchData(data_idx, current_time);
        }

        Pointer<HierarchySideDataOpsReal<NDIM, double>> hier_velocity_data_ops =
            new HierarchySideDataOpsReal<NDIM, double>(patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());
        Pointer<HierarchyCellDataOpsReal<NDIM, double>> hier_pressure_data_ops =
            new HierarchyCellDataOpsReal<NDIM, double>(patch_hierarchy, 0, patch_hierarchy->getFinestLevelNumber());

        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            u_init->setDataOnPatchHierarchy(u_current_idx, u_var, patch_hierarchy, current_time);
        }
        else
        {
            hier_velocity_data_ops->setToScalar(u_current_idx, 0.0, false);
        }
        hier_velocity_data_ops->setToScalar(u_scratch_idx, 0.0, false);
        hier_velocity_data_ops->setToScalar(f_scratch_idx, 0.0, false);

        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        std::vector<Pointer<CoarsenSchedule<NDIM>>> u_synch_scheds(finest_ln + 1);
        std::vector<Pointer<RefineSchedule<NDIM>>> u_ghost_fill_scheds(finest_ln + 1);
        std::vector<Pointer<RefineSchedule<NDIM>>> f_prolongation_scheds(finest_ln + 1);

        ib_method_ops->initializePatchHierarchy(patch_hierarchy,
                                                gridding_algorithm,
                                                u_current_idx,
                                                u_synch_scheds,
                                                u_ghost_fill_scheds,
                                                0,
                                                current_time,
                                                true);
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();

        ib_method_ops->preprocessIntegrateData(current_time, new_time, /*num_cycles*/ 1);
        ib_method_ops->updateFixedLEOperators();

        std::vector<std::vector<int>> num_dofs_per_proc(finest_ln + 1);
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
                num_dofs_per_proc[ln], u_dof_index_idx, p_dof_index_idx, level);
        }

        Mat A = nullptr;
        ib_method_ops->constructLagrangianForceJacobian(A, MATAIJ, new_time);
        Mat J = nullptr;
        ib_method_ops->constructInterpOp(J, IBKernel::IB_4, num_dofs_per_proc[finest_ln], u_dof_index_idx, new_time);

        HierarchyMathOps hier_math_ops("hier_math_ops", patch_hierarchy);
        const int wgt_cc_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();
        const int wgt_sc_idx = hier_math_ops.getSideWeightPatchDescriptorIndex();

        Pointer<SAMRAIVectorReal<NDIM, double>> eul_sol_vec =
            new SAMRAIVectorReal<NDIM, double>("eul_sol_vec", patch_hierarchy, 0, finest_ln);
        eul_sol_vec->addComponent(u_var, u_sol_idx, wgt_sc_idx, hier_velocity_data_ops);
        eul_sol_vec->addComponent(p_var, p_sol_idx, wgt_cc_idx, hier_pressure_data_ops);
        eul_sol_vec->allocateVectorData();

        Pointer<SAMRAIVectorReal<NDIM, double>> eul_rhs_vec =
            new SAMRAIVectorReal<NDIM, double>("eul_rhs_vec", patch_hierarchy, 0, finest_ln);
        eul_rhs_vec->addComponent(f_var, f_rhs_idx, wgt_sc_idx, hier_velocity_data_ops);
        eul_rhs_vec->addComponent(g_var, g_rhs_idx, wgt_cc_idx, hier_pressure_data_ops);
        eul_rhs_vec->allocateVectorData();

        hier_velocity_data_ops->copyData(u_sol_idx, u_current_idx);
        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            p_init->setDataOnPatchHierarchy(p_sol_idx, p_var, patch_hierarchy, current_time);
        }
        else
        {
            hier_pressure_data_ops->setToScalar(p_sol_idx, 0.0, false);
        }
        hier_velocity_data_ops->setToScalar(f_rhs_idx, 0.0, false);
        hier_pressure_data_ops->setToScalar(g_rhs_idx, 0.0, false);

        const double lambda = 0.0;
        PoissonSpecifications U_problem_coefs("stokes_ib_solver_components::U_problem_coefs");
        U_problem_coefs.setCConstant(rho / dt + lambda);
        U_problem_coefs.setDConstant(-mu);

        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);
        if (periodic_shift.min() <= 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                const std::string bc_coefs_name = "u_bc_coefs_" + std::to_string(d);
                const std::string bc_coefs_db_name = "VelocityBcCoefs_" + std::to_string(d);
                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
        }

        Pointer<StaggeredStokesOperator> stokes_op =
            new StaggeredStokesOperator("stokes_ib_solver_components::stokes_op", false);
        stokes_op->setVelocityPoissonSpecifications(U_problem_coefs);
        stokes_op->setPhysicalBcCoefs(u_bc_coefs, nullptr);
        stokes_op->setTimeInterval(current_time, new_time);
        stokes_op->setSolutionTime(new_time);

        StaggeredStokesIBOperator::Context ctx;
        ctx.ib_implicit_ops = ib_method_ops;
        ctx.stokes_op = stokes_op;
        ctx.u_phys_bdry_op = nullptr;
        ctx.hier_velocity_data_ops = hier_velocity_data_ops;
        ctx.u_synch_scheds = u_synch_scheds;
        ctx.u_ghost_fill_scheds = u_ghost_fill_scheds;
        ctx.f_prolongation_scheds = f_prolongation_scheds;
        ctx.patch_level = patch_hierarchy->getPatchLevel(finest_ln);
        ctx.u_idx = u_scratch_idx;
        ctx.f_idx = f_scratch_idx;
        ctx.u_current_idx = u_current_idx;
        ctx.u_dof_index_idx = u_dof_index_idx;
        ctx.p_dof_index_idx = p_dof_index_idx;
        ctx.use_fixed_le_operators = use_fixed_le_operators;
        ctx.time_stepping_type = IBAMR::string_to_enum<TimeSteppingType>(
            input_db->getStringWithDefault("IB_TIME_STEPPING", "MIDPOINT_RULE"));

        StaggeredStokesIBOperator nonlinear_op("stokes_ib_solver_components::nonlinear_op", false);
        nonlinear_op.setOperatorContext(ctx);
        nonlinear_op.setTimeInterval(current_time, new_time);
        nonlinear_op.setSolutionTime(new_time);
        nonlinear_op.initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);

        Pointer<SAMRAIVectorReal<NDIM, double>> nonlinear_probe = eul_sol_vec->cloneVector("nonlinear_probe");
        nonlinear_probe->allocateVectorData();
        nonlinear_probe->setToScalar(0.0);
        set_divergence_free_probe_velocity(nonlinear_probe->getComponentDescriptorIndex(0), patch_hierarchy);

        Pointer<SAMRAIVectorReal<NDIM, double>> f_probe = eul_rhs_vec->cloneVector("f_probe");
        f_probe->allocateVectorData();
        f_probe->setToScalar(0.0);
        nonlinear_op.apply(*nonlinear_probe, *f_probe);

        double nonlinear_side_norm = std::numeric_limits<double>::quiet_NaN();
        double nonlinear_cell_norm = std::numeric_limits<double>::quiet_NaN();
        const bool expect_trivial_nonlinear = (std::abs(structure_spec.spring_stiffness) <= 1.0e-14) &&
                                              (std::abs(rho) <= 1.0e-14) && (std::abs(mu) <= 1.0e-14);
        if (!side_l2_norm_is_finite(
                hier_velocity_data_ops, f_probe->getComponentDescriptorIndex(0), wgt_sc_idx, nonlinear_side_norm) ||
            !cell_l2_norm_is_finite(
                hier_pressure_data_ops, f_probe->getComponentDescriptorIndex(1), wgt_cc_idx, nonlinear_cell_norm))
        {
            ++test_failures;
            pout << "nonlinear operator produced non-finite norm" << std::endl;
        }
        else if (nonlinear_side_norm <= 1.0e-14 && nonlinear_cell_norm <= 1.0e-14)
        {
            if (!expect_trivial_nonlinear)
            {
                ++test_failures;
                pout << "nonlinear operator action is trivial" << std::endl;
            }
        }
        else if (expect_trivial_nonlinear)
        {
            ++test_failures;
            pout << "nonlinear operator action is nontrivial when SPRING_STIFFNESS, RHO, and MU are zero" << std::endl;
        }

        Pointer<StaggeredStokesIBJacobianOperator> jac_op =
            new StaggeredStokesIBJacobianOperator("stokes_ib_solver_components::jacobian_op");
        jac_op->setOperatorContext(ctx);
        jac_op->setTimeInterval(current_time, new_time);
        jac_op->setSolutionTime(new_time);
        jac_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);
        jac_op->formJacobian(*eul_sol_vec);

        Pointer<SAMRAIVectorReal<NDIM, double>> v = eul_sol_vec->cloneVector("v");
        v->allocateVectorData();
        v->setToScalar(0.0);
        hier_velocity_data_ops->setToScalar(v->getComponentDescriptorIndex(0), 1.0, false);
        hier_pressure_data_ops->setToScalar(v->getComponentDescriptorIndex(1), -0.25, false);

        Pointer<SAMRAIVectorReal<NDIM, double>> jv = eul_rhs_vec->cloneVector("jv");
        jv->allocateVectorData();
        jv->setToScalar(0.0);
        jac_op->apply(*v, *jv);

        const double fd_rel_tol = input_db->getDoubleWithDefault("FD_REL_TOL", 5.0e-2);
        Pointer<PETScMFFDJacobianOperator> mffd_jac_op =
            new PETScMFFDJacobianOperator("stokes_ib_solver_components::mffd_jacobian_op", "ib_jac_mffd_");
        mffd_jac_op->setOperator(Pointer<GeneralOperator>(&nonlinear_op, false));
        mffd_jac_op->setTimeInterval(current_time, new_time);
        mffd_jac_op->setSolutionTime(new_time);
        mffd_jac_op->initializeOperatorState(*eul_sol_vec, *eul_rhs_vec);
        mffd_jac_op->formJacobian(*eul_sol_vec);

        Pointer<SAMRAIVectorReal<NDIM, double>> fd_jv = eul_rhs_vec->cloneVector("fd_jv");
        fd_jv->allocateVectorData();
        fd_jv->setToScalar(0.0);
        mffd_jac_op->apply(*v, *fd_jv);

        Pointer<SAMRAIVectorReal<NDIM, double>> diff = eul_rhs_vec->cloneVector("diff");
        diff->allocateVectorData();
        diff->subtract(fd_jv, jv);

        double jv_side_norm = std::numeric_limits<double>::quiet_NaN();
        double jv_cell_norm = std::numeric_limits<double>::quiet_NaN();
        double diff_side_norm = std::numeric_limits<double>::quiet_NaN();
        double diff_cell_norm = std::numeric_limits<double>::quiet_NaN();
        const bool jv_finite =
            side_l2_norm_is_finite(
                hier_velocity_data_ops, jv->getComponentDescriptorIndex(0), wgt_sc_idx, jv_side_norm) &&
            cell_l2_norm_is_finite(
                hier_pressure_data_ops, jv->getComponentDescriptorIndex(1), wgt_cc_idx, jv_cell_norm);
        const bool diff_finite =
            side_l2_norm_is_finite(
                hier_velocity_data_ops, diff->getComponentDescriptorIndex(0), wgt_sc_idx, diff_side_norm) &&
            cell_l2_norm_is_finite(
                hier_pressure_data_ops, diff->getComponentDescriptorIndex(1), wgt_cc_idx, diff_cell_norm);
        if (!jv_finite || !diff_finite)
        {
            ++test_failures;
            pout << "jacobian norms are non-finite" << std::endl;
        }
        else if (jv_side_norm <= 1.0e-14 && jv_cell_norm <= 1.0e-14)
        {
            ++test_failures;
            pout << "jacobian action is trivial" << std::endl;
        }
        else
        {
            const double rel_error =
                std::sqrt(diff_side_norm * diff_side_norm + diff_cell_norm * diff_cell_norm) /
                std::max(std::sqrt(jv_side_norm * jv_side_norm + jv_cell_norm * jv_cell_norm), 1.0e-14);
            const bool fd_relative_error_valid = rel_error <= fd_rel_tol;
            pout << "fd_relative_error_valid = " << (fd_relative_error_valid ? "true" : "false") << std::endl;
            if (!fd_relative_error_valid)
            {
                ++test_failures;
                pout << "fd_relative_error exceeds tolerance: " << fd_rel_tol << std::endl;
            }
        }
        mffd_jac_op->deallocateOperatorState();

        Pointer<Database> stokes_ib_precond_db =
            input_db->isDatabase("stokes_ib_precond_db") ? input_db->getDatabase("stokes_ib_precond_db") : nullptr;

        Pointer<StaggeredStokesIBLevelRelaxationFACOperator> fac_op = new StaggeredStokesIBLevelRelaxationFACOperator(
            "stokes_ib_solver_components::fac_op", stokes_ib_precond_db, "stokes_ib_pc_");
        Pointer<StaggeredStokesIBJacobianFACPreconditioner> fac_pc = new StaggeredStokesIBJacobianFACPreconditioner(
            "stokes_ib_solver_components::fac_pc", fac_op, stokes_ib_precond_db, "stokes_ib_pc_");
        Pointer<StaggeredStokesPhysicalBoundaryHelper> bc_helper = new StaggeredStokesPhysicalBoundaryHelper();

        fac_pc->setVelocityPoissonSpecifications(U_problem_coefs);
        fac_pc->setPhysicalBcCoefs(u_bc_coefs, nullptr);
        fac_pc->setPhysicalBoundaryHelper(bc_helper);
        fac_pc->setTimeInterval(current_time, new_time);
        fac_pc->setSolutionTime(new_time);
        fac_pc->setHomogeneousBc(true);
        fac_pc->setComponentsHaveNullSpace(false, true);
        fac_pc->setIBTimeSteppingType(ctx.time_stepping_type);
        fac_pc->setIBForceJacobian(A);
        fac_pc->setIBInterpOp(J);
        fac_pc->setIBImplicitStrategy(ib_method_ops);
        fac_pc->initializeSolverState(*eul_sol_vec, *eul_rhs_vec);

        const bool verify_galerkin_operator_borrowing =
            input_db->getBoolWithDefault("VERIFY_GALERKIN_OPERATOR_BORROWING", false);
        auto check_fac_residual_work_vector_cache = [&](double& reuse_error)
        {
            Pointer<SAMRAIVectorReal<NDIM, double>> first_residual =
                eul_rhs_vec->cloneVector("fac_residual_workspace_first");
            Pointer<SAMRAIVectorReal<NDIM, double>> second_residual =
                eul_rhs_vec->cloneVector("fac_residual_workspace_second");
            first_residual->allocateVectorData();
            second_residual->allocateVectorData();

#if defined(PETSC_USE_LOG)
            petsc_vec_creation_count = 0;
            petsc_vec_destruction_count = 0;
#if PETSC_VERSION_GE(3, 20, 0)
            PetscLogHandler log_handler = nullptr;
            PetscErrorCode log_ierr = PetscLogHandlerCreateLegacy(PETSC_COMM_WORLD,
                                                                  ignore_petsc_log_event,
                                                                  ignore_petsc_log_event,
                                                                  count_petsc_vec_creation,
                                                                  count_petsc_vec_destruction,
                                                                  &log_handler);
            IBTK_CHKERRQ(log_ierr);
            log_ierr = PetscLogHandlerStart(log_handler);
            IBTK_CHKERRQ(log_ierr);
#else
            // Older PETSc versions expose the object callbacks directly.
            const auto saved_create = PetscLogPHC;
            const auto saved_destroy = PetscLogPHD;
            PetscLogPHC = count_petsc_vec_creation;
            PetscLogPHD = count_petsc_vec_destruction;
#endif
            // Verify that logging is active before trusting zero event counts.
            Vec logging_probe = nullptr;
            PetscErrorCode ierr = VecCreate(PETSC_COMM_WORLD, &logging_probe);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&logging_probe);
            IBTK_CHKERRQ(ierr);
            const bool logging_active = petsc_vec_creation_count == 1 && petsc_vec_destruction_count == 1;
            if (!logging_active) pout << "PETSc Vec logging calibration failed" << std::endl;
            petsc_vec_creation_count = 0;
            petsc_vec_destruction_count = 0;
#endif
            fac_op->computeResidual(*first_residual, *nonlinear_probe, *jv, 0, finest_ln);
            fac_op->computeResidual(*second_residual, *nonlinear_probe, *jv, 0, finest_ln);
#if defined(PETSC_USE_LOG)
#if PETSC_VERSION_GE(3, 20, 0)
            ierr = PetscLogHandlerStop(log_handler);
            IBTK_CHKERRQ(ierr);
            ierr = PetscLogHandlerDestroy(&log_handler);
            IBTK_CHKERRQ(ierr);
#else
            PetscLogPHC = saved_create;
            PetscLogPHD = saved_destroy;
#endif
            if (!logging_active || petsc_vec_creation_count != 0 || petsc_vec_destruction_count != 0)
            {
                ++test_failures;
                pout << "FAC residual Vec allocation check failed: creations = " << petsc_vec_creation_count
                     << ", destructions = " << petsc_vec_destruction_count << std::endl;
            }
#endif
            const double residual_norm = first_residual->maxNorm();
            first_residual->subtract(first_residual, second_residual);
            reuse_error = std::abs(first_residual->maxNorm());
            free_vector_components(*first_residual);
            free_vector_components(*second_residual);
            return std::isfinite(residual_norm) && residual_norm > 1.0e-12 && std::isfinite(reuse_error) &&
                   reuse_error == 0.0;
        };

        double fac_residual_work_vector_reuse_error = 0.0;
        const bool fac_residual_repeat_valid =
            check_fac_residual_work_vector_cache(fac_residual_work_vector_reuse_error);
        if (!fac_residual_repeat_valid) ++test_failures;

        bool galerkin_operator_available_valid = !verify_galerkin_operator_borrowing;
        bool galerkin_operator_creator_lifetime_valid = !verify_galerkin_operator_borrowing;
        if (verify_galerkin_operator_borrowing)
        {
            Pointer<StaggeredStokesPETScLevelSolver> coarse_level_solver =
                fac_op->getStaggeredStokesPETScLevelSolver(0);
            Mat supplied_operator = nullptr;
            PetscErrorCode ierr = KSPGetOperators(coarse_level_solver->getPETScKSP(), &supplied_operator, nullptr);
            IBTK_CHKERRQ(ierr);
            PetscInt reference_count_before = 0;
            ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(supplied_operator), &reference_count_before);
            IBTK_CHKERRQ(ierr);
            galerkin_operator_available_valid = supplied_operator != nullptr && reference_count_before > 1;
            coarse_level_solver->deallocateSolverState();
            PetscInt n_rows = 0;
            PetscInt n_columns = 0;
            ierr = MatGetSize(supplied_operator, &n_rows, &n_columns);
            IBTK_CHKERRQ(ierr);
            PetscInt reference_count_after = 0;
            ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(supplied_operator), &reference_count_after);
            IBTK_CHKERRQ(ierr);
            // Solver-state teardown retains the installed input in addition to
            // the FAC creator's reference. Explicit clearing releases only the former.
            galerkin_operator_creator_lifetime_valid = n_rows > 0 && n_rows == n_columns && reference_count_after == 2;
            coarse_level_solver->setOperatorMat(nullptr);
            galerkin_operator_creator_lifetime_valid =
                matrix_references(supplied_operator) == 1 && galerkin_operator_creator_lifetime_valid;
            if (!galerkin_operator_available_valid || !galerkin_operator_creator_lifetime_valid) ++test_failures;
        }

        bool fac_reinitialization_valid = true;
        double fac_residual_work_vector_reinitialize_error = 0.0;
        fac_pc->deallocateSolverState();
        fac_pc->initializeSolverState(*eul_sol_vec, *eul_rhs_vec);
        const bool fac_residual_repeat_reinitialize_valid =
            check_fac_residual_work_vector_cache(fac_residual_work_vector_reinitialize_error);
        if (verify_galerkin_operator_borrowing)
        {
            Pointer<StaggeredStokesPETScLevelSolver> coarse_level_solver =
                fac_op->getStaggeredStokesPETScLevelSolver(0);
            Mat supplied_operator = nullptr;
            PetscErrorCode ierr = KSPGetOperators(coarse_level_solver->getPETScKSP(), &supplied_operator, nullptr);
            IBTK_CHKERRQ(ierr);
            PetscInt reference_count = 0;
            ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(supplied_operator), &reference_count);
            IBTK_CHKERRQ(ierr);
            fac_reinitialization_valid = supplied_operator != nullptr && reference_count > 1;
        }
        if (!fac_residual_repeat_reinitialize_valid || !fac_reinitialization_valid) ++test_failures;
        Mat SAJ = fac_op->getEulerianElasticityLevelOp(finest_ln);
        jac_op->setIBCouplingJacobian(SAJ);

        Pointer<PETScKrylovLinearSolver> linear_solver =
            new PETScKrylovLinearSolver("stokes_ib_solver_components::linear_solver", nullptr, "ib_");
        linear_solver->setOperator(jac_op);
        linear_solver->setPreconditioner(fac_pc);
        linear_solver->setTimeInterval(current_time, new_time);
        linear_solver->setSolutionTime(new_time);
        linear_solver->setInitialGuessNonzero(false);
        // Richardson self-scaling makes the FAC action nonlinear. Flexible GMRES
        // uses a right preconditioner and measures the unpreconditioned residual.
        linear_solver->setKSPType("fgmres");

        Pointer<SAMRAIVectorReal<NDIM, double>> linear_sol = eul_sol_vec->cloneVector("linear_sol");
        linear_sol->allocateVectorData();
        linear_sol->setToScalar(0.0);
        linear_solver->initializeSolverState(*linear_sol, *jv);
        PetscErrorCode linear_ierr = KSPSetPCSide(linear_solver->getPETScKSP(), PC_RIGHT);
        IBTK_CHKERRQ(linear_ierr);
        linear_ierr = KSPSetNormType(linear_solver->getPETScKSP(), KSP_NORM_UNPRECONDITIONED);
        IBTK_CHKERRQ(linear_ierr);
        const bool linear_success = linear_solver->solveSystem(*linear_sol, *jv);
        if (!linear_success)
        {
            ++test_failures;
            pout << "krylov linear solve failed" << std::endl;
        }
        pout << "krylov_linear_iterations = " << linear_solver->getNumIterations() << std::endl;
        jac_op->apply(*linear_sol, *diff);
        diff->subtract(diff, jv);
        const double actual_residual = diff->L2Norm();
        const double residual_limit =
            std::max(linear_solver->getAbsoluteTolerance(), linear_solver->getRelativeTolerance() * jv->L2Norm());
        KSPConvergedReason reason;
        PCSide side;
        KSPNormType norm_type;
        KSPType ksp_type;
        linear_ierr = KSPGetConvergedReason(linear_solver->getPETScKSP(), &reason);
        IBTK_CHKERRQ(linear_ierr);
        linear_ierr = KSPGetPCSide(linear_solver->getPETScKSP(), &side);
        IBTK_CHKERRQ(linear_ierr);
        linear_ierr = KSPGetNormType(linear_solver->getPETScKSP(), &norm_type);
        IBTK_CHKERRQ(linear_ierr);
        linear_ierr = KSPGetType(linear_solver->getPETScKSP(), &ksp_type);
        IBTK_CHKERRQ(linear_ierr);
        const bool krylov_linear_residual_valid =
            std::string(ksp_type) == KSPFGMRES && side == PC_RIGHT && norm_type == KSP_NORM_UNPRECONDITIONED &&
            reason > 0 && std::isfinite(linear_solver->getResidualNorm()) && linear_solver->getResidualNorm() >= 0.0 &&
            std::isfinite(actual_residual) && actual_residual <= std::max(1.0e-12, 2.0 * residual_limit);
        if (!krylov_linear_residual_valid)
        {
            pout << "actual_residual = " << actual_residual
                 << ", acceptance_bound = " << std::max(1.0e-12, 2.0 * residual_limit)
                 << ", reported_residual = " << linear_solver->getResidualNorm() << ", pc_side = " << side
                 << ", norm_type = " << norm_type << ", ksp_type = " << ksp_type << ", reason = " << reason
                 << std::endl;
        }
        if (!krylov_linear_residual_valid) ++test_failures;
        pout << "krylov_linear_residual_valid = " << (krylov_linear_residual_valid ? "true" : "false") << std::endl;

        double linear_side_norm = std::numeric_limits<double>::quiet_NaN();
        double linear_cell_norm = std::numeric_limits<double>::quiet_NaN();
        if (!side_l2_norm_is_finite(
                hier_velocity_data_ops, linear_sol->getComponentDescriptorIndex(0), wgt_sc_idx, linear_side_norm) ||
            !cell_l2_norm_is_finite(
                hier_pressure_data_ops, linear_sol->getComponentDescriptorIndex(1), wgt_cc_idx, linear_cell_norm))
        {
            ++test_failures;
            pout << "krylov linear solve produced non-finite norm" << std::endl;
        }
        else if (linear_side_norm <= 1.0e-14 && linear_cell_norm <= 1.0e-14)
        {
            ++test_failures;
            pout << "krylov linear solve action is trivial" << std::endl;
        }

        linear_solver->deallocateSolverState();
        fac_pc->deallocateSolverState();

        if (verify_galerkin_operator_borrowing)
        {
            pout << "galerkin_operator_available_valid = " << (galerkin_operator_available_valid ? "true" : "false")
                 << std::endl;
            pout << "galerkin_operator_creator_lifetime_valid = "
                 << (galerkin_operator_creator_lifetime_valid ? "true" : "false") << std::endl;
            pout << "galerkin_operator_reinitialization_valid = " << (fac_reinitialization_valid ? "true" : "false")
                 << std::endl;
        }
        pout << "FAC Vec allocation measurement: calibrated and enforced with PETSC_USE_LOG; unavailable without it."
             << std::endl;
        pout << "fac_residual_work_vector_reuse_exact = "
             << (fac_residual_work_vector_reuse_error == 0.0 ? "true" : "false") << std::endl;
        pout << "fac_residual_repeat_valid = " << (fac_residual_repeat_valid ? "true" : "false") << std::endl;
        pout << "fac_residual_work_vector_reinitialize_exact = "
             << (fac_residual_work_vector_reinitialize_error == 0.0 ? "true" : "false") << std::endl;
        pout << "fac_residual_repeat_reinitialize_valid = "
             << (fac_residual_repeat_reinitialize_valid ? "true" : "false") << std::endl;

        jac_op->deallocateOperatorState();
        nonlinear_op.deallocateOperatorState();

        ib_method_ops->postprocessIntegrateData(current_time, new_time, /*num_cycles*/ 1);

        for (auto vec : { nonlinear_probe, f_probe, v, jv, fd_jv, diff, linear_sol }) free_vector_components(*vec);

        deallocate_vector_data(*eul_sol_vec);
        deallocate_vector_data(*eul_rhs_vec);
        free_vector_components(*eul_sol_vec);
        free_vector_components(*eul_rhs_vec);

        for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int data_idx : allocated_patch_data_indices)
            {
                if (level->checkAllocated(data_idx)) level->deallocatePatchData(data_idx);
            }
        }

        PetscErrorCode ierr = MatDestroy(&A);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&J);
        IBTK_CHKERRQ(ierr);

        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

        pout << "test_failures = " << test_failures << std::endl;
    }

    return test_failures;
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
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "components.log");
    const std::string test_case = app->getInputDatabase()->getStringWithDefault("test_case", "interpolation");
    if (test_case == "interpolation") return run_interpolation(app, input_file);
    if (test_case == "operators") return run_operators(app);
    if (test_case == "level_borrowing") return run_level_operator(app, false);
    if (test_case == "level_augmentation") return run_level_operator(app, true);
    if (test_case == "set_operator_initialized") return run_initialized_matrix_setter(app, false);
    if (test_case == "set_augmentation_initialized") return run_initialized_matrix_setter(app, true);
    if (test_case == "level_state") return run_level_state(app);
    if (test_case == "foundation") return run_foundation(app);
    TBOX_ERROR("Unknown component test case: " << test_case << '\n');
    return 1;
}
