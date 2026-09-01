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
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>

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
#include <map>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
using Kernel = void (*)(double, double*);
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
void
three_point_probe(double r, double* w)
{
    w[0] = 0.125 + r / 32.0;
    w[1] = 0.375 - r / 32.0;
    w[2] = 0.5;
}

int
check_enums_and_kernels()
{
    const std::array<DeltaFunctionType, 18> types = { PIECEWISE_CONSTANT,
                                                      PIECEWISE_LINEAR,
                                                      BSPLINE_3,
                                                      BSPLINE_4,
                                                      BSPLINE_5,
                                                      BSPLINE_6,
                                                      COMPOSITE_BSPLINE_23,
                                                      COMPOSITE_BSPLINE_32,
                                                      COMPOSITE_BSPLINE_34,
                                                      COMPOSITE_BSPLINE_43,
                                                      COMPOSITE_BSPLINE_45,
                                                      COMPOSITE_BSPLINE_54,
                                                      COMPOSITE_BSPLINE_56,
                                                      COMPOSITE_BSPLINE_65,
                                                      IB_3,
                                                      IB_4,
                                                      IB_5,
                                                      IB_6 };
    bool enums = true;
    for (const auto type : types)
        enums = enums && IBAMR::string_to_enum<DeltaFunctionType>(IBAMR::enum_to_string(type)) == type;
    enums = enums && BSPLINE_1 == PIECEWISE_CONSTANT && BSPLINE_2 == PIECEWISE_LINEAR &&
            IBAMR::string_to_enum<DeltaFunctionType>("bspline_1") == PIECEWISE_CONSTANT &&
            IBAMR::string_to_enum<DeltaFunctionType>("bspline_2") == PIECEWISE_LINEAR &&
            IBAMR::enum_to_string(BSPLINE_1) == "PIECEWISE_CONSTANT" &&
            IBAMR::enum_to_string(BSPLINE_2) == "PIECEWISE_LINEAR" &&
            IBAMR::string_to_enum<DeltaFunctionType>("invalid") == UNKNOWN_DELTA_FUNCTION_TYPE &&
            IBAMR::enum_to_string(UNKNOWN_DELTA_FUNCTION_TYPE) == "UNKNOWN_DELTA_FUNCTION_TYPE";

    struct Sample
    {
        Kernel kernel;
        int width;
        double r;
        std::array<double, 6> expected;
    };
    const double a = (2.0 - std::sqrt(2.0)) / 8.0, b = (2.0 + std::sqrt(2.0)) / 8.0;
    const double K6 = (59.0 - std::sqrt(261.0)) / 60.0;
    // Independent samples of the one-dimensional kernel definitions. The IB_5
    // values use lagrangian_ib_5_delta in lagrangian_delta.f.m4 evaluated at r-i;
    // in particular, both support endpoints vanish when r=2.5.
    const std::array<Sample, 10> samples = {
        { { PETScMatUtilities::piecewise_constant_delta_fcn, 1, -0.25, { 1 } },
          { PETScMatUtilities::piecewise_linear_delta_fcn, 2, 0.25, { 0.75, 0.25 } },
          { PETScMatUtilities::bspline_3_delta_fcn, 4, 1.5, { 0, 0.5, 0.5, 0 } },
          { PETScMatUtilities::bspline_4_delta_fcn, 4, 1.5, { 1.0 / 48, 23.0 / 48, 23.0 / 48, 1.0 / 48 } },
          { PETScMatUtilities::bspline_5_delta_fcn, 6, 2.5, { 0, 1.0 / 24, 11.0 / 24, 11.0 / 24, 1.0 / 24, 0 } },
          { PETScMatUtilities::bspline_6_delta_fcn,
            6,
            2.5,
            { 1.0 / 3840, 237.0 / 3840, 1682.0 / 3840, 1682.0 / 3840, 237.0 / 3840, 1.0 / 3840 } },
          { PETScMatUtilities::ib_3_delta_fcn, 4, 1.5, { 0, 0.5, 0.5, 0 } },
          { PETScMatUtilities::ib_4_delta_fcn, 4, 1.5, { a, b, b, a } },
          { PETScMatUtilities::ib_5_delta_fcn,
            6,
            2.5,
            { 0, 0.0612224005711746881, 0.438777599428825312, 0.438777599428825312, 0.0612224005711746881, 0 } },
          { PETScMatUtilities::ib_6_delta_fcn,
            6,
            3.0,
            { 0, -1.0 / 16 + K6 / 8, 0.25, 5.0 / 8 - K6 / 4, 0.25, -1.0 / 16 + K6 / 8 } } }
    };
    const std::array<int, 10> widths = { PETScMatUtilities::piecewise_constant_delta_stencil,
                                         PETScMatUtilities::piecewise_linear_delta_stencil,
                                         PETScMatUtilities::bspline_3_delta_stencil,
                                         PETScMatUtilities::bspline_4_delta_stencil,
                                         PETScMatUtilities::bspline_5_delta_stencil,
                                         PETScMatUtilities::bspline_6_delta_stencil,
                                         PETScMatUtilities::ib_3_delta_stencil,
                                         PETScMatUtilities::ib_4_delta_stencil,
                                         PETScMatUtilities::ib_5_delta_stencil,
                                         PETScMatUtilities::ib_6_delta_stencil };
    bool kernels = true;
    for (unsigned int n = 0; n < samples.size(); ++n)
    {
        const auto& sample = samples[n];
        std::array<double, 8> weights;
        weights.fill(-123.0);
        sample.kernel(sample.r, weights.data());
        kernels = kernels && sample.width == widths[n];
        for (int i = 0; i < sample.width; ++i)
            kernels = kernels && std::isfinite(weights[i]) && std::abs(weights[i] - sample.expected[i]) < 1.0e-12;
        for (unsigned int i = sample.width; i < weights.size(); ++i) kernels = kernels && weights[i] == -123.0;
    }
    // IB_5: pointwise Fortran definition, on both sides of the nearest-center
    // change. IB_6: lagrangian_ib_6_interp2d's pm3,...,pp2 recurrence with
    // ic_lower=0 and X/dx=r+0.5, hence its coordinate is 1-X/dx+2.5=3-r.
    // These values were evaluated independently at high precision; reversing
    // the off-center weight order must not pass as a symmetric-kernel check.
    const std::array<Sample, 4> off_center_samples = { { { PETScMatUtilities::ib_5_delta_fcn,
                                                           6,
                                                           2.25,
                                                           { 0.000539644595320609716,
                                                             0.128737522475479593,
                                                             0.514244366143986938,
                                                             0.333140121904304905,
                                                             0.0233383448809079538,
                                                             0 } },
                                                         { PETScMatUtilities::ib_5_delta_fcn,
                                                           6,
                                                           2.75,
                                                           { 0,
                                                             0.0233383448809079538,
                                                             0.333140121904304905,
                                                             0.514244366143986938,
                                                             0.128737522475479593,
                                                             0.000539644595320609716 } },
                                                         { PETScMatUtilities::ib_6_delta_fcn,
                                                           6,
                                                           2.25,
                                                           { 0.00965617417165844278,
                                                             0.174648694040214713,
                                                             0.431221688477088836,
                                                             0.325168575099164853,
                                                             0.0591221373512527211,
                                                             0.000182730860620434541 } },
                                                         { PETScMatUtilities::ib_6_delta_fcn,
                                                           6,
                                                           2.75,
                                                           { 0.000182730860620434541,
                                                             0.0591221373512527211,
                                                             0.325168575099164853,
                                                             0.431221688477088836,
                                                             0.174648694040214713,
                                                             0.00965617417165844278 } } } };
    for (const auto& sample : off_center_samples)
    {
        std::array<double, 8> weights;
        weights.fill(-123.0);
        sample.kernel(sample.r, weights.data());
        for (int i = 0; i < sample.width; ++i)
            kernels = kernels && std::isfinite(weights[i]) && std::abs(weights[i] - sample.expected[i]) < 1.0e-12;
        kernels = kernels && weights[6] == -123.0 && weights[7] == -123.0;
    }
    // Partition of unity and linear reproduction across the callback's full
    // lower-stencil displacement interval, including the IB_5 center switch.
    for (Kernel kernel : { PETScMatUtilities::ib_5_delta_fcn, PETScMatUtilities::ib_6_delta_fcn })
        for (int k = 0; k <= 64; ++k)
        {
            const double r = 2.0 + k / 64.0;
            std::array<double, 6> weights;
            kernel(r, weights.data());
            double sum = 0.0, first_moment = 0.0;
            for (unsigned int i = 0; i < weights.size(); ++i)
            {
                sum += weights[i];
                first_moment += i * weights[i];
            }
            kernels = kernels && std::isfinite(sum) && std::isfinite(first_moment) && std::abs(sum - 1.0) < 1.0e-12 &&
                      std::abs(first_moment - r) < 1.0e-12;
        }
    pout << "enum_aliases_valid = " << (enums ? "true" : "false") << '\n';
    pout << "kernel_values_valid = " << (kernels ? "true" : "false") << '\n';
    return !enums + !kernels;
}

void
generate_probes(const unsigned int& structure, const int& level, int& count, std::vector<IBTK::Point>& positions, void*)
{
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
expected_weight(int width, int offset, double distance)
{
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
check_matrix(Mat matrix, Vec positions, Pointer<SideData<NDIM, int>> dofs, int component_width, int transverse_width)
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
                expected[column] = expected_weight(width[0], i, x) * expected_weight(width[1], j, y);
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
    PetscErrorCode ierr;
    int failures = 0;
    {
        Pointer<AppInitializer> app = new AppInitializer(argc, argv, "interpolation.log");
        failures += check_enums_and_kernels();
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
        initializer->setStructureNamesOnLevel(0, { "probes" });
        initializer->registerInitStructureFunction(generate_probes);
        method->registerLInitStrategy(initializer);
        gridding->makeCoarsestLevel(hierarchy, 0.0);
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        // This fixture deliberately uses a single patch/rank; all tested stencils
        // lie in its interior, so neither ghost-numbering nor boundary extensions
        // can obscure the centering/column assertions.
        if (IBTK_MPI::getNodes() != 1 || level->getNumberOfPatches() != 1)
            TBOX_ERROR("Interpolation fixture requires one patch on one rank\n");
        VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = variables->getContext("interpolation");
        Pointer<SideVariable<NDIM, double>> velocity = new SideVariable<NDIM, double>("velocity");
        Pointer<SideVariable<NDIM, int>> indices = new SideVariable<NDIM, int>("indices");
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

        const std::array<Kernel, 4> kernel = { PETScMatUtilities::piecewise_constant_delta_fcn,
                                               PETScMatUtilities::piecewise_linear_delta_fcn,
                                               three_point_probe,
                                               PETScMatUtilities::ib_4_delta_fcn };
        bool placement = true, scalar_equivalence = true, lifecycle = true;
        // Repeat setup/use/cleanup on the same IBMethod to exercise scratch state
        // invalidation, with fixed-operator updates at both new and half times.
        for (int step = 0; step < 2; ++step)
        {
            const double current_time = step * 0.125, new_time = current_time + 0.125;
            method->preprocessIntegrateData(current_time, new_time, 1);
            method->updateFixedLEOperators();
            Vec X = method->getLDataManager()->getLData("X", 0)->getVec();
            for (int cw = 1; cw <= 4; ++cw)
                for (int tw = 1; tw <= 4; ++tw)
                {
                    Mat matrix = nullptr;
                    IBImplicitStrategy& strategy = *method;
                    strategy.constructInterpOp(matrix, kernel[cw - 1], cw, kernel[tw - 1], tw, counts, dof, new_time);
                    placement = check_matrix(matrix, X, dofs, cw, tw) && placement;
                    if (cw == tw)
                    {
                        Mat scalar = nullptr;
                        strategy.constructInterpOp(scalar, kernel[cw - 1], cw, counts, dof, new_time);
                        PetscBool equal;
                        ierr = MatEqual(matrix, scalar, &equal);
                        IBTK_CHKERRQ(ierr);
                        scalar_equivalence = scalar_equivalence && equal;
                        PETScMatUtilities::constructPatchLevelSCInterpOp(
                            scalar, kernel[cw - 1], cw, X, counts, dof, level);
                        ierr = MatEqual(matrix, scalar, &equal);
                        IBTK_CHKERRQ(ierr);
                        scalar_equivalence = scalar_equivalence && equal;
                        ierr = MatDestroy(&scalar);
                        IBTK_CHKERRQ(ierr);
                    }
                    ierr = MatDestroy(&matrix);
                    IBTK_CHKERRQ(ierr);
                }
            Mat half = nullptr;
            method->constructInterpOp(half, kernel[0], 1, counts, dof, current_time + 0.0625);
            lifecycle = check_matrix(half, X, dofs, 1, 1) && lifecycle;
            ierr = MatDestroy(&half);
            IBTK_CHKERRQ(ierr);
            method->postprocessIntegrateData(current_time, new_time, 1);
            method->postprocessData();
        }
        failures += !placement + !scalar_equivalence + !lifecycle;
        pout << "stencil_columns_and_weights_valid = " << (placement ? "true" : "false") << '\n';
        pout << "scalar_component_transverse_equivalence_valid = " << (scalar_equivalence ? "true" : "false") << '\n';
        pout << "standalone_lifecycle_valid = " << (lifecycle ? "true" : "false") << '\n';
        level->deallocatePatchData(u);
        level->deallocatePatchData(dof);
        variables->removePatchDataIndex(u);
        variables->removePatchDataIndex(dof);
        pout << "test_failures = " << failures << std::endl;
    }
    return failures;
}
