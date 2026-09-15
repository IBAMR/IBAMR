// Copyright (c) 2026 by the IBAMR developers
// This file is part of IBAMR and is distributed under the 3-clause BSD license.
#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/ib_kernels.h>

#include <CartesianPatchGeometry.h>
#include <SideGeometry.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>

#include "coupling.h"
#include "fixture.h"

#include <ibtk/app_namespaces.h>

namespace
{
template <class Evaluator>
void
benchmark(const std::string& name,
          const Evaluator& evaluator,
          const int cells,
          const int markers,
          const int iterations,
          const int repeats,
          const bool shuffled)
{
    using namespace MatrixFreeTest;
    const Pointer<Patch<NDIM>> patch = make_patch(cells);
    const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
    const Box<NDIM>& box = patch->getBox();
    Pointer<SideData<NDIM, double>> field = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(4));
    Pointer<SideData<NDIM, double>> spread = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(4));
    fill_field(*field);
    const std::vector<double> positions = make_positions(*patch, markers, shuffled);
    const std::vector<double> shifts(positions.size(), 0.0);
    std::vector<int> indices(markers);
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<double> values(positions.size()), force(positions.size());
    std::array<std::vector<double>, NDIM> component_values, component_force;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        component_values[axis].resize(markers);
        component_force[axis].resize(markers);
        for (int point = 0; point < markers; ++point)
        {
            force[NDIM * point + axis] = component_force[axis][point] =
                0.5 + 0.03125 * ((static_cast<int>((positions[NDIM * point + axis] - geometry->getXLower()[axis]) /
                                                   geometry->getDx()[axis]) +
                                  axis) %
                                 13);
        }
    }
    const Experimental::SideCoupling coupling(*patch, *field);
    FortranInterpolate* const fortran_gather = get_fortran_interpolate(name);
    FortranSpread* const fortran_scatter = get_fortran_spread(name);
    const IBKernelTensorProduct legacy_kernel(name);
    std::array<std::array<double, NDIM>, NDIM> x_lower, x_upper;
    std::array<Box<NDIM>, NDIM> side_boxes;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        side_boxes[axis] = SideGeometry<NDIM>::toSideBox(box, axis);
        for (int d = 0; d < NDIM; ++d)
        {
            x_lower[axis][d] = geometry->getXLower()[d] - (axis == d ? 0.5 * geometry->getDx()[d] : 0.0);
            x_upper[axis][d] = geometry->getXUpper()[d] + (axis == d ? 0.5 * geometry->getDx()[d] : 0.0);
        }
    }
    const IntVector<NDIM>& ghosts = field->getGhostCellWidth();
    const std::array<const char*, 8> labels = { "cpp_loop_gather",     "fortran_loop_gather", "cpp_patch_gather",
                                                "LEInteractor_gather", "cpp_loop_spread",     "fortran_loop_spread",
                                                "cpp_patch_spread",    "LEInteractor_spread" };
    const std::array<std::function<void()>, 8> operations = {
        [&]
        {
            [&]<std::size_t... Axis>(std::index_sequence<Axis...>)
            {
                (coupling.template interpolateAxis<Axis>(
                     evaluator, field->getPointer(Axis), positions, indices, shifts, component_values[Axis].data()),
                 ...);
            }(std::make_index_sequence<NDIM>{});
        },
        [&]
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM>& side = side_boxes[axis];
                fortran_gather(geometry->getDx(),
                               x_lower[axis].data(),
                               x_upper[axis].data(),
                               1,
                               side.lower()(0),
                               side.upper()(0),
                               side.lower()(1),
                               side.upper()(1),
#if (NDIM == 3)
                               side.lower()(2),
                               side.upper()(2),
#endif
                               ghosts(0),
                               ghosts(1),
#if (NDIM == 3)
                               ghosts(2),
#endif
                               field->getPointer(axis),
                               indices.data(),
                               shifts.data(),
                               markers,
                               positions.data(),
                               component_values[axis].data());
            }
        },
        [&] { couple<false>(evaluator, *patch, *field, positions, indices, {}, values.data()); },
        [&] { LEInteractor::interpolate(values, NDIM, positions, NDIM, field, patch, box, legacy_kernel); },
        [&]
        {
            [&]<std::size_t... Axis>(std::index_sequence<Axis...>)
            {
                (coupling.template spreadAxis<Axis>(
                     evaluator, spread->getPointer(Axis), positions, indices, shifts, component_force[Axis].data()),
                 ...);
            }(std::make_index_sequence<NDIM>{});
        },
        [&]
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM>& side = side_boxes[axis];
                fortran_scatter(geometry->getDx(),
                                x_lower[axis].data(),
                                x_upper[axis].data(),
                                1,
                                indices.data(),
                                shifts.data(),
                                markers,
                                positions.data(),
                                component_force[axis].data(),
                                side.lower()(0),
                                side.upper()(0),
                                side.lower()(1),
                                side.upper()(1),
#if (NDIM == 3)
                                side.lower()(2),
                                side.upper()(2),
#endif
                                ghosts(0),
                                ghosts(1),
#if (NDIM == 3)
                                ghosts(2),
#endif
                                spread->getPointer(axis));
            }
        },
        [&] { couple<true>(evaluator, *patch, *spread, positions, indices, {}, force.data()); },
        [&] { LEInteractor::spread(spread, force, NDIM, positions, NDIM, patch, box, legacy_kernel); }
    };
    const auto consume = [&](const int operation)
    {
        double checksum = 0.0;
        if (operation < 4)
        {
            for (int point = 0; point < markers; ++point)
            {
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    checksum += (1.0 + 0.03125 * (point % 7)) *
                                (operation < 2 ? component_values[axis][point] : values[NDIM * point + axis]);
                }
            }
        }
        else
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const double* const data = spread->getPointer(axis);
                const int size = spread->getArrayData(axis).getBox().size();
                for (int i = 0; i < size; ++i)
                {
                    checksum += data[i] * (1.0 + 0.03125 * (i % 17));
                }
            }
        }
        return checksum;
    };
    const auto capture = [&](const int operation)
    {
        std::vector<double> result;
        if (operation < 4)
        {
            result.resize(positions.size());
            for (int point = 0; point < markers; ++point)
            {
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    result[NDIM * point + axis] =
                        operation < 2 ? component_values[axis][point] : values[NDIM * point + axis];
                }
            }
        }
        else
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const double* const data = spread->getPointer(axis);
                result.insert(result.end(), data, data + spread->getArrayData(axis).getBox().size());
            }
        }
        return result;
    };
    // Validate full outputs before timing, including the direct Fortran ABI calls.
    for (int operation = 0; operation < 8; ++operation)
    {
        spread->fillAll(0.0);
        operations[operation]();
        const std::vector<double> actual = capture(operation);
        spread->fillAll(0.0);
        operations[operation < 4 ? 0 : 4]();
        const std::vector<double> reference = capture(operation < 4 ? 0 : 4);
        for (std::size_t i = 0; i < actual.size(); ++i)
        {
            if (!std::isfinite(actual[i]) ||
                std::abs(actual[i] - reference[i]) > 1.0e-10 * std::max(1.0, std::abs(reference[i])))
            {
                TBOX_ERROR("Benchmark output mismatch for " << name << ' ' << labels[operation] << '\n');
            }
        }
    }
    for (int warmup = 0; warmup < 3; ++warmup)
    {
        for (const std::function<void()>& operation : operations)
        {
            spread->fillAll(0.0);
            operation();
        }
    }
    // Rotate implementation order across repeats; reset and checksums are outside timing.
    for (int repeat = 0; repeat < repeats; ++repeat)
    {
        for (int slot = 0; slot < 8; ++slot)
        {
            const int operation = (slot + repeat) % 8;
            spread->fillAll(0.0);
            const std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            for (int iteration = 0; iteration < iterations; ++iteration)
            {
                operations[operation]();
            }
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            std::cout << NDIM << ',' << name << ',' << cells << ',' << markers << ',' << shuffled << ',' << iterations
                      << ',' << repeat << ',' << labels[operation] << ',' << seconds << ',' << consume(operation)
                      << '\n';
        }
    }
}
} // namespace

int
main(int argc, char** argv)
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    if (argc != 6)
    {
        std::cerr << "Usage: benchmark cells markers iterations repeats shuffled(0|1)\n";
        return 1;
    }
    const int cells = std::atoi(argv[1]), markers = std::atoi(argv[2]), iterations = std::atoi(argv[3]),
              repeats = std::atoi(argv[4]);
    const bool shuffled = std::atoi(argv[5]);
    if (cells < 8 || markers < 1 || iterations < 1 || repeats < 3)
    {
        return 1;
    }
    std::cout << std::setprecision(17)
              << "dimension,kernel,cells,markers,shuffled,iterations,repeat,operation,seconds,checksum\n";
    benchmark(
        "IB_4", IBKernelEvaluatorTensorProduct{ IBKernels::IB4{} }, cells, markers, iterations, repeats, shuffled);
    benchmark("BSPLINE_3",
              IBKernelEvaluatorTensorProduct{ IBKernels::BSpline<3>{} },
              cells,
              markers,
              iterations,
              repeats,
              shuffled);
    benchmark("BSPLINE_6",
              IBKernelEvaluatorTensorProduct{ IBKernels::BSpline<6>{} },
              cells,
              markers,
              iterations,
              repeats,
              shuffled);
    return 0;
}
