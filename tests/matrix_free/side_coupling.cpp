// Copyright (c) 2026 by the IBAMR developers
// This file is part of IBAMR and is distributed under the 3-clause BSD license.
#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/ib_kernels.h>

#include <CartesianPatchGeometry.h>
#include <SideIterator.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <numeric>

#include "coupling.h"
#include "fixture.h"

#include <ibtk/app_namespaces.h>

namespace
{
template <class Evaluator>
void
check_case(const std::string& name, const Evaluator& evaluator, const bool compare_fortran)
{
    using namespace MatrixFreeTest;
    const Pointer<Patch<NDIM>> patch = make_patch(8);
    const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
    const Box<NDIM>& box = patch->getBox();
    Pointer<SideData<NDIM, double>> field = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(4));
    Pointer<SideData<NDIM, double>> spread = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(4));
    Pointer<SideData<NDIM, double>> legacy = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(4));
    fill_field(*field);
    std::vector<double> positions = make_positions(*patch, 6, false);
    // Overlap two markers; include face/center ties and both patch edges.
    for (int d = 0; d < NDIM; ++d)
    {
        positions[NDIM + d] = positions[d];
        positions[2 * NDIM + d] = geometry->getXLower()[d] + (d % 2 ? 3.0 : 2.5) * geometry->getDx()[d];
        positions[3 * NDIM + d] = geometry->getXUpper()[d] - 0.125 * geometry->getDx()[d];
    }
    const std::vector<double> saved_positions = positions;
    std::vector<int> indices(6);
    std::iota(indices.begin(), indices.end(), 0);
    std::vector<double> force(positions.size()), values(positions.size()), reference_values(positions.size(), 0.0);
    for (std::size_t i = 0; i < force.size(); ++i)
    {
        force[i] = 0.3 + 0.1 * static_cast<double>(i);
    }
    const std::vector<double> saved_force = force;
    couple<false>(evaluator, *patch, *field, positions, indices, {}, values.data());
    spread->fillAll(0.25);
    couple<true>(evaluator, *patch, *spread, positions, indices, {}, force.data());
    double volume = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        volume *= geometry->getDx()[d];
    }
    double spread_error = 0.0, adjoint_grid = 0.0, spread_norm = 0.0;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (SideIterator<NDIM> i(field->getGhostBox(), axis); i; i++)
        {
            double expected = 0.0;
            for (int point : indices)
            {
                double weight = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    const double grid_x = geometry->getXLower()[d] +
                                          (i()(d) - box.lower()(d) + (d == axis ? 0.0 : 0.5)) * geometry->getDx()[d];
                    weight *=
                        reference_weight(name, axis, d, (positions[NDIM * point + d] - grid_x) / geometry->getDx()[d]);
                }
                reference_values[NDIM * point + axis] += weight * (*field)(i());
                expected += weight * force[NDIM * point + axis];
            }
            const double actual = ((*spread)(i()) - 0.25) * volume;
            spread_error = std::max(spread_error, std::abs(actual - expected));
            adjoint_grid += actual * (*field)(i());
            spread_norm += actual * actual;
        }
    }
    double gather_error = 0.0, adjoint_markers = 0.0;
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        gather_error = std::max(gather_error, std::abs(values[i] - reference_values[i]));
        adjoint_markers += values[i] * force[i];
    }
    double legacy_gather_error = 0.0, legacy_spread_error = 0.0;
    if (compare_fortran)
    {
        std::vector<double> legacy_values(values.size());
        LEInteractor::interpolate(legacy_values, NDIM, positions, NDIM, field, patch, box, name);
        legacy->fillAll(0.25);
        LEInteractor::spread(legacy, force, NDIM, positions, NDIM, patch, box, name);
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            legacy_gather_error = std::max(legacy_gather_error, std::abs(values[i] - legacy_values[i]));
        }
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (SideIterator<NDIM> i(field->getGhostBox(), axis); i; i++)
            {
                legacy_spread_error = std::max(legacy_spread_error, volume * std::abs((*spread)(i()) - (*legacy)(i())));
            }
        }
    }
    if (positions != saved_positions || force != saved_force)
    {
        TBOX_ERROR("Coupling modified its input markers.\n");
    }
    plog << name << " value " << values[2 * NDIM] << " spread_norm " << std::sqrt(spread_norm) << " reference_errors "
         << gather_error << ' ' << spread_error << " adjoint_error " << std::abs(adjoint_grid - adjoint_markers);
    if (compare_fortran)
    {
        plog << " Fortran_gather_error " << legacy_gather_error
             << (NDIM == 3 && name == "IB_5" ? " known_Fortran_spread_defect " : " Fortran_spread_error ")
             << legacy_spread_error;
    }
    else
    {
        plog << " Fortran_comparison unavailable";
    }
    plog << '\n';
}

void
check_indexed_clipping()
{
    using namespace MatrixFreeTest;
    const Pointer<Patch<NDIM>> patch = make_patch(8);
    const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
    SideData<NDIM, double> field(patch->getBox(), 1, IntVector<NDIM>(0));
    SideData<NDIM, double> spread(patch->getBox(), 1, IntVector<NDIM>(0));
    fill_field(field);
    std::vector<double> positions = make_positions(*patch, 4, false);
    const std::array<int, 3> indices = { 3, 1, 3 };
    std::array<double, 3 * NDIM> shifts{};
    for (int d = 0; d < NDIM; ++d)
    {
        shifts[d] = -0.73125 * geometry->getDx()[d];
        shifts[2 * NDIM + d] = 0.25 * geometry->getDx()[d];
    }
    const IBKernelEvaluatorTensorProduct kernel{ IBKernels::BSpline<3>{}, IBKernels::BSpline<2>{} };
    std::vector<double> force(positions.size(), 0.75), values(positions.size(), -17.0);
    // Duplicate indices with different shifts are meaningful for spread. Gather uses unique indices.
    couple<false, float>(kernel,
                         *patch,
                         field,
                         positions,
                         std::span(indices).first(2),
                         std::span(shifts).first(2 * NDIM),
                         values.data());
    spread.fillAll(0.0);
    couple<true>(kernel, *patch, spread, positions, indices, shifts, force.data());
    double volume = 1.0, gather_error = 0.0, spread_error = 0.0;
    for (int d = 0; d < NDIM; ++d)
    {
        volume *= geometry->getDx()[d];
    }
    for (int axis = 0; axis < NDIM; ++axis)
    {
        std::array<double, 2> expected_values{};
        for (SideIterator<NDIM> i(field.getGhostBox(), axis); i; i++)
        {
            double expected_spread = 0.0;
            for (std::size_t point = 0; point < indices.size(); ++point)
            {
                double weight = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    const double grid_x =
                        geometry->getXLower()[d] +
                        (i()(d) - patch->getBox().lower()(d) + (d == axis ? 0.0 : 0.5)) * geometry->getDx()[d];
                    weight *=
                        reference_weight("COMPOSITE_BSPLINE_32",
                                         axis,
                                         d,
                                         (positions[NDIM * indices[point] + d] + shifts[NDIM * point + d] - grid_x) /
                                             geometry->getDx()[d]);
                }
                expected_spread += weight * 0.75;
                if (point < 2)
                {
                    expected_values[point] += weight * field(i());
                }
            }
            spread_error = std::max(spread_error, std::abs(volume * spread(i()) - expected_spread));
        }
        for (std::size_t point = 0; point < 2; ++point)
        {
            gather_error =
                std::max(gather_error, std::abs(values[NDIM * indices[point] + axis] - expected_values[point]));
        }
        if (values[axis] != -17.0 || values[2 * NDIM + axis] != -17.0)
        {
            TBOX_ERROR("Gather changed an unselected marker.\n");
        }
    }
    plog << "indexed_clipping float_gather_error " << gather_error << " spread_error " << spread_error
         << " selected_value " << values[NDIM] << '\n';
}
} // namespace

int
main(int argc, char** argv)
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    PIO::logOnlyNodeZero("output");
    plog << std::setprecision(12) << std::scientific;
    check_case("IB_4", IBKernelEvaluatorTensorProduct{ IBKernels::IB4{} }, true);
    check_case("IB_5", IBKernelEvaluatorTensorProduct{ IBKernels::IB5{} }, true);
    check_case("BSPLINE_3", IBKernelEvaluatorTensorProduct{ IBKernels::BSpline<3>{} }, true);
    check_case("BSPLINE_6", IBKernelEvaluatorTensorProduct{ IBKernels::BSpline<6>{} }, true);
    check_case("COMPOSITE_BSPLINE_32",
               IBKernelEvaluatorTensorProduct{ IBKernels::BSpline<3>{}, IBKernels::BSpline<2>{} },
               true);
    check_case("COMPOSITE_BSPLINE_23",
               IBKernelEvaluatorTensorProduct{ IBKernels::BSpline<2>{}, IBKernels::BSpline<3>{} },
               true);
    check_case("COSINE_4", IBKernelEvaluatorTensorProduct{ MatrixFreeTest::CosineKernel{} }, false);
    check_indexed_clipping();
    return 0;
}
