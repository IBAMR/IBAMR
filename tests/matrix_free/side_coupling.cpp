// Copyright (c) 2026 by the IBAMR developers
// This file is part of IBAMR and is distributed under the 3-clause BSD license.
#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/IBOperator.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IndexUtilities.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/ib_kernel_evaluators.h>

#include <tbox/InputManager.h>
#include <tbox/MemoryDatabase.h>

#include <CartesianPatchGeometry.h>
#include <PatchDescriptor.h>
#include <SideIterator.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <tuple>
#include <vector>

#include "cartesian_coupling.h"
#include "coupling.h"
#include "fixture.h"

#include <ibtk/app_namespaces.h>

namespace
{
using Mode = IBTK::Experimental::TensorProductMode;

// Deliberately nonseparable, with nonzero support endpoints and const-only input.
class EndpointEvaluator
{
public:
    EndpointEvaluator() : d_coefficients(1 << NDIM)
    {
        for (std::size_t k = 0; k < d_coefficients.size(); ++k)
        {
            d_coefficients[k] = k + 1;
        }
    }
    template <int Axis>
    static constexpr std::array<std::size_t, NDIM> get_stencil_widths()
    {
        std::array<std::size_t, NDIM> widths;
        widths.fill(2);
        return widths;
    }
    template <int Axis, class Output, class Input>
    Output evaluate(const std::array<Input, NDIM>&) const
    {
        Output weights{};
        for (std::size_t k = 0; k < weights.size(); ++k)
        {
            weights[k] = d_coefficients[k];
        }
        return weights;
    }
    template <int Axis, class Output, class Input>
    Output evaluate(std::array<Input, NDIM>&) const = delete;

private:
    std::vector<double> d_coefficients;
};
struct CopyOnlyEndpoint : EndpointEvaluator
{
    CopyOnlyEndpoint() = default;
    CopyOnlyEndpoint(const CopyOnlyEndpoint&) = default;
    CopyOnlyEndpoint(CopyOnlyEndpoint&&) = delete;
};
struct MoveOnlyEndpoint : EndpointEvaluator
{
    MoveOnlyEndpoint() = default;
    MoveOnlyEndpoint(const MoveOnlyEndpoint&) = delete;
    MoveOnlyEndpoint(MoveOnlyEndpoint&&) = default;
};
struct ConstFactors : IBKernelEvaluatorTensorProduct<IBKernelEvaluators::BSpline<2>>
{
    ConstFactors() : IBKernelEvaluatorTensorProduct(IBKernelEvaluators::BSpline<2>{})
    {
    }
    using IBKernelEvaluatorTensorProduct::evaluateFactors;
    template <int Axis, class Coefficient, class Input>
    auto evaluateFactors(std::array<Input, NDIM>&) const = delete;
};

// A non-templated consumer retains and forwards a shared immutable handle.
class OperatorConsumer
{
public:
    explicit OperatorConsumer(IBOperator op) : d_operator(std::move(op))
    {
    }
    const IBOperator& getOperator() const
    {
        return d_operator;
    }

private:
    IBOperator d_operator;
};

void
check_operator_contracts()
{
    const OperatorConsumer consumer = []
    {
        CopyOnlyEndpoint evaluator;
        IBOperator op(evaluator);
        return OperatorConsumer(op);
    }();
    const IBOperator moved(MoveOnlyEndpoint{});
    const Pointer<Patch<NDIM>> patch = MatrixFreeTest::make_patch(8);
    const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
    SideData<NDIM, double> field(patch->getBox(), 1, IntVector<NDIM>(4));
    MatrixFreeTest::fill_field(field);
    const std::array<int, 1> indices{ 0 };
    double maximum_error = 0.0, tie_value = 0.0;
    for (int neighbor : { -1, 0, 1 })
    {
        std::array<double, NDIM> positions;
        hier::Index<NDIM> lower;
        for (int d = 0; d < NDIM; ++d)
        {
            const double tie = geometry->getXLower()[d] + (4.0 + (d == 0 ? 0.0 : 0.5)) * geometry->getDx()[d];
            positions[d] = neighbor == 0 ? tie : std::nextafter(tie, neighbor < 0 ? -INFINITY : INFINITY);
        }
        // Use SAMRAI's physical-cell mapping; direct/assembled parity is checked
        // independently with actual matrices in interpolation_utilities.
        const hier::Index<NDIM> cell = IndexUtilities::getCellIndex(positions, geometry, patch->getBox());
        for (int d = 0; d < NDIM; ++d)
        {
            const double center =
                (cell(d) - patch->getBox().lower()(d) + 0.5) * geometry->getDx()[d] + geometry->getXLower()[d];
            lower(d) = cell(d) - (d != 0 && positions[d] <= center ? 1 : 0);
        }
        double expected = 0.0;
        for (int n = 0; n < (1 << NDIM); ++n)
        {
            hier::Index<NDIM> index = lower;
            for (int d = 0; d < NDIM; ++d)
            {
                index(d) += (n >> d) & 1;
            }
            expected += (n + 1) * field(SideIndex<NDIM>(index, 0, SideIndex<NDIM>::Lower));
        }
        for (const IBOperator* op : { &consumer.getOperator(), &moved })
        {
            double actual = -1.0;
            op->interpolate<DataCentering::SIDE>(*patch, field, 0, 0, 0, positions, indices, {}, &actual);
            TBOX_ASSERT(std::isfinite(actual));
            maximum_error = std::max(maximum_error, std::abs(actual - expected));
            if (neighbor == 0)
            {
                tie_value = actual;
            }
        }
        const IBOperator factors{ ConstFactors{} };
        double actual = 0.0;
        factors.interpolate<DataCentering::SIDE>(*patch, field, 0, 0, 0, positions, indices, {}, &actual);
        factors.spread<DataCentering::SIDE>(*patch, field, 0, 0, 0, positions, {}, {}, nullptr);
    }
    TBOX_ASSERT(maximum_error < 1.0e-12);
    plog << "endpoint_ties max_error " << maximum_error << " tie_value " << tie_value << '\n';
}

// Check the real operator against physical hat-function weights, without using
// the shared stencil helpers or the assembled operator's ratio convention.
void
check_signed_ratios()
{
    const hier::Box<NDIM> base_box(hier::Index<NDIM>(0), hier::Index<NDIM>(23));
    hier::BoxArray<NDIM> domain(1);
    domain[0] = base_box;
    std::array<double, NDIM> x_lower, x_upper, position;
    x_lower.fill(0.3);
    x_upper.fill(2.7);
    position.fill(1.13);
    const Pointer<CartesianGridGeometry<NDIM>> grid =
        new CartesianGridGeometry<NDIM>("signed_ratio_geometry", x_lower.data(), x_upper.data(), domain, false);
    tbox::Array<tbox::Array<bool>> boundaries(NDIM);
    for (int d = 0; d < NDIM; ++d)
    {
        boundaries[d].resizeArray(2);
        boundaries[d][0] = boundaries[d][1] = true;
    }
    const IBOperator op{ IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::BSpline<2>{} } };
    const std::array<int, 1> indices{ 0 };
    const double force = 1.7;
    for (const int ratio : { -2, 3 })
    {
        const hier::Box<NDIM> box = ratio < 0 ? hier::Box<NDIM>::coarsen(base_box, IntVector<NDIM>(-ratio)) :
                                                hier::Box<NDIM>::refine(base_box, IntVector<NDIM>(ratio));
        const Pointer<Patch<NDIM>> patch = new Patch<NDIM>(box, new PatchDescriptor<NDIM>());
        tbox::Array<tbox::Array<bool>> periodic_boundaries(NDIM);
        for (int d = 0; d < NDIM; ++d)
        {
            periodic_boundaries[d].resizeArray(2);
            periodic_boundaries[d][0] = periodic_boundaries[d][1] = false;
        }
        grid->setGeometryDataOnPatch(*patch, IntVector<NDIM>(ratio), boundaries, periodic_boundaries);
        const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
        SideData<NDIM, double> field(box, 1, IntVector<NDIM>(0)), spread(box, 1, IntVector<NDIM>(0));
        const double dx = (x_upper[0] - x_lower[0]) / box.numberCells(0);
        const double volume = std::pow(dx, NDIM);
        if (!std::isfinite(volume) || volume <= 0.0)
        {
            TBOX_ERROR("Invalid volume in signed-ratio regression.\n");
        }
        double gather_error = 0.0, spread_error = 0.0, mass_error = 0.0, adjoint_error = 0.0;
        TBOX_ASSERT(geometry->getRatio() == IntVector<NDIM>(ratio));
        for (int d = 0; d < NDIM; ++d)
        {
            TBOX_ASSERT(geometry->getDx()[d] > 0.0);
            TBOX_ASSERT(std::abs(geometry->getDx()[d] - dx) < 1.0e-14);
        }
        for (int axis = 0; axis < NDIM; ++axis)
        {
            double expected = 0.0;
            for (SideIterator<NDIM> i(box, axis); i; i++)
            {
                double value = 1.0, weight = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    const double x = x_lower[d] + (i()(d) + (d == axis ? 0.0 : 0.5)) * dx;
                    value += (d + 1) * x * x;
                    weight *= std::max(0.0, 1.0 - std::abs(position[d] - x) / dx);
                }
                field(i()) = value;
                expected += weight * value;
            }
            double actual = 0.0;
            op.interpolate<DataCentering::SIDE>(*patch, field, axis, 0, axis, position, indices, {}, &actual, 1, grid);
            if (!std::isfinite(actual) || !std::isfinite(expected))
            {
                TBOX_ERROR("Nonfinite gather in signed-ratio regression.\n");
            }
            spread.fillAll(0.0);
            op.spread<DataCentering::SIDE>(*patch, spread, axis, 0, axis, position, indices, {}, &force, 1, grid);
            double mass = 0.0, pairing = 0.0;
            for (SideIterator<NDIM> i(box, axis); i; i++)
            {
                double weight = 1.0;
                for (int d = 0; d < NDIM; ++d)
                {
                    const double x = x_lower[d] + (i()(d) + (d == axis ? 0.0 : 0.5)) * dx;
                    weight *= std::max(0.0, 1.0 - std::abs(position[d] - x) / dx);
                }
                if (!std::isfinite(spread(i())) || !std::isfinite(weight) ||
                    !std::isfinite(volume * spread(i()) - force * weight))
                {
                    TBOX_ERROR("Nonfinite spread in signed-ratio regression.\n");
                }
                spread_error = std::max(spread_error, std::abs(volume * spread(i()) - force * weight));
                mass += volume * spread(i());
                pairing += volume * spread(i()) * field(i());
            }
            if (!std::isfinite(mass) || !std::isfinite(pairing) || !std::isfinite(actual - expected) ||
                !std::isfinite(mass - force) || !std::isfinite(pairing - force * actual))
            {
                TBOX_ERROR("Nonfinite reduction in signed-ratio regression.\n");
            }
            gather_error = std::max(gather_error, std::abs(actual - expected));
            mass_error = std::max(mass_error, std::abs(mass - force));
            adjoint_error = std::max(adjoint_error, std::abs(pairing - force * actual));
        }
        plog << "ratio " << ratio << " dx " << dx << " errors " << gather_error << ' ' << spread_error << ' '
             << mass_error << ' ' << adjoint_error << std::endl;
        TBOX_ASSERT(gather_error < 1.0e-11 && spread_error < 1.0e-11 && mass_error < 1.0e-11 &&
                    adjoint_error < 1.0e-11);
    }
}

template <class Evaluator>
void
check_case(const std::string& name, const Evaluator& evaluator, const bool fortran_available, const int ghosts)
{
    using namespace MatrixFreeTest;
    const Pointer<Patch<NDIM>> patch = make_patch(8);
    const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
    const Box<NDIM>& box = patch->getBox();
    Pointer<SideData<NDIM, double>> field = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(ghosts));
    Pointer<SideData<NDIM, double>> spread = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(ghosts));
    Pointer<SideData<NDIM, double>> legacy = new SideData<NDIM, double>(box, 1, IntVector<NDIM>(ghosts));
    const bool compare_fortran = fortran_available && ghosts > 0;
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
    const auto compare_mode = [&]<Mode Application>()
    {
        SideData<NDIM, double> other_spread(box, 1, IntVector<NDIM>(ghosts));
        other_spread.fillAll(0.25);
        std::vector<double> other_values(values.size());
        couple<false, double, Application>(evaluator, *patch, *field, positions, indices, {}, other_values.data());
        couple<true, double, Application>(evaluator, *patch, other_spread, positions, indices, {}, force.data());
        std::array<double, 3> errors{};
        double grid_pairing = 0.0, marker_pairing = 0.0;
        for (std::size_t i = 0; i < values.size(); ++i)
        {
            errors[0] = std::max(errors[0], std::abs(other_values[i] - values[i]));
            marker_pairing += other_values[i] * force[i];
        }
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (SideIterator<NDIM> i(field->getGhostBox(), axis); i; i++)
            {
                errors[1] = std::max(errors[1], volume * std::abs(other_spread(i()) - (*spread)(i())));
                grid_pairing += volume * (other_spread(i()) - 0.25) * (*field)(i());
            }
        }
        errors[2] = std::abs(grid_pairing - marker_pairing);
        if (*std::max_element(errors.begin(), errors.end()) > 1.0e-11)
        {
            TBOX_ERROR("Tensor-product application mismatch for " << name << '\n');
        }
        return errors;
    };
    const std::array<double, 3> expanded_errors = compare_mode.template operator()<Mode::EXPANDED>();
    const std::array<double, 3> factorized_errors = compare_mode.template operator()<Mode::FACTORIZED>();
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
    plog << name;
    if (name.starts_with("COMPOSITE_BSPLINE_"))
    {
        plog << " normal " << name[name.size() - 2] << " tangential " << name.back();
    }
    plog << " ghosts " << ghosts << " value " << values[2 * NDIM] << " spread_norm " << std::sqrt(spread_norm)
         << " reference_errors " << gather_error << ' ' << spread_error << " adjoint_error "
         << std::abs(adjoint_grid - adjoint_markers);
    if (compare_fortran)
    {
        plog << " Fortran_gather_error " << legacy_gather_error << " Fortran_spread_error " << legacy_spread_error;
    }
    else
    {
        plog << (fortran_available ? " Fortran_comparison not_run_clipped" : " Fortran_comparison unavailable");
    }
    plog << " expanded_difference " << expanded_errors[0] << ' ' << expanded_errors[1] << " factorized_difference "
         << factorized_errors[0] << ' ' << factorized_errors[1] << " other_adjoint_error "
         << std::max(expanded_errors[2], factorized_errors[2]) << '\n';
}

template <Mode Application>
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
    const IBKernelEvaluatorTensorProduct kernel{ IBKernelEvaluators::BSpline<3>{}, IBKernelEvaluators::BSpline<2>{} };
    std::vector<double> force(positions.size(), 0.75), values(positions.size(), -17.0);
    // Duplicate indices with different shifts are meaningful for spread. Gather uses unique indices.
    couple<false, float, Application>(kernel,
                                      *patch,
                                      field,
                                      positions,
                                      std::span(indices).first(2),
                                      std::span(shifts).first(2 * NDIM),
                                      values.data());
    spread.fillAll(0.0);
    couple<true, float, Application>(kernel, *patch, spread, positions, indices, shifts, force.data());
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
    plog << "indexed_clipping "
         << (Application == Mode::EXPANDED   ? "expanded" :
             Application == Mode::FACTORIZED ? "factorized" :
                                               "contracted")
         << " float_gather_error " << gather_error << " spread_error " << spread_error << " selected_value "
         << values[NDIM] << '\n';
}
void
check_factor_ownership()
{
    const auto factors = []
    {
        std::array<double, NDIM> r;
        r.fill(2.125);
        return IBKernelEvaluatorTensorProduct{
            IBKernelEvaluators::BSpline<6>{}, IBKernelEvaluators::BSpline<5>{}
        }.evaluateFactors<0>(r);
    }();
    auto copied_factors = factors;
    std::get<0>(copied_factors)[0] += 1.0;
    std::array<double, NDIM> r;
    r.fill(2.125);
    const IBKernelEvaluatorTensorProduct kernel{ IBKernelEvaluators::BSpline<6>{}, IBKernelEvaluators::BSpline<5>{} };
    const auto expanded =
        kernel.evaluate<0, IBKernelEvaluators::Weights<double, detail::ib_kernel_stencil_size<decltype(kernel), 0>()>>(
            r);
    double error = 0.0;
    for (std::size_t n = 0; n < expanded.size(); ++n)
    {
        double weight = std::get<0>(factors)[n % 6] * std::get<1>(factors)[(n / 6) % 5];
#if (NDIM == 3)
        weight *= std::get<2>(factors)[n / 30];
#endif
        error = std::max(error, std::abs(weight - expanded[n]));
    }
    plog << "owned_factors product_error " << error << " independent_copy_delta "
         << std::get<0>(copied_factors)[0] - std::get<0>(factors)[0] << '\n';
}
} // namespace

int
main(int argc, char** argv)
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    PIO::logOnlyNodeZero("output");
    plog << std::setprecision(12) << std::scientific;
    if (argc > 1)
    {
        Pointer<MemoryDatabase> input = new MemoryDatabase("input");
        InputManager::getManager()->parseInputFile(argv[1], input);
        if (input->getBoolWithDefault("test_signed_ratios", false))
        {
            check_signed_ratios();
            return 0;
        }
        if (input->getBoolWithDefault("test_operator_contracts", false))
        {
            check_operator_contracts();
            return 0;
        }
        if (input->getBoolWithDefault("test_centerings", false))
        {
            MatrixFreeTest::check_cartesian_coupling();
            return 0;
        }
    }
    for (int ghosts : { 4, 0 })
    {
        check_case("IB_4", IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::IB4{} }, true, ghosts);
        check_case("IB_5", IBKernelEvaluatorTensorProduct{ IBKernelEvaluators::IB5{} }, true, ghosts);
        MatrixFreeTest::for_each_bspline([&](const std::string& name, const auto& evaluator)
                                         { check_case(name, evaluator, name != "COMPOSITE_BSPLINE_12", ghosts); });
        check_case("COSINE_4", MatrixFreeTest::CartesianCosineKernel{}, false, ghosts);
    }
    check_indexed_clipping<Mode::EXPANDED>();
    check_indexed_clipping<Mode::FACTORIZED>();
    check_indexed_clipping<Mode::CONTRACTED>();
    check_factor_ownership();
    return 0;
}
