// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/IBKernel.h>
#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/IBKernelTensorProduct.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/ib_kernel_evaluators.h>

#include <tbox/Utilities.h>

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <locale>
#include <set>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "../tests.h"

using IBTK::IBKernel;
using IBTK::IBKernelTensorProduct;

bool ib_kernel_static_initialization_valid();

struct MissingExtentWeights
{
    float operator[](std::size_t) const;
};

struct RuntimeExtentWeights
{
    float operator[](std::size_t) const;
};

struct FractionalExtentWeights
{
    float operator[](std::size_t) const;
};

template <>
struct IBTK::IBKernelWeightsTraits<MissingExtentWeights>
{
    using value_type = float;
};

template <>
struct IBTK::IBKernelWeightsTraits<RuntimeExtentWeights>
{
    using value_type = float;
    static std::size_t extent;
};

template <>
struct IBTK::IBKernelWeightsTraits<FractionalExtentWeights>
{
    using value_type = float;
    static constexpr double extent = 2.5;
};

static_assert(!IBTK::IBKernelWeights<MissingExtentWeights>);
static_assert(!IBTK::IBKernelWeights<RuntimeExtentWeights>);
static_assert(!IBTK::IBKernelWeights<FractionalExtentWeights>);

namespace
{
constexpr std::size_t ENCODED_BLOCK_COUNT = 2;

struct ScalarWidth
{
    static constexpr std::size_t get_stencil_width()
    {
        return 1;
    }
};

struct NonConstScalar : ScalarWidth
{
    template <class Output, class Input>
    Output evaluate(const Input&);
};

struct BorrowedScalar : ScalarWidth
{
    template <class Output, class Input>
    const Output& evaluate(const Input&) const;
};

struct EmptyScalar : ScalarWidth
{
    static constexpr std::size_t get_stencil_width()
    {
        return 0;
    }
};

struct ImmovableScalar : ScalarWidth
{
    ImmovableScalar() = default;
    ImmovableScalar(const ImmovableScalar&) = delete;
    ImmovableScalar(ImmovableScalar&&) = delete;
    template <class Output, class Input>
    Output evaluate(const Input&) const;
};

struct ExplicitConstructionScalar : ScalarWidth
{
    explicit ExplicitConstructionScalar(double value) : value(value)
    {
    }
    explicit ExplicitConstructionScalar(const ExplicitConstructionScalar&) = default;
    explicit ExplicitConstructionScalar(ExplicitConstructionScalar&&) = default;
    template <class Output, class Input>
    Output evaluate(const Input&) const
    {
        return Output{ value };
    }
    double value;
};

struct RvalueOnlyScalar : ScalarWidth
{
    template <class Output>
    Output evaluate(double&&) const;
};

struct FloatOnlyScalar : ScalarWidth
{
    template <class Output, class Input>
    requires std::same_as<IBTK::ib_kernel_weights_value_t<Output>, float> Output evaluate(const Input&) const
    {
        return Output{ 0.5f };
    }
};

struct MoveOnlyScalar : ScalarWidth
{
    MoveOnlyScalar() = default;
    MoveOnlyScalar(const MoveOnlyScalar&) = delete;
    MoveOnlyScalar(MoveOnlyScalar&&) = default;
    template <class Output, class Input>
    Output evaluate(const Input&) const;
};

struct NonConstantWidth
{
    static std::size_t get_stencil_width()
    {
        return 1;
    }
};

struct FractionalWidth
{
    static constexpr double get_stencil_width()
    {
        return 2.5;
    }
};

struct IntWidth
{
    static constexpr int get_stencil_width()
    {
        return 2;
    }
};

struct WrongTypeWidths
{
    template <int Axis>
    static constexpr std::array<int, NDIM> get_stencil_widths()
    {
        return {};
    }
    template <int Axis, class Output, class Input>
    Output evaluate(const std::array<Input, NDIM>&) const;
};

struct NonConstantWidths
{
    template <int Axis>
    static std::array<std::size_t, NDIM> get_stencil_widths()
    {
        std::array<std::size_t, NDIM> widths;
        widths.fill(1);
        return widths;
    }
    template <int Axis, class Output, class Input>
    Output evaluate(const std::array<Input, NDIM>&) const;
};

template <class NormalEvaluator, class TransverseEvaluator>
concept HasTensorProduct = requires
{
    typename IBTK::IBKernelEvaluatorTensorProduct<NormalEvaluator, TransverseEvaluator>;
};

template <std::size_t Width, std::size_t Count>
struct TensorShape
{
    template <int Axis>
    static constexpr std::array<std::size_t, NDIM> get_stencil_widths()
    {
        std::array<std::size_t, NDIM> widths;
        widths.fill(Width);
        return widths;
    }
    template <int Axis, class Output, class Input>
    requires(IBTK::IBKernelWeightsTraits<Output>::extent == Count) Output
        evaluate(const std::array<Input, NDIM>&) const;
};

static_assert(IBTK::IBKernelEvaluatorScalar<IBTK::IBKernelEvaluators::IB3>);
static_assert(IBTK::IBKernelEvaluatorScalar<IBTK::IBKernelEvaluators::IB4, long double, float>);
static_assert(IBTK::IBKernelEvaluatorScalar<IBTK::IBKernelEvaluators::IB5, float, long double>);
static_assert(IBTK::IBKernelEvaluatorScalar<IBTK::IBKernelEvaluators::IB6>);
static_assert(IBTK::IBKernelEvaluatorScalar<IBTK::IBKernelEvaluators::BSpline<9>>);
static_assert(!IBTK::IBKernelEvaluatorScalar<NonConstScalar>);
static_assert(!IBTK::IBKernelEvaluatorScalar<BorrowedScalar>);
static_assert(!IBTK::IBKernelEvaluatorScalar<EmptyScalar>);
static_assert(!IBTK::IBKernelEvaluatorScalar<int>);
static_assert(!IBTK::IBKernelEvaluatorScalar<RvalueOnlyScalar>);
static_assert(IBTK::IBKernelEvaluatorScalar<ImmovableScalar>);
static_assert(
    !std::is_constructible_v<IBTK::IBKernelEvaluatorTensorProduct<ImmovableScalar>, ImmovableScalar, ImmovableScalar>);
static_assert(requires {
    IBTK::IBKernelEvaluatorTensorProduct{ ExplicitConstructionScalar{ 0.5 } };
    IBTK::IBKernelEvaluatorTensorProduct{ ExplicitConstructionScalar{ 0.5 }, ExplicitConstructionScalar{ 0.25 } };
});
static_assert(!IBTK::IBKernelScalarStencil<NonConstantWidth>);
static_assert(!IBTK::IBKernelScalarStencil<FractionalWidth>);
static_assert(!IBTK::IBKernelScalarStencil<IntWidth>);
static_assert(!IBTK::IBKernelScalarStencil<int>);
static_assert(!IBTK::IBKernelEvaluatorCartesian<WrongTypeWidths>);
static_assert(!IBTK::IBKernelEvaluatorCartesian<NonConstantWidths>);
template <class Product, int Axis>
concept CanEvaluateAlongAxis = requires(const Product& product, const std::array<double, NDIM>& r)
{
    product.template evaluate<
        Axis,
        IBTK::IBKernelEvaluators::Weights<double, IBTK::detail::ib_kernel_stencil_size<Product, 0>()>>(r);
};
static_assert(CanEvaluateAlongAxis<IBTK::IBKernelEvaluatorTensorProduct<IBTK::IBKernelEvaluators::IB4>, 0>);
static_assert(!CanEvaluateAlongAxis<IBTK::IBKernelEvaluatorTensorProduct<IBTK::IBKernelEvaluators::IB4>, -1>);
static_assert(!CanEvaluateAlongAxis<IBTK::IBKernelEvaluatorTensorProduct<IBTK::IBKernelEvaluators::IB4>, NDIM>);
static_assert(IBTK::IBKernelEvaluatorScalar<FloatOnlyScalar, float, float>);
static_assert(!IBTK::IBKernelEvaluatorScalar<FloatOnlyScalar>);
static_assert(HasTensorProduct<FloatOnlyScalar, FloatOnlyScalar>);
static_assert(IBTK::IBKernelEvaluatorCartesian<IBTK::IBKernelEvaluatorTensorProduct<FloatOnlyScalar>, double, float>);
static_assert(!IBTK::IBKernelEvaluatorCartesian<IBTK::IBKernelEvaluatorTensorProduct<FloatOnlyScalar>, double, double>);
static_assert(!std::is_constructible_v<IBTK::IBKernelEvaluatorTensorProduct<MoveOnlyScalar>, MoveOnlyScalar>);
static_assert(
    std::is_constructible_v<IBTK::IBKernelEvaluatorTensorProduct<MoveOnlyScalar>, MoveOnlyScalar, MoveOnlyScalar>);
static_assert(!HasTensorProduct<int, IBTK::IBKernelEvaluators::IB4>);
static_assert(!HasTensorProduct<IBTK::IBKernelEvaluators::IB4&, IBTK::IBKernelEvaluators::IB4&>);
static_assert(IBTK::IBKernelEvaluatorCartesian<IBTK::IBKernelEvaluatorTensorProduct<IBTK::IBKernelEvaluators::IB4>>);
static_assert(IBTK::IBKernelEvaluatorCartesian<TensorShape<1, 1>>);
static_assert(!IBTK::IBKernelEvaluatorCartesian<TensorShape<0, 1>>);
static_assert(!IBTK::IBKernelEvaluatorCartesian<TensorShape<2, 1>>);
static_assert(!IBTK::IBKernelEvaluatorCartesian<int>);

template <class Kernel, class Output>
concept EvaluatesInto = requires(const Kernel& kernel, const double& r)
{
    {
        kernel.template evaluate<Output>(r)
    } -> std::same_as<Output>;
};
static_assert(!EvaluatesInto<IBTK::IBKernelEvaluators::IB4, IBTK::IBKernelEvaluators::Weights<double, 3>>);

std::vector<std::string>
spellings(const std::string& name)
{
    std::string lower = name, mixed = name;
    std::use_facet<std::ctype<char>>(std::locale::classic()).tolower(lower.data(), lower.data() + lower.size());
    for (std::size_t i = 0; i < name.size(); i += 2)
    {
        mixed[i] = lower[i];
    }
    return { name, lower, mixed };
}

// Test-only inspection of the compact layout, without adding a public key API.
std::array<std::uint64_t, ENCODED_BLOCK_COUNT>
numeric_blocks(const IBKernel& kernel)
{
    std::array<std::uint64_t, ENCODED_BLOCK_COUNT> blocks{};
    static_assert(sizeof(kernel) == sizeof(blocks), "identity must contain the expected number of blocks");
    std::memcpy(blocks.data(), &kernel, sizeof(kernel));
    return blocks;
}

IBKernelTensorProduct
copy_product(IBKernelTensorProduct kernel)
{
    return kernel;
}

template <class Evaluator, std::size_t N>
double
sample_error(const Evaluator& evaluator, double r, const IBTK::IBKernelEvaluators::Weights<double, N>& expected)
{
    const IBTK::IBKernelEvaluators::Weights<double, N> weights =
        evaluator.template evaluate<IBTK::IBKernelEvaluators::Weights<double, N>>(r);
    const IBTK::IBKernelEvaluators::Weights<float, N> float_weights =
        evaluator.template evaluate<IBTK::IBKernelEvaluators::Weights<float, N>>(r);
    const IBTK::IBKernelEvaluators::Weights<long double, N> extended_weights =
        evaluator.template evaluate<IBTK::IBKernelEvaluators::Weights<long double, N>>(r);
    static_assert(std::tuple_size<decltype(weights)>::value == N, "Natural stencil size changed");
    double error = 0.0;
    for (std::size_t i = 0; i < N; ++i)
    {
        const double entry_error = std::abs(weights[i] - expected[i]);
        TBOX_ASSERT(std::abs(float_weights[i] - expected[i]) <= 128 * std::numeric_limits<float>::epsilon());
        TBOX_ASSERT(std::abs(extended_weights[i] - expected[i]) <= 1.0e-12L);
        if (!(entry_error <= 1.0e-12))
        {
            TBOX_ERROR("Kernel sample error = " << entry_error << '\n');
        }
        error = std::max(error, entry_error);
    }
    return error;
}

template <class Evaluator>
double
moment_error(const Evaluator& evaluator)
{
    constexpr int width = Evaluator::get_stencil_width();
    double error = 0.0;
    for (int k = 0; k < 64; ++k)
    {
        const double r = 0.5 * width - 1.0 + k / 64.0;
        const IBTK::IBKernelEvaluators::Weights<double, width> w =
            evaluator.template evaluate<IBTK::IBKernelEvaluators::Weights<double, width>>(r);
        double sum = 0.0, moment = 0.0;
        for (int i = 0; i < width; ++i)
        {
            sum += w[i];
            moment += i * w[i];
        }
        const double sum_error = std::abs(sum - 1.0), first_moment_error = std::abs(moment - r);
        if (!(sum_error <= 1.0e-12 && first_moment_error <= 1.0e-12))
        {
            TBOX_ERROR("Kernel moment errors = " << sum_error << ", " << first_moment_error << '\n');
        }
        error = std::max({ error, sum_error, first_moment_error });
    }
    return error;
}

int
check_kernels()
{
    using namespace IBTK;
    const long double displacement = 0.8L;
    const float coefficient_displacement = displacement;
    const IBKernelEvaluators::Weights<float, 2> selected =
        IBKernelEvaluators::BSpline<2>{}.evaluate<IBKernelEvaluators::Weights<float, 2>>(displacement);
    TBOX_ASSERT(selected[0] == 1 - coefficient_displacement && selected[1] == coefficient_displacement);
    const double a = (2.0 - std::sqrt(2.0)) / 8.0, b = (2.0 + std::sqrt(2.0)) / 8.0;
    const double K6 = (59.0 - std::sqrt(261.0)) / 60.0;
    double error = std::max(
        { sample_error(IBKernelEvaluators::BSpline<1>{}, -0.25, IBKernelEvaluators::Weights<double, 1>{ 1.0 }),
          sample_error(IBKernelEvaluators::BSpline<2>{}, 0.25, IBKernelEvaluators::Weights<double, 2>{ 0.75, 0.25 }),
          sample_error(
              IBKernelEvaluators::BSpline<3>{}, 1.0, IBKernelEvaluators::Weights<double, 3>{ 0.125, 0.75, 0.125 }),
          sample_error(IBKernelEvaluators::BSpline<4>{},
                       1.5,
                       IBKernelEvaluators::Weights<double, 4>{ 1.0 / 48, 23.0 / 48, 23.0 / 48, 1.0 / 48 }),
          sample_error(IBKernelEvaluators::BSpline<5>{},
                       1.5,
                       IBKernelEvaluators::Weights<double, 5>{ 1.0 / 24, 11.0 / 24, 11.0 / 24, 1.0 / 24, 0 }),
          sample_error(IBKernelEvaluators::BSpline<6>{},
                       2.5,
                       IBKernelEvaluators::Weights<double, 6>{
                           1.0 / 3840, 237.0 / 3840, 1682.0 / 3840, 1682.0 / 3840, 237.0 / 3840, 1.0 / 3840 }),
          sample_error(
              IBKernelEvaluators::IB3{}, 1.0, IBKernelEvaluators::Weights<double, 3>{ 1.0 / 6, 2.0 / 3, 1.0 / 6 }),
          sample_error(IBKernelEvaluators::IB4{}, 1.5, IBKernelEvaluators::Weights<double, 4>{ a, b, b, a }),
          sample_error(IBKernelEvaluators::IB5{},
                       1.5,
                       IBKernelEvaluators::Weights<double, 5>{ 0.0612224005711746881,
                                                               0.438777599428825312,
                                                               0.438777599428825312,
                                                               0.0612224005711746881,
                                                               0 }),
          sample_error(IBKernelEvaluators::IB6{},
                       3.0,
                       IBKernelEvaluators::Weights<double, 6>{
                           0, -1.0 / 16 + K6 / 8, 0.25, 5.0 / 8 - K6 / 4, 0.25, -1.0 / 16 + K6 / 8 }) });
    // Exact rational values from the truncated-power definition at r = 11/4.
    error = std::max(error,
                     sample_error(IBKernelEvaluators::BSpline<7>{},
                                  2.75,
                                  IBKernelEvaluators::Weights<double, 7>{ 729.0 / 2949120,
                                                                          112546.0 / 2949120,
                                                                          963327.0 / 2949120,
                                                                          1434812.0 / 2949120,
                                                                          422087.0 / 2949120,
                                                                          15618.0 / 2949120,
                                                                          1.0 / 2949120 }));
    // Exact rational values from the truncated-power definition at r = 7/2.
    error = std::max(error,
                     sample_error(IBKernelEvaluators::BSpline<8>{},
                                  3.5,
                                  IBKernelEvaluators::Weights<double, 8>{ 1.0 / 645120,
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
                       sample_error(IBKernelEvaluators::IB5{},
                                    2.25,
                                    IBKernelEvaluators::Weights<double, 5>{ 0.000539644595320609716,
                                                                            0.128737522475479593,
                                                                            0.514244366143986938,
                                                                            0.333140121904304905,
                                                                            0.0233383448809079538 }),
                       sample_error(IBKernelEvaluators::IB5{},
                                    1.75,
                                    IBKernelEvaluators::Weights<double, 5>{ 0.0233383448809079538,
                                                                            0.333140121904304905,
                                                                            0.514244366143986938,
                                                                            0.128737522475479593,
                                                                            0.000539644595320609716 }),
                       sample_error(IBKernelEvaluators::IB6{},
                                    2.25,
                                    IBKernelEvaluators::Weights<double, 6>{ 0.00965617417165844278,
                                                                            0.174648694040214713,
                                                                            0.431221688477088836,
                                                                            0.325168575099164853,
                                                                            0.0591221373512527211,
                                                                            0.000182730860620434541 }),
                       sample_error(IBKernelEvaluators::IB6{},
                                    2.75,
                                    IBKernelEvaluators::Weights<double, 6>{ 0.000182730860620434541,
                                                                            0.0591221373512527211,
                                                                            0.325168575099164853,
                                                                            0.431221688477088836,
                                                                            0.174648694040214713,
                                                                            0.00965617417165844278 }) });
    const double moments = std::max({ moment_error(IBKernelEvaluators::BSpline<2>{}),
                                      moment_error(IBKernelEvaluators::BSpline<3>{}),
                                      moment_error(IBKernelEvaluators::BSpline<4>{}),
                                      moment_error(IBKernelEvaluators::BSpline<5>{}),
                                      moment_error(IBKernelEvaluators::BSpline<6>{}),
                                      moment_error(IBKernelEvaluators::BSpline<7>{}),
                                      moment_error(IBKernelEvaluators::BSpline<8>{}),
                                      moment_error(IBKernelEvaluators::IB3{}),
                                      moment_error(IBKernelEvaluators::IB4{}),
                                      moment_error(IBKernelEvaluators::IB5{}),
                                      moment_error(IBKernelEvaluators::IB6{}) });

    return error > 1.0e-12 || moments > 1.0e-12;
}

template <int Axis>
double
tensor_product_error()
{
    using namespace IBTK;
    const IBKernelEvaluatorTensorProduct product{ IBKernelEvaluators::IB4{}, IBKernelEvaluators::IB3{} };
    std::array<double, NDIM> r;
    r.fill(1.0);
    r[Axis] = 1.5;
    const IBKernelEvaluators::Weights<double, NDIM == 2 ? 12 : 36> weights =
        product.template evaluate<Axis, IBKernelEvaluators::Weights<double, NDIM == 2 ? 12 : 36>>(r);
    constexpr std::array<std::size_t, NDIM> widths = product.template get_stencil_widths<Axis>();
    static_assert(weights.size() == (NDIM == 2 ? 12 : 36), "Natural tensor stencil size");
    const double a = (2.0 - std::sqrt(2.0)) / 8.0, b = (2.0 + std::sqrt(2.0)) / 8.0;
    const IBKernelEvaluators::Weights<double, 4> normal = { a, b, b, a };
    const IBKernelEvaluators::Weights<double, 3> tangent = { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
    double error = 0.0;
    for (std::size_t entry = 0; entry < weights.size(); ++entry)
    {
        std::size_t index = entry;
        double expected = 1.0;
        for (int d = 0; d < NDIM; ++d)
        {
            const int j = index % widths[d];
            index /= widths[d];
            expected *= d == Axis ? normal[j] : tangent[j];
        }
        const double entry_error = std::abs(weights[entry] - expected);
        if (!(entry_error <= 1.0e-12))
        {
            TBOX_ERROR("Tensor weight error = " << entry_error << '\n');
        }
        error = std::max(error, entry_error);
    }
    return error;
}

double
check_tensor_products()
{
    using namespace IBTK;
    TBOX_ASSERT(check_kernels() == 0);
    const IBKernelEvaluatorTensorProduct float_only{ FloatOnlyScalar{} };
    const IBKernelEvaluators::Weights<float, 1> float_product =
        float_only.template evaluate<0, IBKernelEvaluators::Weights<float, 1>>(std::array<double, NDIM>{});
    TBOX_ASSERT(float_product[0] == std::ldexp(1.0f, -NDIM));
    const IBKernelEvaluatorTensorProduct explicit_copy{ ExplicitConstructionScalar{ 0.5 } };
    const IBKernelEvaluatorTensorProduct explicit_moves{ ExplicitConstructionScalar{ 0.25 },
                                                         ExplicitConstructionScalar{ 0.5 } };
    const IBKernelEvaluators::Weights<double, 1> copied =
        explicit_copy.template evaluate<0, IBKernelEvaluators::Weights<double, 1>>(std::array<double, NDIM>{});
    const IBKernelEvaluators::Weights<double, 1> moved =
        explicit_moves.template evaluate<NDIM - 1, IBKernelEvaluators::Weights<double, 1>>(std::array<double, NDIM>{});
    TBOX_ASSERT(copied[0] == std::ldexp(1.0, -NDIM));
    TBOX_ASSERT(moved[0] == std::ldexp(1.0, -NDIM - 1));
    double error = std::max(tensor_product_error<0>(), tensor_product_error<1>());
#if NDIM == 3
    error = std::max(error, tensor_product_error<2>());
#endif
    const IBKernelEvaluatorTensorProduct bspline3{ IBKernelEvaluators::BSpline<3>{} };
    const IBKernelEvaluatorTensorProduct bspline5{ IBKernelEvaluators::BSpline<5>{} };
    std::array<double, NDIM> r;
    r.fill(1.0);
    const IBKernelEvaluators::Weights<double, NDIM == 2 ? 9 : 27> weights3 =
        bspline3.template evaluate<0, IBKernelEvaluators::Weights<double, NDIM == 2 ? 9 : 27>>(r);
    r.fill(1.5);
    const IBKernelEvaluators::Weights<double, NDIM == 2 ? 25 : 125> weights5 =
        bspline5.template evaluate<NDIM - 1, IBKernelEvaluators::Weights<double, NDIM == 2 ? 25 : 125>>(r);
    static_assert(weights3.size() == (NDIM == 2 ? 9 : 27), "Natural three-point tensor stencil size");
    static_assert(weights5.size() == (NDIM == 2 ? 25 : 125), "Natural five-point tensor stencil size");
    const double error3 = std::abs(weights3[weights3.size() / 2] - std::pow(0.75, NDIM));
    const double error5 = std::abs(weights5[0] - std::pow(1.0 / 24.0, NDIM));
    if (!(error3 <= 1.0e-12 && error5 <= 1.0e-12))
    {
        TBOX_ERROR("B-spline tensor weight errors = " << error3 << ", " << error5 << '\n');
    }
    error = std::max({ error, error3, error5 });
    return error;
}
} // namespace

int
main(int argc, char* argv[])
{
    static_assert(std::is_trivially_copyable<IBKernel>::value, "identity copies must be trivial");
    static_assert(std::is_standard_layout<IBKernel>::value, "test inspects the two-block layout");
    static_assert(std::is_convertible<std::string, IBKernelTensorProduct>::value,
                  "strings must convert to tensor products");
    static_assert(std::is_convertible<const char*, IBKernelTensorProduct>::value,
                  "C strings must convert to tensor products");
    static_assert(std::is_convertible<IBKernel, IBKernelTensorProduct>::value,
                  "scalar kernels must convert to tensor products");
    static_assert(!std::is_default_constructible<IBKernel>::value, "scalar kernels require explicit construction");
    static_assert(!std::is_default_constructible<IBKernelTensorProduct>::value,
                  "tensor products require explicit construction");
    static_assert(!std::is_constructible<IBKernelTensorProduct, std::vector<IBKernel>>::value,
                  "dynamic factor containers are not a public construction surface");
    static_assert(!std::is_constructible<IBKernelTensorProduct, std::array<IBKernel, NDIM>>::value,
                  "exact arrays are not a public construction surface");
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const std::string input_file = argc > 1 ? argv[1] : "";
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Logger::Appender> abort_appender = new TestAppender();
    SAMRAI::tbox::Logger::getInstance()->setAbortAppender(abort_appender);
    if (input_file.find("invalid_scalar") != std::string::npos)
    {
        SAMRAI::tbox::PIO::logOnlyNodeZero("output");
        const IBKernel invalid("BAD NAME");
        return 0;
    }

    // Ranks construct shared names in opposite orders with disjoint extras.
    // Compare the actual integers via typed MPI transport, not decoded names
    // or host-order bytes. The serial fixtures check the same value path.
    constexpr const char* COMMON_NAMES[] = {
        "IB_4", "APPLICATION_KERNEL", "ABCDEFGHIJKLMNOPQRSTUVWX", "custom_unknown"
    };
    constexpr std::size_t COMMON_NAME_COUNT = std::size(COMMON_NAMES);
    std::array<std::uint64_t, ENCODED_BLOCK_COUNT * COMMON_NAME_COUNT> common_blocks{};
    for (std::size_t j = 0; j < COMMON_NAME_COUNT; ++j)
    {
        const std::size_t i = rank % 2 ? COMMON_NAME_COUNT - 1 - j : j;
        const IBKernel extra((rank % 2 ? "ODD_EXTRA_" : "EVEN_EXTRA_") + std::to_string(j));
        const IBKernel common(COMMON_NAMES[i]);
        const std::array<std::uint64_t, ENCODED_BLOCK_COUNT> blocks = numeric_blocks(common);
        for (std::size_t block = 0; block < ENCODED_BLOCK_COUNT; ++block)
        {
            common_blocks[ENCODED_BLOCK_COUNT * i + block] = blocks[block];
        }
    }
    auto root_blocks = common_blocks;
    MPI_Bcast(root_blocks.data(), static_cast<int>(root_blocks.size()), MPI_UINT64_T, 0, MPI_COMM_WORLD);
    TBOX_ASSERT(common_blocks == root_blocks);

    // Independently tabulated values for the canonical letter, digit, and
    // underscore encoding with fixed-width block padding.
    struct Encoding
    {
        const char* name;
        std::array<std::uint64_t, ENCODED_BLOCK_COUNT> blocks;
    };
    const Encoding encodings[] = {
        { "A", { { UINT64_C(238572050223552512), 0 } } },
        { "_", { { UINT64_C(8827165858271442944), 0 } } },
        { "ABCDEFGHIJKL", { { UINT64_C(251642104107238734), 0 } } },
        { "ABCDEFGHIJKLM", { { UINT64_C(251642104107238734), UINT64_C(3101436652906182656) } } },
        { "ABCDEFGHIJKLMNOPQRSTUVWX", { { UINT64_C(251642104107238734), UINT64_C(3191881425781291314) } } },
        { "ABCDEFGHIJKLMNOPQRSTUVW_", { { UINT64_C(251642104107238734), UINT64_C(3191881425781291327) } } },
        { "USER_DEFINED", { { UINT64_C(5130207666400688530), 0 } } },
        { "IB_4", { { UINT64_C(2165952653010967808), 0 } } }
    };
    for (const auto& expected : encodings)
    {
        const IBKernel kernel(expected.name), copy(kernel);
        TBOX_ASSERT(numeric_blocks(kernel) == expected.blocks);
        TBOX_ASSERT(copy == kernel && numeric_blocks(copy) == expected.blocks);
        TBOX_ASSERT(copy.getName() == expected.name);
    }
    TBOX_ASSERT(IBKernel("ABCDEFGHIJKL") != IBKernel("ABCDEFGHIJKLM"));
    TBOX_ASSERT(IBKernel("ABCDEFGHIJKLMNOPQRSTUVWX") != IBKernel("ABCDEFGHIJKLMNOPQRSTUVW_"));
    TBOX_ASSERT(IBKernel("A") != IBKernel("A_") && IBKernel("A") != IBKernel("AA"));
    const char* const scalar_c_string = "piecewise_constant";
    TBOX_ASSERT(IBKernel{ scalar_c_string } == IBKernel(IBKernel::BSPLINE_1));
    TBOX_ASSERT(!IBKernel::is_valid_name("ABCDEFGHIJKLMNOPQRSTUVWXY"));
    for (const std::string& invalid :
         { std::string("A B"), std::string("A-B"), std::string("A\0B", 3), std::string("A\x80", 2) })
    {
        TBOX_ASSERT(!IBKernel::is_valid_name(invalid));
    }

    const IBKernel standard_values[] = { IBKernel::BSPLINE_1, IBKernel::BSPLINE_2, IBKernel::BSPLINE_3,
                                         IBKernel::BSPLINE_4, IBKernel::BSPLINE_5, IBKernel::BSPLINE_6,
                                         IBKernel::IB_3,      IBKernel::IB_4,      IBKernel::IB_4_W8,
                                         IBKernel::IB_5,      IBKernel::IB_6,      IBKernel::PIECEWISE_CUBIC };
    const char* scalar_names[] = { "BSPLINE_1", "BSPLINE_2", "BSPLINE_3", "BSPLINE_4", "BSPLINE_5", "BSPLINE_6",
                                   "IB_3",      "IB_4",      "IB_4_W8",   "IB_5",      "IB_6",      "PIECEWISE_CUBIC" };
    const std::vector<IBKernel>& catalog = IBKernel::get_standard_kernels();
    const std::size_t n_standard_values = sizeof(standard_values) / sizeof(standard_values[0]);
    TBOX_ASSERT(catalog.size() == n_standard_values);
    TBOX_ASSERT(std::is_sorted(catalog.begin(), catalog.end()));
    TBOX_ASSERT(std::set<IBKernel>(catalog.begin(), catalog.end()).size() == catalog.size());
    for (std::size_t i = 0; i < n_standard_values; ++i)
    {
        TBOX_ASSERT(standard_values[i] == IBKernel(scalar_names[i]));
        TBOX_ASSERT(catalog[i] == standard_values[i]);
        TBOX_ASSERT(catalog[i].getName() != "USER_DEFINED" && catalog[i].getName() != "PIECEWISE_CONSTANT" &&
                    catalog[i].getName() != "PIECEWISE_LINEAR");
        for (std::size_t j = 0; j < n_standard_values; ++j)
        {
            TBOX_ASSERT((standard_values[i] < standard_values[j]) == (i < j));
        }
    }
    TBOX_ASSERT(ib_kernel_static_initialization_valid());
    for (const char* name : { "INVALID", "UNSUPPORTED", "UNKNOWN", "UNKNOWN_CUSTOM" })
    {
        TBOX_ASSERT(IBKernel(name).getName() == name);
    }
    for (const char* name : scalar_names)
    {
        for (const std::string& spelling : spellings(name))
        {
            const IBKernel kernel(spelling);
            TBOX_ASSERT(kernel.getName() == name);
            TBOX_ASSERT(kernel == IBKernel(name));
            TBOX_ASSERT(IBKernelTensorProduct(spelling) == IBKernelTensorProduct({ kernel }));
            TBOX_ASSERT(IBKernelTensorProduct(spelling).size() == 1);
        }
    }

    const std::pair<const char*, const char*> aliases[] = { { "PIECEWISE_CONSTANT", "BSPLINE_1" },
                                                            { "PIECEWISE_LINEAR", "BSPLINE_2" } };
    for (const auto& alias : aliases)
    {
        for (const std::string& spelling : spellings(alias.first))
        {
            TBOX_ASSERT(IBKernel(spelling).getName() == alias.second);
            TBOX_ASSERT(IBKernelTensorProduct(spelling) == IBKernelTensorProduct(alias.second));
        }
    }

    struct Composite
    {
        const char* name;
        const char* normal;
        const char* tangential;
    };
    const Composite composites[] = {
        { "COMPOSITE_BSPLINE_12", "BSPLINE_1", "BSPLINE_2" }, { "COMPOSITE_BSPLINE_21", "BSPLINE_2", "BSPLINE_1" },
        { "COMPOSITE_BSPLINE_23", "BSPLINE_2", "BSPLINE_3" }, { "COMPOSITE_BSPLINE_32", "BSPLINE_3", "BSPLINE_2" },
        { "COMPOSITE_BSPLINE_34", "BSPLINE_3", "BSPLINE_4" }, { "COMPOSITE_BSPLINE_43", "BSPLINE_4", "BSPLINE_3" },
        { "COMPOSITE_BSPLINE_45", "BSPLINE_4", "BSPLINE_5" }, { "COMPOSITE_BSPLINE_54", "BSPLINE_5", "BSPLINE_4" },
        { "COMPOSITE_BSPLINE_56", "BSPLINE_5", "BSPLINE_6" }, { "COMPOSITE_BSPLINE_65", "BSPLINE_6", "BSPLINE_5" },
        { "DISCONTINUOUS_LINEAR", "BSPLINE_2", "BSPLINE_1" }
    };
    for (const auto& expected : composites)
    {
        for (const std::string& spelling : spellings(expected.name))
        {
            const auto product = IBKernelTensorProduct(spelling);
            TBOX_ASSERT(product.size() == 2);
            TBOX_ASSERT(product[0].getName() == expected.normal);
            TBOX_ASSERT(product[1].getName() == expected.tangential);
            TBOX_ASSERT(product == IBKernelTensorProduct({ IBKernel(expected.normal), IBKernel(expected.tangential) }));
            TBOX_ASSERT(product != IBKernelTensorProduct({ IBKernel(expected.tangential), IBKernel(expected.normal) }));
            TBOX_ASSERT(IBKernel(spelling).getName() == std::string(expected.name));
        }
    }

    const char* higher_order_scalars[] = { "BSPLINE_7", "BSPLINE_8", "BSPLINE_12", "BSPLINE_16" };
    for (const char* scalar_name : higher_order_scalars)
    {
        for (const std::string& spelling : spellings(scalar_name))
        {
            const IBKernel kernel(spelling);
            TBOX_ASSERT(kernel.getName() == scalar_name);
            TBOX_ASSERT(IBKernel(spelling.c_str()) == kernel);
            TBOX_ASSERT(IBKernelTensorProduct(spelling) == IBKernelTensorProduct({ kernel }));
            TBOX_ASSERT(IBKernelTensorProduct(spelling.c_str()) == IBKernelTensorProduct({ kernel }));
        }
    }
    const Composite higher_order_composites[] = { { "COMPOSITE_BSPLINE_78", "BSPLINE_7", "BSPLINE_8" },
                                                  { "COMPOSITE_BSPLINE_87", "BSPLINE_8", "BSPLINE_7" },
                                                  { "COMPOSITE_BSPLINE_29", "BSPLINE_2", "BSPLINE_9" },
                                                  { "COMPOSITE_BSPLINE_88", "BSPLINE_8", "BSPLINE_8" },
                                                  { "COMPOSITE_BSPLINE_7_8", "BSPLINE_7", "BSPLINE_8" },
                                                  { "COMPOSITE_BSPLINE_8_7", "BSPLINE_8", "BSPLINE_7" },
                                                  { "COMPOSITE_BSPLINE_2_9", "BSPLINE_2", "BSPLINE_9" },
                                                  { "COMPOSITE_BSPLINE_8_8", "BSPLINE_8", "BSPLINE_8" },
                                                  { "COMPOSITE_BSPLINE_12_11", "BSPLINE_12", "BSPLINE_11" },
                                                  { "COMPOSITE_BSPLINE_11_12", "BSPLINE_11", "BSPLINE_12" },
                                                  { "COMPOSITE_BSPLINE_10_11", "BSPLINE_10", "BSPLINE_11" } };
    for (const auto& expected : higher_order_composites)
    {
        for (const std::string& spelling : spellings(expected.name))
        {
            const auto product = IBKernelTensorProduct(spelling);
            const IBKernel normal(expected.normal), tangential(expected.tangential);
            TBOX_ASSERT(IBKernelTensorProduct::is_valid_name(spelling));
            TBOX_ASSERT(product == IBKernelTensorProduct(spelling.c_str()));
            TBOX_ASSERT(product == IBKernelTensorProduct({ normal, tangential }));
            TBOX_ASSERT(product[0] == normal);
            TBOX_ASSERT(product.size() == (normal == tangential ? 1 : 2));
            if (!product.isIsotropic())
            {
                TBOX_ASSERT(product[1] == tangential);
            }
        }
    }

    std::string name = "ApplicationKernel";
    const IBKernel custom(name);
    name = "changed";
    TBOX_ASSERT(custom.getName() == "APPLICATIONKERNEL");
    for (const std::string& spelling : spellings("APPLICATIONKERNEL"))
    {
        TBOX_ASSERT(custom == IBKernel(spelling));
    }
    TBOX_ASSERT(custom != IBKernel("USER_DEFINED"));
    TBOX_ASSERT(IBKernel("IB_4_custom").getName() == "IB_4_CUSTOM");
    TBOX_ASSERT(IBKernel("composite_bspline_custom").getName() == "COMPOSITE_BSPLINE_CUSTOM");
    TBOX_ASSERT(IBKernelTensorProduct("composite_bspline_custom") ==
                IBKernelTensorProduct({ IBKernel("COMPOSITE_BSPLINE_CUSTOM") }));
    TBOX_ASSERT(IBKernelTensorProduct(custom.getName()) == IBKernelTensorProduct({ custom }));
    IBKernel normal("ApplicationKernel"), tangential("AnotherKernel");
    const IBKernelTensorProduct product({ normal, tangential });
    normal = IBKernel("changed");
    tangential = IBKernel("changed");
    TBOX_ASSERT(product[0] == custom && product[1] == IBKernel("AnotherKernel"));
    TBOX_ASSERT(IBKernelTensorProduct({ custom }) == IBKernelTensorProduct({ IBKernel("applicationkernel") }));
    const IBKernelTensorProduct c_string_name{ "DISCONTINUOUS_LINEAR" };
    TBOX_ASSERT(c_string_name == IBKernelTensorProduct({ IBKernel::BSPLINE_2, IBKernel::BSPLINE_1 }));
    TBOX_ASSERT(IBKernelTensorProduct({ IBKernel("IB_4"), IBKernel("IB_3") }) !=
                IBKernelTensorProduct({ IBKernel("IB_3"), IBKernel("IB_4") }));
    TBOX_ASSERT(IBKernelTensorProduct({ IBKernel("piecewise_constant"), IBKernel("BSPLINE_1") }) ==
                IBKernelTensorProduct({ IBKernel("BSPLINE_1"), IBKernel("BSPLINE_1") }));

    TBOX_ASSERT(!IBKernel::is_valid_name("") && !IBKernelTensorProduct::is_valid_name(""));

    const IBKernel runtime_normal("IB_4"), runtime_tangential("IB_3");
    const IBKernelTensorProduct scalar(IBKernel::IB_4), single({ runtime_normal });
    const IBKernelTensorProduct repeated = { IBKernel::IB_4, IBKernel::IB_4 };
    const IBKernelTensorProduct pair({ runtime_normal, runtime_tangential });
    const IBKernelTensorProduct ordered({ runtime_normal, runtime_tangential });
    TBOX_ASSERT(scalar.size() == 1 && single.size() == 1 && repeated.size() == 1 && pair.size() == 2);
    TBOX_ASSERT(scalar == single && single == repeated && scalar == repeated);
    TBOX_ASSERT(scalar == IBKernel::IB_4 && IBKernel::IB_4 == scalar);
    TBOX_ASSERT(pair != IBKernel::IB_4 && IBKernel::IB_4 != pair);
    const std::string scalar_name = "IB_4";
    TBOX_ASSERT(scalar == scalar_name && scalar_name == scalar);
    TBOX_ASSERT(scalar == "IB_4" && "IB_4" == scalar);
    TBOX_ASSERT(pair != scalar_name && scalar_name != pair);
    TBOX_ASSERT(pair != "IB_4" && "IB_4" != pair);
    TBOX_ASSERT(single.isIsotropic() && repeated.isIsotropic() && !pair.isIsotropic());
    const IBKernelTensorProduct ordered_products[] = { IBKernelTensorProduct(IBKernel::BSPLINE_1),
                                                       IBKernelTensorProduct(
                                                           { IBKernel::BSPLINE_1, IBKernel::BSPLINE_2 }),
                                                       IBKernelTensorProduct(IBKernel::BSPLINE_2) };
    TBOX_ASSERT(std::is_sorted(std::begin(ordered_products), std::end(ordered_products)));
    TBOX_ASSERT(std::set<IBKernelTensorProduct>(std::begin(ordered_products), std::end(ordered_products)).size() ==
                std::size(ordered_products));
    TBOX_ASSERT(!(single < repeated) && !(repeated < single));
    const IBKernelTensorProduct prefix(IBKernel::IB_4);
    const IBKernelTensorProduct longer({ IBKernel::IB_4, IBKernel::IB_3 });
    TBOX_ASSERT(prefix < longer && longer > prefix);
    TBOX_ASSERT(prefix <= repeated && prefix >= repeated);
    TBOX_ASSERT(copy_product(IBKernel::IB_4) == scalar);
    TBOX_ASSERT(copy_product(std::string("IB_4")) == scalar);
    TBOX_ASSERT(copy_product("IB_4") == scalar);
    std::ostringstream diagnostic;
    diagnostic << single << ' ' << pair;
    TBOX_ASSERT(diagnostic.str() == "(IB_4) (IB_4,IB_3)");
    for (std::size_t slot = 0; slot < ordered.size(); ++slot)
    {
        const IBKernel expected = slot == 0 ? runtime_normal : runtime_tangential;
        TBOX_ASSERT(ordered[slot] == expected);
    }

    const double tensor_error = check_tensor_products();
    if (rank == 0)
    {
        std::ofstream out("output");
        TBOX_ASSERT(static_cast<bool>(out));
        for (const char* scalar_name : higher_order_scalars)
        {
            out << scalar_name << " = " << IBKernelTensorProduct(scalar_name) << '\n';
        }
        for (const auto& composite : higher_order_composites)
        {
            out << composite.name << " = " << IBKernelTensorProduct(composite.name) << '\n';
        }
        out << "tensor_product_max_error = " << tensor_error << '\n';
    }
    return !(tensor_error <= 1.0e-12);
}
