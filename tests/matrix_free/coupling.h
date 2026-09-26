// Copyright (c) 2026 by the IBAMR developers
// This file is part of IBAMR and is distributed under the 3-clause BSD license.
#ifndef included_matrix_free_coupling
#define included_matrix_free_coupling

#include <ibtk/config.h>

#include <ibtk/IBKernelEvaluatorTensorProduct.h>
#include <ibtk/ib_kernel_evaluators.h>

#include <SideCoupling.h>

#include <string>

namespace MatrixFreeTest
{
/*! \brief Visit BS2-6 and both adjacent-order composite orientations for orders 1-6. */
template <class Visitor>
void for_each_bspline(const Visitor& visit);

/*! \brief Evaluate an application-defined four-point cosine kernel. */
class CosineKernel
{
public:
    /*! \brief Return the actual support width. */
    static constexpr std::size_t get_stencil_width();

    /*! \brief Evaluate directly in the caller's coefficient precision. */
    template <class Output, class Input>
    Output evaluate(Input r) const;
};

/*! \brief Custom Cartesian evaluator without the optional factorized interface. */
class CartesianCosineKernel
{
public:
    template <int Axis>
    static constexpr std::array<std::size_t, NDIM> get_stencil_widths();

    template <int Axis, class Output, class Input>
    Output evaluate(const std::array<Input, NDIM>& r) const;
};

/*! \brief Apply all vector components to interleaved marker data. */
template <bool Spread,
          class Coefficient = double,
          IBTK::Experimental::TensorProductMode Mode = IBTK::Experimental::TensorProductMode::CONTRACTED,
          class Evaluator>
void couple(const Evaluator& evaluator,
            const SAMRAI::hier::Patch<NDIM>& patch,
            SAMRAI::pdat::SideData<NDIM, double>& field,
            std::span<const double> positions,
            std::span<const int> indices,
            std::span<const double> shifts,
            std::conditional_t<Spread, const double*, double*> values);
} // namespace MatrixFreeTest

#include "coupling-inl.h"
#endif
