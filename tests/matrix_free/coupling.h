// Copyright (c) 2026 by the IBAMR developers
// This file is part of IBAMR and is distributed under the 3-clause BSD license.
#ifndef included_matrix_free_coupling
#define included_matrix_free_coupling

#include <ibtk/config.h>

#include <SideCoupling.h>

namespace MatrixFreeTest
{
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

/*! \brief Apply all vector components to interleaved marker data. */
template <bool Spread, class Coefficient = double, class Evaluator>
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
