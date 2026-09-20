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

#ifndef included_IBTK_PETScMatUtilities_inl
#define included_IBTK_PETScMatUtilities_inl

#include <ibtk/config.h>

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/Utilities.h>

#include <Box.h>
#include <Index.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideIndex.h>

#include <algorithm>
#include <array>
#include <limits>
#include <utility>

namespace IBTK
{
template <IBKernelEvaluatorCartesian<double, PetscScalar> Evaluator>
inline void
PETScMatUtilities::constructPatchLevelSCInterpOp(Mat& mat,
                                                 const Evaluator& evaluator,
                                                 Vec X_vec,
                                                 const std::vector<int>& num_dofs_per_proc,
                                                 int dof_index_idx,
                                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level)
{
    constexpr std::array<std::array<std::size_t, NDIM>, NDIM> widths = {
        Evaluator::template get_stencil_widths<0>(),
        Evaluator::template get_stencil_widths<1>()
#if (NDIM == 3)
            ,
        Evaluator::template get_stencil_widths<2>()
#endif
    };
    constexpr std::size_t max_width = [&widths]()
    {
        std::size_t result = 0;
        for (const auto& component_widths : widths)
        {
            for (std::size_t width : component_widths)
            {
                result = std::max(result, width);
            }
        }
        return result;
    }();
    static_assert(max_width <= static_cast<std::size_t>(std::numeric_limits<int>::max()),
                  "Interpolation stencil width exceeds the SAMRAI index range");
    std::array<std::array<int, NDIM>, NDIM> index_widths;
    for (int axis = 0; axis < NDIM; ++axis)
    {
        for (int d = 0; d < NDIM; ++d)
        {
            index_widths[axis][d] = static_cast<int>(widths[axis][d]);
        }
    }
    SCInterpOpData data(mat, X_vec, index_widths, num_dofs_per_proc, dof_index_idx, patch_level);
    constructSCInterpOpAxis<0>(data, evaluator);
    constructSCInterpOpAxis<1>(data, evaluator);
#if (NDIM == 3)
    constructSCInterpOpAxis<2>(data, evaluator);
#endif
    data.assemble();
}

template <int Axis, IBKernelEvaluatorCartesian<double, PetscScalar> Evaluator>
inline void
PETScMatUtilities::constructSCInterpOpAxis(SCInterpOpData& data, const Evaluator& evaluator)
{
    using namespace SAMRAI;
    constexpr std::size_t nvalues = detail::ib_kernel_stencil_size<Evaluator, Axis>();
    using Weights = IBKernelEvaluators::Weights<PetscScalar, nvalues>;
    static_assert(nvalues <= static_cast<std::size_t>(std::numeric_limits<PetscInt>::max()));
    for (int point = 0; point < data.d_n_local_points; ++point)
    {
        const double* const X = &data.d_positions[NDIM * point];
        const hier::Box<NDIM>& box = data.d_stencil_boxes[point][Axis];
        const hier::Index<NDIM>& lower = box.lower();
        std::array<double, NDIM> r;
        for (int d = 0; d < NDIM; ++d)
        {
            const double x_lower =
                (static_cast<double>(lower(d) - data.d_domain_lower(d)) + (d == Axis ? 0.0 : 0.5)) * data.d_dx[d] +
                data.d_x_lower[d];
            r[d] = (X[d] - x_lower) / data.d_dx[d];
        }
        const Weights values = evaluator.template evaluate<Axis, Weights>(std::as_const(r));

        std::array<PetscInt, nvalues> columns;

        const tbox::Pointer<pdat::SideData<NDIM, int>>& indices = data.d_dof_index_data[point];
        int entry = 0;
        bool stencil_outside_domain = false;
        for (typename hier::Box<NDIM>::Iterator b(box); b; b++, ++entry)
        {
            columns[entry] = (*indices)(pdat::SideIndex<NDIM>(b(), Axis, pdat::SideIndex<NDIM>::Lower));
            stencil_outside_domain = stencil_outside_domain || columns[entry] < 0;
        }
        if (stencil_outside_domain)
        {
            // A negative column is outside the domain (e.g. at a non-periodic physical boundary) and has no
            // DOF; MatSetValues() below silently drops that entry instead of adding it, so the row's weights
            // sum to less than one. No folding or renormalization is attempted.
            IBTK_DO_ONCE(TBOX_WARNING("PETScMatUtilities::constructPatchLevelSCInterpOp():\n"
                                      << "  an interpolation stencil extends past a non-periodic physical "
                                         "boundary; the weights of the stencil points outside the domain are "
                                         "dropped, not folded or renormalized.\n"););
        }
        const PetscInt row = data.d_row_lower + NDIM * point + Axis;
        // Periodic stencil points can share a column; sum their contributions.
        const int ierr = MatSetValues(
            data.d_mat, 1, &row, static_cast<PetscInt>(nvalues), columns.data(), values.data(), ADD_VALUES);
        IBTK_CHKERRQ(ierr);
    }
}

} // namespace IBTK

#endif
