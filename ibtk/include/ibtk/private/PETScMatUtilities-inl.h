// ---------------------------------------------------------------------
//
// Copyright (c) 2011 - 2026 by the IBAMR developers
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

#include <Box.h>
#include <Index.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideIndex.h>

#include <array>
#include <tuple>
#include <type_traits>
#include <utility>

namespace IBTK
{
template <class Evaluator>
inline void
PETScMatUtilities::constructPatchLevelSCInterpOp(Mat& mat,
                                                 const Evaluator& evaluator,
                                                 Vec& X_vec,
                                                 const std::vector<int>& num_dofs_per_proc,
                                                 int dof_index_idx,
                                                 SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level)
{
    constexpr std::array<std::array<int, NDIM>, NDIM> widths = {
        Evaluator::template get_stencil_widths<0>(),
        Evaluator::template get_stencil_widths<1>()
#if (NDIM == 3)
            ,
        Evaluator::template get_stencil_widths<2>()
#endif
    };
    static_assert(
        [widths]
        {
            for (const std::array<int, NDIM>& axis_widths : widths)
                for (int width : axis_widths)
                    if (width <= 0) return false;
            return true;
        }(),
        "Evaluator stencil widths must be positive");
    SCInterpOpData data(mat, X_vec, widths, num_dofs_per_proc, dof_index_idx, patch_level);
    construct_sc_interp_op_axis<0>(data, evaluator);
    construct_sc_interp_op_axis<1>(data, evaluator);
#if (NDIM == 3)
    construct_sc_interp_op_axis<2>(data, evaluator);
#endif
    data.assemble();
}

template <int Axis, class Evaluator>
inline void
PETScMatUtilities::construct_sc_interp_op_axis(SCInterpOpData& data, const Evaluator& evaluator)
{
    using namespace SAMRAI;
    constexpr int nvalues = []
    {
        int size = 1;
        for (int width : Evaluator::template get_stencil_widths<Axis>()) size *= width;
        return size;
    }();
    using Values = decltype(std::declval<const Evaluator&>().template evaluate<Axis>(
        std::declval<const std::array<double, NDIM>&>()));
    static_assert(std::is_same_v<Values, std::array<double, nvalues>>,
                  "A const evaluator must return std::array<double, N> with N equal to the stencil size");
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
        const std::array<double, nvalues> values = evaluator.template evaluate<Axis>(r);

        std::array<int, nvalues> columns;

        tbox::Pointer<hier::Patch<NDIM>> patch = data.d_level->getPatch(data.d_patch_numbers[point]);
        tbox::Pointer<pdat::SideData<NDIM, int>> indices = patch->getPatchData(data.d_dof_index_idx);
#if !defined(NDEBUG)
        TBOX_ASSERT(indices->getDepth() == 1);
#endif
        int entry = 0;
        for (typename hier::Box<NDIM>::Iterator b(box); b; b++, ++entry)
            columns[entry] = (*indices)(pdat::SideIndex<NDIM>(b(), Axis, pdat::SideIndex<NDIM>::Lower));
        const int row = data.d_row_lower + NDIM * point + Axis;
        const int ierr = MatSetValues(data.d_mat, 1, &row, nvalues, columns.data(), values.data(), INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
}

} // namespace IBTK

#endif
