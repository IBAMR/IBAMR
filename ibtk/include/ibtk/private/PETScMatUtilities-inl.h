// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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

#include <ibtk/IBKernelEvaluators.h>
#include <ibtk/IBKernelTensorProductEvaluator.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/PETScMatUtilities.h>

#include <Box.h>
#include <Index.h>
#include <Patch.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideIndex.h>

#include <array>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

namespace IBTK
{
// Geometry and borrowed position access for typed interpolation builders.
struct PETScMatUtilities::SCInterpOpData
{
    /*! \brief Allocate the matrix and determine stencil boxes and local patches. */
    SCInterpOpData(Mat& mat,
                   Vec X,
                   int normal_width,
                   int tangential_width,
                   const std::vector<int>& num_dofs_per_proc,
                   int dof_index_idx,
                   SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level);
    /*! \brief Restore the borrowed position array. */
    ~SCInterpOpData();
    /*! \brief Disallow copying borrowed array access. */
    SCInterpOpData(const SCInterpOpData&) = delete;
    /*! \brief Disallow assigning borrowed array access. */
    SCInterpOpData& operator=(const SCInterpOpData&) = delete;
    /*! \brief Finish matrix assembly. */
    void assemble();

    //! Caller-owned matrix handle.
    Mat& d_mat;
    //! Borrowed position vector and its array access.
    Vec d_X;
    double* d_positions = nullptr;
    //! Grid spacings and physical domain origin.
    std::array<double, NDIM> d_dx, d_x_lower;
    //! Lower index of the physical domain.
    SAMRAI::hier::Index<NDIM> d_domain_lower;
    //! Number of local IB points and first local matrix row.
    int d_n_local_points = 0, d_row_lower = 0;
    //! Local patches and component stencil boxes for each IB point.
    std::vector<int> d_patch_numbers;
    std::vector<std::vector<SAMRAI::hier::Box<NDIM>>> d_stencil_boxes;
    //! Borrowed hierarchy data used to read global column indices.
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> d_level;
    int d_dof_index_idx;
};

template <class Evaluator>
inline void
PETScMatUtilities::register_sc_interp_kernel(const IBKernelTensorProduct& kernel, Evaluator evaluator)
{
    auto& builders = get_sc_interp_op_builders();
    if (builders.find(kernel) != builders.end())
    {
        TBOX_ERROR("PETScMatUtilities::register_sc_interp_kernel(): kernel " << kernel << " is already registered\n");
    }
    builders.emplace(kernel, make_sc_interp_op_builder(std::move(evaluator)));
}

template <class Evaluator>
inline PETScMatUtilities::SCInterpOpBuilder
PETScMatUtilities::make_sc_interp_op_builder(Evaluator evaluator)
{
    const auto owned_evaluator = std::make_shared<const Evaluator>(std::move(evaluator));
    return [owned_evaluator](Mat& mat,
                             Vec& X,
                             const std::vector<int>& num_dofs,
                             int dof_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level)
    { construct_sc_interp_op(mat, *owned_evaluator, X, num_dofs, dof_idx, level); };
}

template <class Evaluator>
inline void
PETScMatUtilities::construct_sc_interp_op(Mat& mat,
                                          const Evaluator& evaluator,
                                          Vec& X_vec,
                                          const std::vector<int>& num_dofs_per_proc,
                                          int dof_index_idx,
                                          SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level)
{
    constexpr auto widths = Evaluator::template get_stencil_widths<0, NDIM>();
    SCInterpOpData data(mat, X_vec, widths[0], widths[1], num_dofs_per_proc, dof_index_idx, patch_level);
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
    for (int point = 0; point < data.d_n_local_points; ++point)
    {
        const double* const X = &data.d_positions[NDIM * point];
        const hier::Box<NDIM>& box = data.d_stencil_boxes[point][Axis];
        const auto& lower = box.lower();
        std::array<double, NDIM> r;
        for (int d = 0; d < NDIM; ++d)
        {
            const double x_lower =
                (static_cast<double>(lower(d) - data.d_domain_lower(d)) + (d == Axis ? 0.0 : 0.5)) * data.d_dx[d] +
                data.d_x_lower[d];
            r[d] = (X[d] - x_lower) / data.d_dx[d];
        }
        const auto values = evaluator.template evaluate<Axis>(r);
        constexpr int nvalues = std::tuple_size<decltype(values)>::value;
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
