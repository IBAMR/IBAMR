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

#ifndef included_IBTK_IBOperatorBuilder_inl
#define included_IBTK_IBOperatorBuilder_inl

#include <ibtk/config.h>

#include <ibtk/IBOperatorBuilder.h>
#include <ibtk/PETScMatUtilities.h>

#include <algorithm>
#include <array>
#include <concepts>
#include <memory>
#include <utility>

namespace IBTK
{
template <class Evaluator>
class IBOperatorBuilder::Model final : public IBOperatorBuilder::Concept
{
public:
    template <class Argument>
    requires std::constructible_from<Evaluator, Argument&&> explicit Model(Argument&& evaluator)
        : d_evaluator(std::forward<Argument>(evaluator))
    {
    }

    int getMinimumGhostWidth() const override
    {
        std::size_t width = 0;
        for (const std::array<std::size_t, NDIM>& widths : {
                 Evaluator::template get_stencil_widths<0>(), Evaluator::template get_stencil_widths<1>()
#if (NDIM == 3)
                                                                  ,
                     Evaluator::template get_stencil_widths<2>()
#endif
             })
        {
            width = std::max(width, *std::max_element(widths.begin(), widths.end()));
        }
        return static_cast<int>((width + 1) / 2 + 1);
    }

    void
    constructInterpolationMatrixSide(Mat& mat,
                                     Vec X_vec,
                                     const std::vector<int>& num_dofs_per_proc,
                                     int dof_index_idx,
                                     SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level) const override
    {
        PETScMatUtilities::constructPatchLevelSCInterpOp(
            mat, d_evaluator, X_vec, num_dofs_per_proc, dof_index_idx, patch_level);
    }

private:
    const Evaluator d_evaluator;
};

template <class Evaluator>
requires(IBKernelEvaluatorCartesian<std::remove_cvref_t<Evaluator>, double, PetscScalar>&&
             std::constructible_from<std::remove_cvref_t<Evaluator>, Evaluator&&>)
    IBOperatorBuilder::IBOperatorBuilder(Evaluator&& evaluator)
    : d_operations(std::make_shared<Model<std::remove_cvref_t<Evaluator>>>(std::forward<Evaluator>(evaluator)))
{
}

inline int
IBOperatorBuilder::getMinimumGhostWidth() const
{
    return d_operations->getMinimumGhostWidth();
}

inline void
IBOperatorBuilder::constructInterpolationMatrixSide(
    Mat& mat,
    Vec X_vec,
    const std::vector<int>& num_dofs_per_proc,
    int dof_index_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level) const
{
    d_operations->constructInterpolationMatrixSide(mat, X_vec, num_dofs_per_proc, dof_index_idx, patch_level);
}
} // namespace IBTK
#endif
