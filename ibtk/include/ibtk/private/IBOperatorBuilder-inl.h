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

#include <memory>
#include <utility>

namespace IBTK
{
template <class Evaluator>
class IBOperatorBuilder::Model final : public IBOperatorBuilder::Concept
{
public:
    explicit Model(Evaluator evaluator) : d_evaluator(std::move(evaluator))
    {
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

template <IBKernelEvaluatorCartesian<double, PetscScalar> Evaluator>
IBOperatorBuilder::IBOperatorBuilder(Evaluator evaluator)
    : d_operations(std::make_shared<Model<Evaluator>>(std::move(evaluator)))
{
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
