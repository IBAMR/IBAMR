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

#ifndef included_IBTK_SCInterpOpRegistry_inl
#define included_IBTK_SCInterpOpRegistry_inl

#include <ibtk/config.h>

#include <ibtk/PETScMatUtilities.h>
#include <ibtk/SCInterpOpRegistry.h>

#include <memory>
#include <utility>

namespace IBTK
{
template <class Evaluator>
inline void
SCInterpOpRegistry::register_kernel(const IBKernelTensorProduct& kernel, Evaluator evaluator)
{
    auto& builders = get_builders();
    if (builders.find(kernel) != builders.end())
    {
        TBOX_ERROR("SCInterpOpRegistry::register_kernel(): kernel " << kernel << " is already registered\n");
    }
    builders.emplace(kernel, make_builder(std::move(evaluator)));
}

template <class Evaluator>
inline SCInterpOpRegistry::Builder
SCInterpOpRegistry::make_builder(Evaluator evaluator)
{
    const auto owned_evaluator = std::make_shared<const Evaluator>(std::move(evaluator));
    return [owned_evaluator](Mat& mat,
                             Vec& X,
                             const std::vector<int>& num_dofs,
                             int dof_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level)
    { PETScMatUtilities::constructPatchLevelSCInterpOp(mat, *owned_evaluator, X, num_dofs, dof_idx, level); };
}

} // namespace IBTK
#endif
