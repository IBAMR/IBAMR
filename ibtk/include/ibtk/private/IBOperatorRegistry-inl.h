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

#ifndef included_IBTK_IBOperatorRegistry_inl
#define included_IBTK_IBOperatorRegistry_inl

#include <ibtk/config.h>

#include <ibtk/IBOperatorRegistry.h>
#include <ibtk/PETScMatUtilities.h>

#include <memory>
#include <utility>

namespace IBTK
{
template <class Evaluator>
inline void
IBOperatorRegistry::register_interpolation_matrix_sc(const IBKernelTensorProduct& kernel, Evaluator evaluator)
{
    for (std::size_t d = 0; d < kernel.size(); ++d)
    {
        if (kernel[d] == IBKernel::UNKNOWN)
        {
            TBOX_ERROR("IBOperatorRegistry::register_interpolation_matrix_sc(): unspecified kernel " << kernel << '\n');
        }
    }
    std::map<IBKernelTensorProduct, Builder>& builders = get_builders();
    if (is_supplied_kernel(kernel) || builders.find(kernel) != builders.end())
    {
        TBOX_ERROR("IBOperatorRegistry::register_interpolation_matrix_sc(): kernel " << kernel
                                                                                     << " is already registered\n");
    }
    builders.emplace(kernel, make_builder(std::move(evaluator)));
}

template <class Evaluator>
inline IBOperatorRegistry::Builder
IBOperatorRegistry::make_builder(Evaluator evaluator)
{
    // Compile the complete matrix operation with this evaluator type. Runtime
    // selection calls it once per matrix, not once per stencil coefficient.
    const auto owned_evaluator = std::make_shared<const Evaluator>(std::move(evaluator));
    return [owned_evaluator](Mat& mat,
                             Vec X,
                             const std::vector<int>& num_dofs,
                             int dof_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> level)
    { PETScMatUtilities::constructPatchLevelSCInterpOp(mat, *owned_evaluator, X, num_dofs, dof_idx, level); };
}

} // namespace IBTK
#endif
