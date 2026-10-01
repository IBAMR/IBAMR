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

#ifndef included_IBTK_CartGridPatchwiseFunction_inl
#define included_IBTK_CartGridPatchwiseFunction_inl

#include <ibtk/config.h>

#include <ibtk/CartGridPatchwiseFunction.h>

#include <utility>

namespace IBTK
{
template <PatchwiseCallback Function>
CartGridPatchwiseFunction<Function>::CartGridPatchwiseFunction(std::string object_name, Function function)
    : CartGridFunction(std::move(object_name)), d_function(std::move(function))
{
}

template <PatchwiseCallback Function>
bool
CartGridPatchwiseFunction<Function>::isTimeDependent() const
{
    return true;
}

template <PatchwiseCallback Function>
void
CartGridPatchwiseFunction<Function>::setDataOnPatch(const int data_idx,
                                                    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                                    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
                                                    const double data_time,
                                                    const bool initial_time,
                                                    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> patch_level)
{
    std::invoke(d_function, data_idx, var, patch, data_time, initial_time, patch_level);
}

template <typename Function>
SAMRAI::tbox::Pointer<CartGridFunction>
make_cart_grid_patchwise_function(std::string object_name, Function&& function)
    requires(PatchwiseCallback<std::decay_t<Function>>&& std::constructible_from<std::decay_t<Function>, Function&&>)
{
    return new CartGridPatchwiseFunction<std::decay_t<Function>>(std::move(object_name),
                                                                 std::forward<Function>(function));
}
} // namespace IBTK

#endif // #ifndef included_IBTK_CartGridPatchwiseFunction_inl
