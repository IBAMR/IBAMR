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

#ifndef included_IBTK_CartGridPointwiseFunction_inl
#define included_IBTK_CartGridPointwiseFunction_inl

#include <ibtk/config.h>

#include <ibtk/CartGridPointwiseFunction.h>

#include <CartesianPatchGeometry.h>
#include <Variable.h>

#include <functional>
#include <utility>

namespace IBTK
{
template <>
inline TensorStorage
string_to_enum<TensorStorage>(const std::string& value)
{
    if (strcasecmp(value.c_str(), "FULL") == 0)
    {
        return TensorStorage::FULL;
    }
    if (strcasecmp(value.c_str(), "SYMMETRIC") == 0)
    {
        return TensorStorage::SYMMETRIC;
    }
    TBOX_ERROR("Unknown TensorStorage: " << value << '\n');
    return TensorStorage::FULL;
}

template <>
inline std::string
enum_to_string<TensorStorage>(const TensorStorage value)
{
    switch (value)
    {
    case TensorStorage::FULL:
        return "FULL";
    case TensorStorage::SYMMETRIC:
        return "SYMMETRIC";
    default:
        TBOX_ERROR("Unknown TensorStorage\n");
    }
    return "";
}

template <PointwiseValue Value, Centering Layout, PointwiseCallback<Value> Function>
CartGridPointwiseFunction<Value, Layout, Function>::CartGridPointwiseFunction(std::string object_name,
                                                                              Function function)
    requires(!std::same_as<Value, MatrixNd>)
    : CartGridFunction(std::move(object_name)), d_function(std::move(function)), d_tensor_storage(TensorStorage::FULL)
{
}

template <PointwiseValue Value, Centering Layout, PointwiseCallback<Value> Function>
CartGridPointwiseFunction<Value, Layout, Function>::CartGridPointwiseFunction(
    std::string object_name,
    Function function,
    const TensorStorage storage) requires std::same_as<Value, MatrixNd>
    : CartGridFunction(std::move(object_name)), d_function(std::move(function)), d_tensor_storage(storage)
{
    switch (storage)
    {
    case TensorStorage::FULL:
    case TensorStorage::SYMMETRIC:
        break;
    default:
        TBOX_ERROR("CartGridPointwiseFunction: invalid tensor storage\n");
    }
}

template <PointwiseValue Value, Centering Layout, PointwiseCallback<Value> Function>
bool
CartGridPointwiseFunction<Value, Layout, Function>::isTimeDependent() const
{
    return true;
}

template <PointwiseValue Value, Centering Layout, PointwiseCallback<Value> Function>
void
CartGridPointwiseFunction<Value, Layout, Function>::setDataOnPatch(
    const int data_idx,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> /*var*/,
    SAMRAI::tbox::Pointer<SAMRAI::hier::Patch<NDIM>> patch,
    const double data_time,
    const bool /*initial_time*/,
    SAMRAI::tbox::Pointer<SAMRAI::hier::PatchLevel<NDIM>> /*patch_level*/)
{
    using namespace SAMRAI;
    if (!patch)
    {
        TBOX_ERROR("CartGridPointwiseFunction: a patch is required\n");
    }
    const tbox::Pointer<hier::PatchData<NDIM>> patch_data = patch->getPatchData(data_idx);
#if !defined(NDEBUG)
    if (!dynamic_cast<Data*>(patch_data.getPointer()))
    {
        TBOX_ERROR("CartGridPointwiseFunction: incompatible or missing patch data\n");
    }
#endif
    Data* const data = static_cast<Data*>(patch_data.getPointer());
    const tbox::Pointer<geom::CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
    if (!geometry)
    {
        TBOX_ERROR("CartGridPointwiseFunction: Cartesian geometry is required\n");
    }
    applyPointwise(*data, patch->getBox(), *geometry, data_time);
}

template <PointwiseValue Value, Centering Layout, PointwiseCallback<Value> Function>
void
CartGridPointwiseFunction<Value, Layout, Function>::applyPointwise(
    Data& data,
    const SAMRAI::hier::Box<NDIM>& box,
    const SAMRAI::geom::CartesianPatchGeometry<NDIM>& geometry,
    const double time)
{
    const int depth = data.getDepth();
    if (depth <= 0)
    {
        TBOX_ERROR("CartGridPointwiseFunction: patch data must have positive depth\n");
    }
    if constexpr (std::same_as<Value, VectorNd>)
    {
        if (depth != NDIM)
        {
            TBOX_ERROR("CartGridPointwiseFunction: VectorNd requires depth NDIM\n");
        }
    }
    else if constexpr (std::same_as<Value, MatrixNd>)
    {
        if (depth != (d_tensor_storage == TensorStorage::FULL ? NDIM * NDIM : NDIM * (NDIM + 1) / 2))
        {
            TBOX_ERROR("CartGridPointwiseFunction: tensor storage does not match patch data depth\n");
        }
    }
    Value q{}, result{};
    if constexpr (std::same_as<Value, VectorXd>)
    {
        q.resize(depth);
        result.resize(depth);
    }
    const auto component = [&](Value& value, const int d) -> double&
    {
        if constexpr (std::same_as<Value, double>)
        {
            return value;
        }
        else if constexpr (std::same_as<Value, MatrixNd>)
        {
            const std::pair<int, int> ij = d_tensor_storage == TensorStorage::FULL ?
                                               std::pair<int, int>{ d / NDIM, d % NDIM } :
                                               voigt_to_tensor_idx(d);
            return value(ij.first, ij.second);
        }
        else
        {
            return value[d];
        }
    };
    const double* const x_lower = geometry.getXLower();
    const double* const dx = geometry.getDx();
    const SAMRAI::hier::Index<NDIM>& index_lower = box.lower();
    const int n_groups = std::same_as<Value, double> ? depth : 1;
    const int n_components = std::same_as<Value, double> ? 1 : depth;
    for (int orientation = 0; orientation < (Layout::is_oriented() ? NDIM : 1); ++orientation)
    {
        if (!Layout::template has_axis<double>(data, orientation))
        {
            continue;
        }
        const int axis = Layout::is_oriented() ? orientation : invalid_index;
        const VectorNd offset = Layout::offset(orientation);
        for (auto it = Layout::begin(box, orientation); it; it++)
        {
            const typename Layout::Index& index = it();
            const SAMRAI::hier::Index<NDIM> cartesian_index = Layout::cartesian_index(index);
            VectorNd x;
            for (int d = 0; d < NDIM; ++d)
            {
                x[d] = x_lower[d] + dx[d] * (cartesian_index(d) - index_lower(d) + offset[d]);
            }
            for (int group = 0; group < n_groups; ++group)
            {
                auto evaluate = [&]() -> decltype(auto)
                {
                    if constexpr (
                        std::is_invocable_r_v<Value, Function&, const Value&, const VectorNd&, double, int, int>)
                    {
                        for (int d = 0; d < n_components; ++d)
                        {
                            component(q, d) = data(index, group + d);
                        }
                        if constexpr (std::same_as<Value, MatrixNd>)
                        {
                            if (d_tensor_storage == TensorStorage::SYMMETRIC)
                            {
                                for (int i = 1; i < NDIM; ++i)
                                {
                                    for (int j = 0; j < i; ++j)
                                    {
                                        q(i, j) = q(j, i);
                                    }
                                }
                            }
                        }
                        return std::invoke(d_function, std::as_const(q), std::as_const(x), time, group, axis);
                    }
                    else
                    {
                        return std::invoke(d_function, std::as_const(x), time, group, axis);
                    }
                };
                decltype(auto) value = evaluate();
                if constexpr (!std::same_as<Value, double>)
                {
                    if (value.rows() != q.rows() || value.cols() != q.cols())
                    {
                        TBOX_ERROR("CartGridPointwiseFunction: callback result has an incompatible shape\n");
                    }
                }
                // Materialize lazy expressions before scattering or overwriting q.
                result = std::forward<decltype(value)>(value);
                if constexpr (std::same_as<Value, MatrixNd>)
                {
                    if (d_tensor_storage == TensorStorage::SYMMETRIC && !result.isApprox(result.transpose()))
                    {
                        TBOX_ERROR(
                            "CartGridPointwiseFunction: symmetric storage requires a symmetric callback result\n");
                    }
                }
                for (int d = 0; d < n_components; ++d)
                {
                    data(index, group + d) = component(result, d);
                }
            }
        }
    }
}

namespace detail
{
template <typename Value, typename Function, typename... Args>
SAMRAI::tbox::Pointer<CartGridFunction>
allocate_cart_grid_pointwise_function(std::string object_name,
                                      SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                      Function&& function,
                                      Args&&... args)
{
    if (!var)
    {
        TBOX_ERROR("CartGridPointwiseFunction: a variable is required\n");
    }
    const DataCentering centering = get_data_centering<double>(*var->getPatchDataFactory());
    if (centering == DataCentering::UNKNOWN)
    {
        TBOX_ERROR("CartGridPointwiseFunction: unsupported variable patch data factory\n");
    }
    return dispatch_data_centering(
        centering,
        [&]<Centering Layout>() -> SAMRAI::tbox::Pointer<CartGridFunction>
        {
            return new CartGridPointwiseFunction<Value, Layout, std::decay_t<Function>>(
                std::move(object_name), std::forward<Function>(function), std::forward<Args>(args)...);
        });
}
} // namespace detail

template <PointwiseValue Value, typename Function>
SAMRAI::tbox::Pointer<CartGridFunction>
make_cart_grid_pointwise_function(std::string object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                  Function&& function)
    requires(!std::same_as<Value, MatrixNd> && PointwiseCallback<std::decay_t<Function>, Value>)
{
    return detail::allocate_cart_grid_pointwise_function<Value>(
        std::move(object_name), var, std::forward<Function>(function));
}

template <PointwiseValue Value, typename Function>
SAMRAI::tbox::Pointer<CartGridFunction>
make_cart_grid_pointwise_function(std::string object_name,
                                  SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM>> var,
                                  Function&& function,
                                  const TensorStorage storage)
    requires(std::same_as<Value, MatrixNd>&& PointwiseCallback<std::decay_t<Function>, Value>)
{
    return detail::allocate_cart_grid_pointwise_function<Value>(
        std::move(object_name), var, std::forward<Function>(function), storage);
}
} // namespace IBTK

#endif // #ifndef included_IBTK_CartGridPointwiseFunction_inl
