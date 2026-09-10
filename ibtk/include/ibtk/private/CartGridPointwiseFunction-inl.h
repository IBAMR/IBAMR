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
#include <CellData.h>
#include <EdgeData.h>
#include <FaceData.h>
#include <NodeData.h>
#include <SideData.h>

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

template <typename Value, typename Function>
CartGridPointwiseFunction<Value, Function>::CartGridPointwiseFunction(std::string object_name, Function function)
    : CartGridFunction(std::move(object_name)), d_function(std::move(function)), d_tensor_storage(TensorStorage::FULL)
{
    static_assert(!std::is_same_v<Value, MatrixNd>, "MatrixNd requires explicit TensorStorage.");
}

template <typename Value, typename Function>
CartGridPointwiseFunction<Value, Function>::CartGridPointwiseFunction(std::string object_name,
                                                                      Function function,
                                                                      const TensorStorage storage)
    : CartGridFunction(std::move(object_name)), d_function(std::move(function)), d_tensor_storage(storage)
{
    static_assert(std::is_same_v<Value, MatrixNd>, "TensorStorage is only applicable to MatrixNd.");
    Values::tensor_depth(storage);
}

template <typename Value, typename Function>
bool
CartGridPointwiseFunction<Value, Function>::isTimeDependent() const
{
    return true;
}

template <typename Value, typename Function>
void
CartGridPointwiseFunction<Value, Function>::setDataOnPatch(
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
    const tbox::Pointer<hier::PatchData<NDIM>> data = patch->getPatchData(data_idx);
    const tbox::Pointer<geom::CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
    if (!data || !geometry)
    {
        TBOX_ERROR("CartGridPointwiseFunction: allocated patch data and Cartesian geometry are required\n");
    }
    // Keep runtime type discovery outside the templated evaluation loop.
    if (tbox::Pointer<pdat::CellData<NDIM, double>> typed_data = data)
    {
        applyPointwise(*typed_data, patch->getBox(), *geometry, data_time);
    }
    else if (tbox::Pointer<pdat::NodeData<NDIM, double>> typed_data = data)
    {
        applyPointwise(*typed_data, patch->getBox(), *geometry, data_time);
    }
    else if (tbox::Pointer<pdat::SideData<NDIM, double>> typed_data = data)
    {
        applyPointwise(*typed_data, patch->getBox(), *geometry, data_time);
    }
    else if (tbox::Pointer<pdat::FaceData<NDIM, double>> typed_data = data)
    {
        applyPointwise(*typed_data, patch->getBox(), *geometry, data_time);
    }
    else if (tbox::Pointer<pdat::EdgeData<NDIM, double>> typed_data = data)
    {
        applyPointwise(*typed_data, patch->getBox(), *geometry, data_time);
    }
    else
    {
        TBOX_ERROR("CartGridPointwiseFunction: unsupported patch data type\n");
    }
}

template <typename Value, typename Function>
template <typename Data>
constexpr bool
CartGridPointwiseFunction<Value, Function>::Centering<Data>::is_staggered()
{
    return !std::is_same_v<Data, SAMRAI::pdat::CellData<NDIM, double>> &&
           !std::is_same_v<Data, SAMRAI::pdat::NodeData<NDIM, double>>;
}

template <typename Value, typename Function>
template <typename Data>
bool
CartGridPointwiseFunction<Value, Function>::Centering<Data>::has_axis(const Data& data, const int axis)
{
    if constexpr (std::is_same_v<Data, SAMRAI::pdat::SideData<NDIM, double>>)
    {
        return data.getDirectionVector()(axis) != 0;
    }
    else
    {
        return true;
    }
}

template <typename Value, typename Function>
template <typename Data>
typename Data::Iterator
CartGridPointwiseFunction<Value, Function>::Centering<Data>::begin(const SAMRAI::hier::Box<NDIM>& box, const int axis)
{
    if constexpr (is_staggered())
    {
        return typename Data::Iterator(box, axis);
    }
    else
    {
        return typename Data::Iterator(box);
    }
}

template <typename Value, typename Function>
template <typename Data>
VectorNd
CartGridPointwiseFunction<Value, Function>::Centering<Data>::offset(const int axis)
{
    VectorNd result;
    if constexpr (std::is_same_v<Data, SAMRAI::pdat::NodeData<NDIM, double>>)
    {
        result.setZero();
    }
    else if constexpr (std::is_same_v<Data, SAMRAI::pdat::EdgeData<NDIM, double>>)
    {
        result.setZero();
        result[axis] = 0.5;
    }
    else
    {
        result.setConstant(0.5);
        if constexpr (is_staggered())
        {
            result[axis] = 0.0;
        }
    }
    return result;
}

template <typename Value, typename Function>
template <typename Data>
template <typename Index>
SAMRAI::hier::Index<NDIM>
CartGridPointwiseFunction<Value, Function>::Centering<Data>::cartesian_index(const Index& index)
{
    if constexpr (std::is_same_v<Data, SAMRAI::pdat::FaceData<NDIM, double>>)
    {
        return index.toCell(1);
    }
    else
    {
        return index;
    }
}

template <typename Value, typename Function>
int
CartGridPointwiseFunction<Value, Function>::Values::tensor_depth(const TensorStorage storage)
{
    switch (storage)
    {
    case TensorStorage::FULL:
        return NDIM * NDIM;
    case TensorStorage::SYMMETRIC:
        return NDIM * (NDIM + 1) / 2;
    default:
        TBOX_ERROR("CartGridPointwiseFunction: invalid tensor storage\n");
    }
    return 0;
}

template <typename Value, typename Function>
void
CartGridPointwiseFunction<Value, Function>::Values::validate_depth(const int depth, const TensorStorage storage)
{
    if (depth <= 0)
    {
        TBOX_ERROR("CartGridPointwiseFunction: patch data must have positive depth\n");
    }
    if constexpr (std::is_same_v<Value, VectorNd>)
    {
        if (depth != NDIM)
        {
            TBOX_ERROR("CartGridPointwiseFunction: VectorNd requires depth NDIM\n");
        }
    }
    else if constexpr (std::is_same_v<Value, MatrixNd>)
    {
        if (depth != (storage == TensorStorage::FULL ? NDIM * NDIM : NDIM * (NDIM + 1) / 2))
        {
            TBOX_ERROR("CartGridPointwiseFunction: tensor storage does not match patch data depth\n");
        }
    }
}

template <typename Value, typename Function>
Value
CartGridPointwiseFunction<Value, Function>::Values::make_value(const int depth)
{
    if constexpr (std::is_same_v<Value, VectorXd>)
    {
        return VectorXd(depth);
    }
    else
    {
        return Value{};
    }
}

template <typename Value, typename Function>
std::pair<int, int>
CartGridPointwiseFunction<Value, Function>::Values::tensor_index(const int component, const TensorStorage storage)
{
    if (storage == TensorStorage::FULL)
    {
        return { component / NDIM, component % NDIM };
    }
    return voigt_to_tensor_idx(component);
}

template <typename Value, typename Function>
template <typename Result>
void
CartGridPointwiseFunction<Value, Function>::Values::assign_result(Value& value, Result&& result, const int depth)
{
    if constexpr (!std::is_same_v<Value, double>)
    {
        const int rows = std::is_same_v<Value, VectorXd> ? depth : NDIM;
        const int cols = std::is_same_v<Value, MatrixNd> ? NDIM : 1;
        if (result.rows() != rows || result.cols() != cols)
        {
            TBOX_ERROR("CartGridPointwiseFunction: callback result has an incompatible shape\n");
        }
    }
    value = std::forward<Result>(result);
}

template <typename Value, typename Function>
template <typename Data, typename Index>
void
CartGridPointwiseFunction<Value, Function>::Values::load(Value& value,
                                                         const Data& data,
                                                         const Index& index,
                                                         const int depth,
                                                         const TensorStorage storage)
{
    if constexpr (std::is_same_v<Value, double>)
    {
        value = data(index, depth);
    }
    else
    {
        for (int d = 0; d < data.getDepth(); ++d)
        {
            if constexpr (std::is_same_v<Value, MatrixNd>)
            {
                const std::pair<int, int> ij = tensor_index(d, storage);
                value(ij.first, ij.second) = data(index, d);
                if (storage == TensorStorage::SYMMETRIC)
                {
                    value(ij.second, ij.first) = value(ij.first, ij.second);
                }
            }
            else
            {
                value[d] = data(index, d);
            }
        }
    }
}

template <typename Value, typename Function>
template <typename Data, typename Index>
void
CartGridPointwiseFunction<Value, Function>::Values::store(const Value& value,
                                                          Data& data,
                                                          const Index& index,
                                                          const int depth,
                                                          const TensorStorage storage)
{
    if constexpr (std::is_same_v<Value, double>)
    {
        data(index, depth) = value;
    }
    else
    {
        if constexpr (std::is_same_v<Value, MatrixNd>)
        {
            if (storage == TensorStorage::SYMMETRIC && !value.isApprox(value.transpose()))
            {
                TBOX_ERROR("CartGridPointwiseFunction: symmetric storage requires a symmetric callback result\n");
            }
        }
        for (int d = 0; d < data.getDepth(); ++d)
        {
            if constexpr (std::is_same_v<Value, MatrixNd>)
            {
                const std::pair<int, int> ij = tensor_index(d, storage);
                data(index, d) = value(ij.first, ij.second);
            }
            else
            {
                data(index, d) = value[d];
            }
        }
    }
}

template <typename Value, typename Function>
template <typename Data>
void
CartGridPointwiseFunction<Value, Function>::applyPointwise(Data& data,
                                                           const SAMRAI::hier::Box<NDIM>& box,
                                                           const SAMRAI::geom::CartesianPatchGeometry<NDIM>& geometry,
                                                           const double time)
{
    Values::validate_depth(data.getDepth(), d_tensor_storage);
    Value q = Values::make_value(data.getDepth());
    Value result = Values::make_value(data.getDepth());
    const double* const x_lower = geometry.getXLower();
    const double* const dx = geometry.getDx();
    const SAMRAI::hier::Index<NDIM>& index_lower = box.lower();
    const int n_groups = std::is_same_v<Value, double> ? data.getDepth() : 1;
    const int n_axes = Centering<Data>::is_staggered() ? NDIM : 1;
    for (int depth = 0; depth < n_groups; ++depth)
    {
        for (int orientation = 0; orientation < n_axes; ++orientation)
        {
            if (!Centering<Data>::has_axis(data, orientation))
            {
                continue;
            }
            const int axis = Centering<Data>::is_staggered() ? orientation : invalid_index;
            const VectorNd offset = Centering<Data>::offset(orientation);
            for (auto it = Centering<Data>::begin(box, orientation); it; it++)
            {
                const auto& index = it();
                const SAMRAI::hier::Index<NDIM> cartesian_index = Centering<Data>::cartesian_index(index);
                VectorNd x;
                for (int d = 0; d < NDIM; ++d)
                {
                    x[d] = x_lower[d] + dx[d] * (cartesian_index(d) - index_lower(d) + offset[d]);
                }
                if constexpr (s_transforms)
                {
                    Values::load(q, data, index, depth, d_tensor_storage);
                    Values::assign_result(
                        result,
                        std::invoke(d_function, std::as_const(q), std::as_const(x), time, depth, axis),
                        data.getDepth());
                }
                else
                {
                    Values::assign_result(
                        result, std::invoke(d_function, std::as_const(x), time, depth, axis), data.getDepth());
                }
                Values::store(result, data, index, depth, d_tensor_storage);
            }
        }
    }
}

template <typename Value, typename Function>
SAMRAI::tbox::Pointer<CartGridFunction>
make_cart_grid_pointwise_function(std::string object_name, Function&& function)
{
    return new CartGridPointwiseFunction<Value, std::decay_t<Function>>(std::move(object_name),
                                                                        std::forward<Function>(function));
}

template <typename Value, typename Function>
SAMRAI::tbox::Pointer<CartGridFunction>
make_cart_grid_pointwise_function(std::string object_name, Function&& function, const TensorStorage storage)
{
    return new CartGridPointwiseFunction<Value, std::decay_t<Function>>(
        std::move(object_name), std::forward<Function>(function), storage);
}
} // namespace IBTK

#endif // #ifndef included_IBTK_CartGridPointwiseFunction_inl
