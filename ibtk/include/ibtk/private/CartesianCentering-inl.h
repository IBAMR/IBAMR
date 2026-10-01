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

#ifndef included_IBTK_CartesianCentering_inl
#define included_IBTK_CartesianCentering_inl

#include <ibtk/config.h>

#include <ibtk/CartesianCentering.h>

#include <cstdlib>
#include <utility>

namespace IBTK
{
template <>
inline DataCentering
string_to_enum<DataCentering>(const std::string& value)
{
    if (strcasecmp(value.c_str(), "CELL") == 0)
    {
        return DataCentering::CELL;
    }
    if (strcasecmp(value.c_str(), "NODE") == 0)
    {
        return DataCentering::NODE;
    }
    if (strcasecmp(value.c_str(), "SIDE") == 0)
    {
        return DataCentering::SIDE;
    }
    if (strcasecmp(value.c_str(), "FACE") == 0)
    {
        return DataCentering::FACE;
    }
    if (strcasecmp(value.c_str(), "EDGE") == 0)
    {
        return DataCentering::EDGE;
    }
    TBOX_ERROR("Unknown DataCentering: " << value << '\n');
    std::abort();
}

template <>
inline std::string
enum_to_string<DataCentering>(const DataCentering value)
{
    switch (value)
    {
    case DataCentering::CELL:
        return "CELL";
    case DataCentering::NODE:
        return "NODE";
    case DataCentering::SIDE:
        return "SIDE";
    case DataCentering::FACE:
        return "FACE";
    case DataCentering::EDGE:
        return "EDGE";
    default:
        TBOX_ERROR("Unknown DataCentering\n");
        std::abort();
    }
}

template <DataCentering C>
constexpr bool
CartesianCentering<C>::is_staggered()
{
    if constexpr (C == DataCentering::CELL || C == DataCentering::NODE)
    {
        return false;
    }
    else
    {
        static_assert(C == DataCentering::SIDE || C == DataCentering::FACE || C == DataCentering::EDGE,
                      "Unsupported Cartesian centering.");
        return true;
    }
}

template <DataCentering C>
constexpr int
CartesianCentering<C>::num_axes()
{
    return is_staggered() ? NDIM : 1;
}

template <DataCentering C>
template <typename T>
bool
CartesianCentering<C>::has_axis(const Data<T>& data, const int axis)
{
    if constexpr (C == DataCentering::SIDE)
    {
        return data.getDirectionVector()(axis) != 0;
    }
    else
    {
        return true;
    }
}

template <DataCentering C>
typename CartesianCentering<C>::template Data<double>::Iterator
CartesianCentering<C>::begin(const SAMRAI::hier::Box<NDIM>& box, const int axis)
{
    if constexpr (is_staggered())
    {
        return typename Data<double>::Iterator(box, axis);
    }
    else
    {
        return typename Data<double>::Iterator(box);
    }
}

template <DataCentering C>
SAMRAI::hier::Box<NDIM>
CartesianCentering<C>::index_box(const SAMRAI::hier::Box<NDIM>& cell_box, const int axis)
{
    if constexpr (C == DataCentering::CELL)
    {
        return cell_box;
    }
    else if constexpr (C == DataCentering::NODE)
    {
        return SAMRAI::pdat::NodeGeometry<NDIM>::toNodeBox(cell_box);
    }
    else if constexpr (C == DataCentering::SIDE)
    {
        return SAMRAI::pdat::SideGeometry<NDIM>::toSideBox(cell_box, axis);
    }
    else if constexpr (C == DataCentering::FACE)
    {
        return SAMRAI::pdat::FaceGeometry<NDIM>::toFaceBox(cell_box, axis);
    }
    else
    {
        static_assert(C == DataCentering::EDGE, "Unsupported Cartesian centering.");
        return SAMRAI::pdat::EdgeGeometry<NDIM>::toEdgeBox(cell_box, axis);
    }
}

template <DataCentering C>
template <typename T>
SAMRAI::pdat::ArrayData<NDIM, T>&
CartesianCentering<C>::array_data(Data<T>& data, const int axis)
{
    if constexpr (is_staggered())
    {
        return data.getArrayData(axis);
    }
    else
    {
        return data.getArrayData();
    }
}

template <DataCentering C>
template <typename T>
const SAMRAI::pdat::ArrayData<NDIM, T>&
CartesianCentering<C>::array_data(const Data<T>& data, const int axis)
{
    if constexpr (is_staggered())
    {
        return data.getArrayData(axis);
    }
    else
    {
        return data.getArrayData();
    }
}

template <DataCentering C>
VectorNd
CartesianCentering<C>::offset(const int axis)
{
    VectorNd result;
    if constexpr (C == DataCentering::CELL)
    {
        result.setConstant(0.5);
    }
    else if constexpr (C == DataCentering::NODE)
    {
        result.setZero();
    }
    else if constexpr (C == DataCentering::SIDE || C == DataCentering::FACE)
    {
        result.setConstant(0.5);
        result[axis] = 0.0;
    }
    else
    {
        static_assert(C == DataCentering::EDGE, "Unsupported Cartesian centering.");
        result.setZero();
        result[axis] = 0.5;
    }
    return result;
}

template <DataCentering C>
SAMRAI::hier::Index<NDIM>
CartesianCentering<C>::cartesian_index(const Index& index)
{
    if constexpr (C == DataCentering::FACE)
    {
        // Undo FaceIndex's coordinate permutation without shifting the normal index.
        return index.toCell(1);
    }
    else
    {
        return index;
    }
}

template <typename T>
std::optional<DataCentering>
find_data_centering(const SAMRAI::hier::PatchDataFactory<NDIM>& factory)
{
    if (dynamic_cast<const typename CartesianCentering<DataCentering::CELL>::template Factory<T>*>(&factory))
    {
        return DataCentering::CELL;
    }
    if (dynamic_cast<const typename CartesianCentering<DataCentering::NODE>::template Factory<T>*>(&factory))
    {
        return DataCentering::NODE;
    }
    if (dynamic_cast<const typename CartesianCentering<DataCentering::SIDE>::template Factory<T>*>(&factory))
    {
        return DataCentering::SIDE;
    }
    if (dynamic_cast<const typename CartesianCentering<DataCentering::FACE>::template Factory<T>*>(&factory))
    {
        return DataCentering::FACE;
    }
    if (dynamic_cast<const typename CartesianCentering<DataCentering::EDGE>::template Factory<T>*>(&factory))
    {
        return DataCentering::EDGE;
    }
    return std::nullopt;
}

template <typename T>
DataCentering
get_data_centering(const SAMRAI::hier::PatchDataFactory<NDIM>& factory)
{
    const std::optional<DataCentering> centering = find_data_centering<T>(factory);
    if (!centering)
    {
        TBOX_ERROR("get_data_centering: unsupported patch data factory for the requested scalar type\n");
    }
    return *centering;
}

template <typename Function>
decltype(auto)
dispatch_data_centering(const DataCentering centering, Function&& function)
{
    switch (centering)
    {
    case DataCentering::CELL:
        return std::forward<Function>(function).template operator()<DataCentering::CELL>();
    case DataCentering::NODE:
        return std::forward<Function>(function).template operator()<DataCentering::NODE>();
    case DataCentering::SIDE:
        return std::forward<Function>(function).template operator()<DataCentering::SIDE>();
    case DataCentering::FACE:
        return std::forward<Function>(function).template operator()<DataCentering::FACE>();
    case DataCentering::EDGE:
        return std::forward<Function>(function).template operator()<DataCentering::EDGE>();
    default:
        TBOX_ERROR("dispatch_data_centering: unsupported data centering\n");
        std::abort();
    }
}
} // namespace IBTK

#endif // #ifndef included_IBTK_CartesianCentering_inl
