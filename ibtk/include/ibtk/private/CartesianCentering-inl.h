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
    return DataCentering::UNKNOWN;
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
        return "UNKNOWN";
    }
}

template <DataCentering C>
constexpr bool
CartesianCentering<C>::is_oriented()
{
    return C == DataCentering::SIDE || C == DataCentering::FACE || C == DataCentering::EDGE;
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
    if constexpr (is_oriented())
    {
        return typename Data<double>::Iterator(box, axis);
    }
    else
    {
        return typename Data<double>::Iterator(box);
    }
}

template <DataCentering C>
VectorNd
CartesianCentering<C>::offset(const int axis)
{
    VectorNd result;
    if constexpr (C == DataCentering::NODE)
    {
        result.setZero();
    }
    else if constexpr (C == DataCentering::EDGE)
    {
        result.setZero();
        result[axis] = 0.5;
    }
    else
    {
        result.setConstant(0.5);
        if constexpr (is_oriented())
        {
            result[axis] = 0.0;
        }
    }
    return result;
}

template <DataCentering C>
SAMRAI::hier::Index<NDIM>
CartesianCentering<C>::cartesian_index(const Index& index)
{
    if constexpr (C == DataCentering::FACE)
    {
        return index.toCell(1);
    }
    else
    {
        return index;
    }
}

template <typename T>
DataCentering
get_data_centering(const SAMRAI::hier::PatchDataFactory<NDIM>& factory)
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
    return DataCentering::UNKNOWN;
}

template <typename Function>
decltype(auto)
dispatch_data_centering(const DataCentering centering, Function&& function)
{
    switch (centering)
    {
    case DataCentering::CELL:
        return std::forward<Function>(function).template operator()<CartesianCentering<DataCentering::CELL>>();
    case DataCentering::NODE:
        return std::forward<Function>(function).template operator()<CartesianCentering<DataCentering::NODE>>();
    case DataCentering::SIDE:
        return std::forward<Function>(function).template operator()<CartesianCentering<DataCentering::SIDE>>();
    case DataCentering::FACE:
        return std::forward<Function>(function).template operator()<CartesianCentering<DataCentering::FACE>>();
    case DataCentering::EDGE:
        return std::forward<Function>(function).template operator()<CartesianCentering<DataCentering::EDGE>>();
    default:
        TBOX_ERROR("dispatch_data_centering: unsupported data centering\n");
        std::abort();
    }
}
} // namespace IBTK

#endif // #ifndef included_IBTK_CartesianCentering_inl
