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

#ifndef included_IBTK_CartesianCentering
#define included_IBTK_CartesianCentering

#include <ibtk/config.h>

#include <ibtk/ibtk_enums.h>
#include <ibtk/ibtk_utilities.h>

#include <ArrayData.h>
#include <Box.h>
#include <CellData.h>
#include <CellDataFactory.h>
#include <CellIndex.h>
#include <EdgeData.h>
#include <EdgeDataFactory.h>
#include <EdgeGeometry.h>
#include <EdgeIndex.h>
#include <FaceData.h>
#include <FaceDataFactory.h>
#include <FaceGeometry.h>
#include <FaceIndex.h>
#include <NodeData.h>
#include <NodeDataFactory.h>
#include <NodeGeometry.h>
#include <NodeIndex.h>
#include <PatchCellDataBasicOps.h>
#include <PatchCellDataNormOpsReal.h>
#include <PatchEdgeDataBasicOps.h>
#include <PatchEdgeDataNormOpsReal.h>
#include <PatchFaceDataBasicOps.h>
#include <PatchFaceDataNormOpsReal.h>
#include <PatchNodeDataBasicOps.h>
#include <PatchNodeDataNormOpsReal.h>
#include <PatchSideDataBasicOps.h>
#include <PatchSideDataNormOpsReal.h>
#include <SideData.h>
#include <SideDataFactory.h>
#include <SideGeometry.h>
#include <SideIndex.h>

#include <optional>
#include <string>
#include <tuple>

namespace IBTK
{
/*! \brief Cartesian patch-data centerings. */
enum class DataCentering
{
    CELL = 0,
    NODE = 1,
    SIDE = 2,
    FACE = 3,
    EDGE = 4
};

template <>
DataCentering string_to_enum<DataCentering>(const std::string& value);

template <>
std::string enum_to_string<DataCentering>(DataCentering value);

/*!
 * \brief Types and geometry for one Cartesian data centering.
 *
 * The enum selects the SAMRAI data, factory, index, patch arithmetic, and patch norm types.
 * Scalar type and data depth are independent of centering; depth is the number
 * of values at each grid location and is obtained from the data object.
 * Side, face, and edge data are staggered, with separate arrays for each
 * coordinate direction. The axis specifies the normal direction for side and
 * face data and the tangent direction for edge data. It must be in [0, NDIM);
 * cell and node operations ignore it.
 */
template <DataCentering C>
struct CartesianCentering
{
    static_assert(C >= DataCentering::CELL && C <= DataCentering::EDGE, "Unsupported Cartesian centering.");

    template <typename T>
    using Data = std::tuple_element_t<static_cast<int>(C),
                                      std::tuple<SAMRAI::pdat::CellData<NDIM, T>,
                                                 SAMRAI::pdat::NodeData<NDIM, T>,
                                                 SAMRAI::pdat::SideData<NDIM, T>,
                                                 SAMRAI::pdat::FaceData<NDIM, T>,
                                                 SAMRAI::pdat::EdgeData<NDIM, T>>>;

    template <typename T>
    using Factory = std::tuple_element_t<static_cast<int>(C),
                                         std::tuple<SAMRAI::pdat::CellDataFactory<NDIM, T>,
                                                    SAMRAI::pdat::NodeDataFactory<NDIM, T>,
                                                    SAMRAI::pdat::SideDataFactory<NDIM, T>,
                                                    SAMRAI::pdat::FaceDataFactory<NDIM, T>,
                                                    SAMRAI::pdat::EdgeDataFactory<NDIM, T>>>;

    template <typename T>
    using PatchOps = std::tuple_element_t<static_cast<int>(C),
                                          std::tuple<SAMRAI::math::PatchCellDataBasicOps<NDIM, T>,
                                                     SAMRAI::math::PatchNodeDataBasicOps<NDIM, T>,
                                                     SAMRAI::math::PatchSideDataBasicOps<NDIM, T>,
                                                     SAMRAI::math::PatchFaceDataBasicOps<NDIM, T>,
                                                     SAMRAI::math::PatchEdgeDataBasicOps<NDIM, T>>>;

    template <typename T>
    using PatchNormOps = std::tuple_element_t<static_cast<int>(C),
                                              std::tuple<SAMRAI::math::PatchCellDataNormOpsReal<NDIM, T>,
                                                         SAMRAI::math::PatchNodeDataNormOpsReal<NDIM, T>,
                                                         SAMRAI::math::PatchSideDataNormOpsReal<NDIM, T>,
                                                         SAMRAI::math::PatchFaceDataNormOpsReal<NDIM, T>,
                                                         SAMRAI::math::PatchEdgeDataNormOpsReal<NDIM, T>>>;

    using Index = std::tuple_element_t<static_cast<int>(C),
                                       std::tuple<SAMRAI::pdat::CellIndex<NDIM>,
                                                  SAMRAI::pdat::NodeIndex<NDIM>,
                                                  SAMRAI::pdat::SideIndex<NDIM>,
                                                  SAMRAI::pdat::FaceIndex<NDIM>,
                                                  SAMRAI::pdat::EdgeIndex<NDIM>>>;

    /*! \brief Whether this centering is side-, face-, or edge-centered. */
    static constexpr bool is_staggered();

    /*! \brief Return the number of coordinate directions with separate arrays: NDIM if staggered, 1 otherwise. */
    static constexpr int num_axes();

    /*!
     * \brief Whether the data allocate the requested coordinate direction.
     *
     * Side data use their direction vector. Face and edge data allocate every
     * direction. Cell and node data ignore axis and always return true.
     */
    template <typename T>
    static bool has_axis(const Data<T>& data, int axis);

    /*! \brief Iterate over the given cell box in the requested coordinate direction. */
    static typename Data<double>::Iterator begin(const SAMRAI::hier::Box<NDIM>& box, int axis);

    /*!
     * \brief Return the box of this centering's indices associated with a cell box.
     *
     * The box is in the index space of the corresponding array of the data object;
     * face boxes use FaceData's permuted index order, as in cartesian_index().
     */
    static SAMRAI::hier::Box<NDIM> index_box(const SAMRAI::hier::Box<NDIM>& cell_box, int axis);

    /*! \brief Return the array storing the data for the requested coordinate direction. */
    template <typename T>
    static SAMRAI::pdat::ArrayData<NDIM, T>& array_data(Data<T>& data, int axis);

    /*! \brief Return the read-only array storing the data for the requested coordinate direction. */
    template <typename T>
    static const SAMRAI::pdat::ArrayData<NDIM, T>& array_data(const Data<T>& data, int axis);

    /*! \brief Return the offset from a cell's lower corner in units of its widths. */
    static VectorNd offset(int axis);

    /*! \brief Return the cell, node, side, face, or edge index coordinates in Cartesian order. */
    static SAMRAI::hier::Index<NDIM> cartesian_index(const Index& index);
};

/*!
 * \brief Identify a Cartesian factory allocating scalar type T, or return an empty optional if there is none.
 *
 * Unlike get_data_centering(), an unsupported centering or scalar type is not an error.
 */
template <typename T>
std::optional<DataCentering> find_data_centering(const SAMRAI::hier::PatchDataFactory<NDIM>& factory);

/*!
 * \brief Identify a Cartesian factory allocating scalar type T.
 *
 * An unsupported centering or scalar type is a fatal error; use
 * find_data_centering() for a non-fatal lookup. A recognized centering does not
 * establish the type currently stored at a patch-data index. Callers must
 * enforce that separate precondition.
 */
template <typename T>
DataCentering get_data_centering(const SAMRAI::hier::PatchDataFactory<NDIM>& factory);

/*!
 * \brief Invoke function.template operator()<C>() once for the selected DataCentering.
 *
 * All centerings must produce the same return type. An invalid enum is a fatal error.
 */
template <typename Function>
decltype(auto) dispatch_data_centering(DataCentering centering, Function&& function);
} // namespace IBTK

#include <ibtk/private/CartesianCentering-inl.h>

#endif // #ifndef included_IBTK_CartesianCentering
