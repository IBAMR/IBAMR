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

#ifndef included_IBTK_data_kind
#define included_IBTK_data_kind

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibtk/config.h>

#include <string>

namespace IBTK
{
namespace detail
{
/*!
 * \brief The centering of double precision patch data.
 */
enum class DataCentering
{
    CELL,
    NODE,
    SIDE,
    EDGE,
    FACE,
    OTHER
};

/*!
 * \brief The centering and depth of double precision patch data.
 */
struct DataKind
{
    DataCentering centering = DataCentering::OTHER;
    int depth = 0;

    bool operator==(const DataKind& other) const;

    bool operator!=(const DataKind& other) const;
};

/*!
 * \brief Return the centering and depth of the patch data with index data_idx, as the patch descriptor has them.
 * Data of another type have the centering OTHER and depth 0.
 */
DataKind get_data_kind(int data_idx);

/*!
 * \brief Return a description of the data kind for error messages, such as "cell-centered data of depth 1".
 */
std::string data_kind_to_string(const DataKind& kind);
} // namespace detail
} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////

#endif // #ifndef included_IBTK_data_kind
