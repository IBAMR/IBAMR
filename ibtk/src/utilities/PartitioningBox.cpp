// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2020 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

/////////////////////////////// INCLUDES /////////////////////////////////////
#include "ibtk/samrai_compatibility_names.h"
#include <ibtk/PartitioningBox.h>
#include <ibtk/ibtk_utilities.h>

#include "SAMRAICartesianPatchGeometry.h"
#include "SAMRAIIntVector.h"
#include "SAMRAIPatch.h"
#include "SAMRAIPatchHierarchy.h"
#include "SAMRAIPatchLevel.h"
#include "SAMRAIPointer.h"
#include "SAMRAIUtilities.h"

#include <algorithm>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <ibtk/namespaces.h> // IWYU pragma: keep

/////////////////////////////// CLASS DEFINITION /////////////////////////////
namespace IBTK
{
PartitioningBox::PartitioningBox()
{
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
    {
        d_bounding_points.first[dim_n] = std::numeric_limits<double>::infinity();
        d_bounding_points.second[dim_n] = std::numeric_limits<double>::infinity();
    }
} // PartitioningBox

PartitioningBox::PartitioningBox(const Point& bottom_point, const Point& top_point)
    : d_bounding_points(bottom_point, top_point)
{
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
    {
        // Do not permit negative volume boxes
        TBOX_ASSERT(bottom()[dim_n] <= top()[dim_n]);
    }
} // PartitioningBox

PartitioningBox::PartitioningBox(const SAMRAICartesianPatchGeometry& patch)
{
    std::copy(patch.getXLower(), patch.getXLower() + NDIM, d_bounding_points.first.data());
    std::copy(patch.getXUpper(), patch.getXUpper() + NDIM, d_bounding_points.second.data());
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
    {
        // Do not permit negative volume boxes
        TBOX_ASSERT(bottom()[dim_n] <= top()[dim_n]);
    }
} // PartitioningBox

double
PartitioningBox::volume() const
{
    // The bounds may be infinite: check for point equality first to avoid
    // evaluating oo - oo
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
    {
        if (d_bounding_points.first[dim_n] == d_bounding_points.second[dim_n]) return 0.0;
    }
    double vol = 1.0;
    for (unsigned int dim_n = 0; dim_n < NDIM; ++dim_n)
    {
        vol *= d_bounding_points.second[dim_n] - d_bounding_points.first[dim_n];
    }
    return vol;
} // volume

PartitioningBoxes::PartitioningBoxes(const SAMRAIPatchHierarchy& hierarchy)
{
    std::vector<IBTK::PartitioningBox> boxes;
    const int finest_level = hierarchy.getFinestLevelNumber();
    SAMRAIPointer<SAMRAIPatchLevel> level = hierarchy.getPatchLevel(finest_level);
    for (SAMRAIPatchLevel::Iterator p(level); p; p++)
    {
        const SAMRAIPatch& patch = *level->getPatch(p());
        SAMRAIPointer<SAMRAICartesianPatchGeometry> patch_geometry = patch.getPatchGeometry();
        boxes.emplace_back(*patch_geometry);
    }

    *this = PartitioningBoxes(boxes.begin(), boxes.end());
}

} // namespace IBTK
//////////////////////////////////////////////////////////////////////////////
