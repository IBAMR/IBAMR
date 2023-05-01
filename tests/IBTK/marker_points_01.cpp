// ---------------------------------------------------------------------
//
// Copyright (c) 2023 - 2023 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/MarkerPatchHierarchy.h>

#include <fstream>

#include <ibtk/app_namespaces.h>

#include "../tests.h"

void
test(MarkerPatch& marker_patch, std::ostream& output)
{
    output << "before pruning:\n";
    for (unsigned int i = 0; i < marker_patch.size(); ++i)
    {
        const auto marker = marker_patch[i];
        output << "  index = " << std::get<0>(marker) << " X = ";
        for (unsigned int d = 0; d < NDIM - 1; ++d)
        {
            output << std::get<1>(marker)[d] << ", ";
        }
        output << std::get<1>(marker)[NDIM - 1] << " V = ";
        for (unsigned int d = 0; d < NDIM - 1; ++d)
        {
            output << std::get<2>(marker)[d] << ", ";
        }
        output << std::get<2>(marker)[NDIM - 1] << " contains = " << marker_patch.contains(std::get<1>(marker)) << '\n';
    }

    auto result = marker_patch.prune();
    output << "after pruning:\n";
    for (unsigned int i = 0; i < marker_patch.size(); ++i)
    {
        const auto marker = marker_patch[i];
        output << "  index = " << std::get<0>(marker) << " X = ";
        for (unsigned int d = 0; d < NDIM - 1; ++d)
        {
            output << std::get<1>(marker)[d] << ", ";
        }
        output << std::get<1>(marker)[NDIM - 1] << " V = ";
        for (unsigned int d = 0; d < NDIM - 1; ++d)
        {
            output << std::get<2>(marker)[d] << ", ";
        }
        output << std::get<2>(marker)[NDIM - 1] << std::endl;
    }

    output << "pruned points:\n";
    for (const auto& index : std::get<0>(result))
    {
        output << "  " << index << std::endl;
    }
}

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    const int N = input_db->getInteger("N");

    IntVector<NDIM> ratio(2);
    // Set up ctor arguments:
    IntVector<NDIM> lower(0);
    IntVector<NDIM> upper(ratio(0) * N);
    Box<NDIM> patch_box(lower, upper);
    std::vector<Box<NDIM> > nonoverlapping_patch_boxes;
    {
        IntVector<NDIM> mid(ratio(0) * N / 2);
        nonoverlapping_patch_boxes.emplace_back(lower, mid);
        nonoverlapping_patch_boxes.emplace_back(mid, upper);
    }
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry =
        new CartesianGridGeometry<NDIM>("CartesianGridGeometry", input_db->getDatabase("CartesianGeometry"), true);

    std::ofstream output("output");
    {
        MarkerPatch marker_patch(patch_box, nonoverlapping_patch_boxes, grid_geometry, ratio);
        marker_patch.insert(9, IBTK::Point(1.0, 0.25), IBTK::Vector(9, 8));
        marker_patch.insert(5, IBTK::Point(0.1, 0.2), IBTK::Vector(5, 4));
        marker_patch.insert(3, IBTK::Point(0.6, 0.7), IBTK::Vector(3, 2));
        marker_patch.insert(4, IBTK::Point(0.2, 0.8), IBTK::Vector(4, 3));
        marker_patch.insert(8, IBTK::Point(0.9, 0.9), IBTK::Vector(8, 7));
        test(marker_patch, output);
    }
}
