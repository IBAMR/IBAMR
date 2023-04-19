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

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <fstream>

#include <ibtk/app_namespaces.h>

#include "../tests.h"

void
test(const MarkerPatch& marker_patch, std::ostream& output)
{
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
}

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv);
    // Skip warnings:
    Logger::getInstance()->setWarning(false);

    const auto rank = IBTK_MPI::getRank();
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM> > error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", NULL, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM> > load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

    // Initialize the AMR patch hierarchy.
    gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
    int tag_buffer = 1;
    int level_number = 0;
    bool done = false;
    while (!done && (gridding_algorithm->levelCanBeRefined(level_number)))
    {
        gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
        done = !patch_hierarchy->finerLevelExists(level_number);
        ++level_number;
    }

    // Set up marker points.
    EigenAlignedVector<IBTK::Point> positions;
    EigenAlignedVector<IBTK::Vector> velocities;
    for (unsigned int i = 0; i < 10; ++i)
    {
        for (unsigned int j = 0; j < 10; ++j)
        {
            positions.emplace_back(double(i) / 10.0, double(j) / 10.0);
            velocities.emplace_back(double(i) * 10 + j, 0.0);
        }
    }
    MarkerPatchHierarchy marker_points("MarkerPoints", patch_hierarchy, positions, velocities);

    std::ostringstream out;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        int local_patch_num = 0;
        Pointer<PatchLevel<NDIM> > current_level = patch_hierarchy->getPatchLevel(ln);
        for (int p = 0; p < current_level->getNumberOfPatches(); ++p)
        {
            if (rank == current_level->getMappingForPatch(p))
            {
                out << "level = " << ln << " patch = " << current_level->getPatch(p)->getBox() << std::endl;
                test(marker_points.getMarkerPatch(ln, local_patch_num), out);

                ++local_patch_num;
            }
        }
    }

    print_strings_on_plog_0(out.str());
}
