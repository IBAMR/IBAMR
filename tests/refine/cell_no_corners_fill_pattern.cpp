// ---------------------------------------------------------------------
//
// Copyright (c) 2019 - 2025 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files

#include <CellGeometry.h>
#include <SAMRAI_config.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

// Set up application namespace declarations
#include <ibtk/app_namespaces.h>

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // this test only works in serial
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cell_no_corners_pattern.log");

    const int N = 10;
    Box<NDIM> dst_box, src_box;
    for (int d = 0; d < NDIM; ++d)
    {
        dst_box.lower(d) = 0;
        dst_box.upper(d) = N - 1;
    }
    const int gcw = 1;
    CellNoCornersFillPattern no_corners_pattern(/*stencil_width*/ gcw, /*overwrite_interior*/ false);
    CellNoCornersFillPattern no_corners_pattern_interior(/*stencil_width*/ gcw, /*overwrite_interior*/ true);

    const int dst_level_num = 1, src_level_num = 1;
    no_corners_pattern.setTargetPatchLevelNumber(dst_level_num);
    no_corners_pattern_interior.setTargetPatchLevelNumber(dst_level_num);

    for (int s0 = -1; s0 <= +1; ++s0)
    {
        src_box.lower(0) = dst_box.lower(0) + s0 * N;
        src_box.upper(0) = dst_box.upper(0) + s0 * N;
        for (int s1 = -1; s1 <= +1; ++s1)
        {
            src_box.lower(1) = dst_box.lower(1) + s1 * N;
            src_box.upper(1) = dst_box.upper(1) + s1 * N;
#if (NDIM == 3)
            for (int s2 = -1; s2 <= +1; ++s2)
            {
                src_box.lower(2) = dst_box.lower(2) + s2 * N;
                src_box.upper(2) = dst_box.upper(2) + s2 * N;
#endif
                tbox::plog << "**************************************************\n";
                tbox::plog << "dst_box: " << dst_box << "\n";
                tbox::plog << "src_box: " << src_box << "\n";

                CellGeometry<NDIM> dst_geometry(dst_box, gcw);
                CellGeometry<NDIM> src_geometry(src_box, gcw);
                Box<NDIM> src_mask = src_box;

                Pointer<CellOverlap<NDIM> > overlap;

                tbox::plog << "default overlap, overwrite_interior == false:\n";
                overlap = dst_geometry.calculateOverlap(dst_geometry,
                                                        src_geometry,
                                                        src_mask,
                                                        /*overwrite_interior*/ false,
                                                        /*src_offset*/ IntVector<NDIM>(0),
                                                        /*retry*/ false);
                overlap->print(tbox::plog);

                tbox::plog << "default overlap, overwrite_interior == true:\n";
                overlap = dst_geometry.calculateOverlap(dst_geometry,
                                                        src_geometry,
                                                        src_mask,
                                                        /*overwrite_interior*/ true,
                                                        /*src_offset*/ IntVector<NDIM>(0),
                                                        /*retry*/ false);
                overlap->print(tbox::plog);

                tbox::plog << "dst_level_num == src_level_num:\n";

                tbox::plog << "no corners pattern, overwrite_interior == false:\n";
                overlap = no_corners_pattern.calculateOverlapOnLevel(dst_geometry,
                                                                     src_geometry,
                                                                     dst_box,
                                                                     src_mask,
                                                                     /*overwrite_interior*/ false,
                                                                     /*src_offset*/ IntVector<NDIM>(0),
                                                                     dst_level_num,
                                                                     src_level_num);
                overlap->print(tbox::plog);

                tbox::plog << "no corners pattern, overwrite_interior == true:\n";
                overlap = no_corners_pattern_interior.calculateOverlapOnLevel(dst_geometry,
                                                                              src_geometry,
                                                                              dst_box,
                                                                              src_mask,
                                                                              /*overwrite_interior*/ true,
                                                                              /*src_offset*/ IntVector<NDIM>(0),
                                                                              dst_level_num,
                                                                              src_level_num);
                overlap->print(tbox::plog);

                tbox::plog << "dst_level_num != target_level_num:\n";

                tbox::plog << "no corners pattern, overwrite_interior == false:\n";
                overlap = no_corners_pattern.calculateOverlapOnLevel(dst_geometry,
                                                                     src_geometry,
                                                                     dst_box,
                                                                     src_mask,
                                                                     /*overwrite_interior*/ false,
                                                                     /*src_offset*/ IntVector<NDIM>(0),
                                                                     dst_level_num - 1,
                                                                     src_level_num - 1);
                overlap->print(tbox::plog);

                tbox::plog << "no corners pattern, overwrite_interior == true:\n";
                overlap = no_corners_pattern_interior.calculateOverlapOnLevel(dst_geometry,
                                                                              src_geometry,
                                                                              dst_box,
                                                                              src_mask,
                                                                              /*overwrite_interior*/ true,
                                                                              /*src_offset*/ IntVector<NDIM>(0),
                                                                              dst_level_num - 1,
                                                                              src_level_num - 1);
                overlap->print(tbox::plog);
#if (NDIM == 3)
            }
#endif
        }
    }
}
