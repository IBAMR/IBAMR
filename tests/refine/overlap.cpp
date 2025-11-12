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
#include <EdgeGeometry.h>
#include <FaceGeometry.h>
#include <NodeGeometry.h>
#include <OuteredgeGeometry.h>
#include <OuterfaceGeometry.h>
#include <OuternodeGeometry.h>
#include <OutersideGeometry.h>
#include <SAMRAI_config.h>
#include <SideGeometry.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
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

                Box<NDIM> src_mask = src_box;

                CellGeometry<NDIM> cell_dst_geometry(dst_box, gcw);
                CellGeometry<NDIM> cell_src_geometry(src_box, gcw);

                Pointer<CellOverlap<NDIM> > cell_overlap;

                tbox::plog << "cell overlap, overwrite_interior == false:\n";
                cell_overlap = cell_dst_geometry.calculateOverlap(cell_dst_geometry,
                                                                  cell_src_geometry,
                                                                  src_mask,
                                                                  /*overwrite_interior*/ false,
                                                                  /*src_offset*/ IntVector<NDIM>(0),
                                                                  /*retry*/ false);
                cell_overlap->getDestinationBoxList().print(tbox::plog);

                tbox::plog << "cell overlap, overwrite_interior == true:\n";
                cell_overlap = cell_dst_geometry.calculateOverlap(cell_dst_geometry,
                                                                  cell_src_geometry,
                                                                  src_mask,
                                                                  /*overwrite_interior*/ true,
                                                                  /*src_offset*/ IntVector<NDIM>(0),
                                                                  /*retry*/ false);
                cell_overlap->getDestinationBoxList().print(tbox::plog);

                EdgeGeometry<NDIM> edge_dst_geometry(dst_box, gcw);
                EdgeGeometry<NDIM> edge_src_geometry(src_box, gcw);

                Pointer<EdgeOverlap<NDIM> > edge_overlap;

                tbox::plog << "edge overlap, overwrite_interior == false:\n";
                edge_overlap = edge_dst_geometry.calculateOverlap(edge_dst_geometry,
                                                                  edge_src_geometry,
                                                                  src_mask,
                                                                  /*overwrite_interior*/ false,
                                                                  /*src_offset*/ IntVector<NDIM>(0),
                                                                  /*retry*/ false);
                for (int axis = 0; axis < NDIM; ++axis) edge_overlap->getDestinationBoxList(axis).print(tbox::plog);

                tbox::plog << "edge overlap, overwrite_interior == true:\n";
                edge_overlap = edge_dst_geometry.calculateOverlap(edge_dst_geometry,
                                                                  edge_src_geometry,
                                                                  src_mask,
                                                                  /*overwrite_interior*/ true,
                                                                  /*src_offset*/ IntVector<NDIM>(0),
                                                                  /*retry*/ false);
                for (int axis = 0; axis < NDIM; ++axis) edge_overlap->getDestinationBoxList(axis).print(tbox::plog);

                FaceGeometry<NDIM> face_dst_geometry(dst_box, gcw);
                FaceGeometry<NDIM> face_src_geometry(src_box, gcw);

                Pointer<FaceOverlap<NDIM> > face_overlap;

                tbox::plog << "face overlap, overwrite_interior == false:\n";
                face_overlap = face_dst_geometry.calculateOverlap(face_dst_geometry,
                                                                  face_src_geometry,
                                                                  src_mask,
                                                                  /*overwrite_interior*/ false,
                                                                  /*src_offset*/ IntVector<NDIM>(0),
                                                                  /*retry*/ false);
                for (int axis = 0; axis < NDIM; ++axis) face_overlap->getDestinationBoxList(axis).print(tbox::plog);

                tbox::plog << "face overlap, overwrite_interior == true:\n";
                face_overlap = face_dst_geometry.calculateOverlap(face_dst_geometry,
                                                                  face_src_geometry,
                                                                  src_mask,
                                                                  /*overwrite_interior*/ true,
                                                                  /*src_offset*/ IntVector<NDIM>(0),
                                                                  /*retry*/ false);
                for (int axis = 0; axis < NDIM; ++axis) face_overlap->getDestinationBoxList(axis).print(tbox::plog);

                NodeGeometry<NDIM> node_dst_geometry(dst_box, gcw);
                NodeGeometry<NDIM> node_src_geometry(src_box, gcw);

                Pointer<NodeOverlap<NDIM> > node_overlap;

                tbox::plog << "node overlap, overwrite_interior == false:\n";
                node_overlap = node_dst_geometry.calculateOverlap(node_dst_geometry,
                                                                  node_src_geometry,
                                                                  src_mask,
                                                                  /*overwrite_interior*/ false,
                                                                  /*src_offset*/ IntVector<NDIM>(0),
                                                                  /*retry*/ false);
                node_overlap->getDestinationBoxList().print(tbox::plog);

                tbox::plog << "node overlap, overwrite_interior == true:\n";
                node_overlap = node_dst_geometry.calculateOverlap(node_dst_geometry,
                                                                  node_src_geometry,
                                                                  src_mask,
                                                                  /*overwrite_interior*/ true,
                                                                  /*src_offset*/ IntVector<NDIM>(0),
                                                                  /*retry*/ false);
                node_overlap->getDestinationBoxList().print(tbox::plog);

                SideGeometry<NDIM> side_dst_geometry(dst_box, gcw, 1);
                SideGeometry<NDIM> side_src_geometry(src_box, gcw, 1);

                Pointer<SideOverlap<NDIM> > side_overlap;

                tbox::plog << "side overlap, overwrite_interior == false:\n";
                side_overlap = side_dst_geometry.calculateOverlap(side_dst_geometry,
                                                                  side_src_geometry,
                                                                  src_mask,
                                                                  /*overwrite_interior*/ false,
                                                                  /*src_offset*/ IntVector<NDIM>(0),
                                                                  /*retry*/ false);
                for (int axis = 0; axis < NDIM; ++axis) side_overlap->getDestinationBoxList(axis).print(tbox::plog);

                tbox::plog << "side overlap, overwrite_interior == true:\n";
                side_overlap = side_dst_geometry.calculateOverlap(side_dst_geometry,
                                                                  side_src_geometry,
                                                                  src_mask,
                                                                  /*overwrite_interior*/ true,
                                                                  /*src_offset*/ IntVector<NDIM>(0),
                                                                  /*retry*/ false);
                for (int axis = 0; axis < NDIM; ++axis) side_overlap->getDestinationBoxList(axis).print(tbox::plog);

                OuteredgeGeometry<NDIM> oedge_dst_geometry(dst_box, gcw);
                OuteredgeGeometry<NDIM> oedge_src_geometry(dst_box, gcw);

                Pointer<EdgeOverlap<NDIM> > oedge_overlap;

                tbox::plog << "outer edge overlap, overwrite_interior == false:\n";
                oedge_overlap = oedge_dst_geometry.calculateOverlap(oedge_dst_geometry,
                                                                    oedge_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ false,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                for (int axis = 0; axis < NDIM; ++axis) oedge_overlap->getDestinationBoxList(axis).print(tbox::plog);

                oedge_overlap = oedge_dst_geometry.calculateOverlap(edge_dst_geometry,
                                                                    oedge_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ false,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                for (int axis = 0; axis < NDIM; ++axis) oedge_overlap->getDestinationBoxList(axis).print(tbox::plog);

                tbox::plog << "outer edge overlap, overwrite_interior == true:\n";
                oedge_overlap = oedge_dst_geometry.calculateOverlap(oedge_dst_geometry,
                                                                    oedge_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ true,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                for (int axis = 0; axis < NDIM; ++axis) oedge_overlap->getDestinationBoxList(axis).print(tbox::plog);

                oedge_overlap = oedge_dst_geometry.calculateOverlap(edge_dst_geometry,
                                                                    oedge_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ true,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                for (int axis = 0; axis < NDIM; ++axis) oedge_overlap->getDestinationBoxList(axis).print(tbox::plog);

                OuterfaceGeometry<NDIM> oface_dst_geometry(dst_box, gcw);
                OuterfaceGeometry<NDIM> oface_src_geometry(dst_box, gcw);

                Pointer<FaceOverlap<NDIM> > oface_overlap;

                tbox::plog << "outer face overlap, overwrite_interior == false:\n";
                oface_overlap = oface_dst_geometry.calculateOverlap(face_dst_geometry,
                                                                    oface_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ false,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                for (int axis = 0; axis < NDIM; ++axis) oface_overlap->getDestinationBoxList(axis).print(tbox::plog);

                tbox::plog << "outer face overlap, overwrite_interior == true:\n";
                oface_overlap = oface_dst_geometry.calculateOverlap(face_dst_geometry,
                                                                    oface_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ true,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                for (int axis = 0; axis < NDIM; ++axis) oface_overlap->getDestinationBoxList(axis).print(tbox::plog);

                OuternodeGeometry<NDIM> onode_dst_geometry(dst_box, gcw);
                OuternodeGeometry<NDIM> onode_src_geometry(dst_box, gcw);

                Pointer<NodeOverlap<NDIM> > onode_overlap;

                tbox::plog << "outer node overlap, overwrite_interior == false:\n";
                onode_overlap = onode_dst_geometry.calculateOverlap(onode_dst_geometry,
                                                                    onode_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ false,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                onode_overlap->getDestinationBoxList().print(tbox::plog);

                onode_overlap = onode_dst_geometry.calculateOverlap(node_dst_geometry,
                                                                    onode_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ false,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                onode_overlap->getDestinationBoxList().print(tbox::plog);

                tbox::plog << "outer node overlap, overwrite_interior == true:\n";
                onode_overlap = onode_dst_geometry.calculateOverlap(onode_dst_geometry,
                                                                    onode_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ true,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                onode_overlap->getDestinationBoxList().print(tbox::plog);

                onode_overlap = onode_dst_geometry.calculateOverlap(node_dst_geometry,
                                                                    onode_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ true,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                onode_overlap->getDestinationBoxList().print(tbox::plog);

                OutersideGeometry<NDIM> oside_dst_geometry(dst_box, gcw);
                OutersideGeometry<NDIM> oside_src_geometry(dst_box, gcw);

                Pointer<SideOverlap<NDIM> > oside_overlap;

                tbox::plog << "outer side overlap, overwrite_interior == false:\n";
                oside_overlap = oside_dst_geometry.calculateOverlap(side_dst_geometry,
                                                                    oside_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ false,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                for (int axis = 0; axis < NDIM; ++axis) oside_overlap->getDestinationBoxList(axis).print(tbox::plog);

                tbox::plog << "outer side overlap, overwrite_interior == true:\n";
                oside_overlap = oside_dst_geometry.calculateOverlap(side_dst_geometry,
                                                                    oside_src_geometry,
                                                                    src_mask,
                                                                    /*overwrite_interior*/ true,
                                                                    /*src_offset*/ IntVector<NDIM>(0),
                                                                    /*retry*/ true);
                for (int axis = 0; axis < NDIM; ++axis) oside_overlap->getDestinationBoxList(axis).print(tbox::plog);

#if (NDIM == 3)
            }
#endif
        }
    }
}
