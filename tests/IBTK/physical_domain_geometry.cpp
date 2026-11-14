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

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/PhysicalBoundaryUtilities.h>

#include <tbox/Pointer.h>

#include <Box.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <IntVector.h>

#include <ibtk/app_namespaces.h>

// Test of IBTK::PhysicalBoundaryUtilities::getPeriodicPhysicalDomain() and
// IBTK::PhysicalBoundaryUtilities::isBoundaryFace().

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "physical_domain_geometry.log");

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

        // Count the cells and the faces on the boundaries of the periodic physical domain, at ratios 1 and 2.
        for (int ratio = 1; ratio <= 2; ++ratio)
        {
            const BoxArray<NDIM> domain =
                PhysicalBoundaryUtilities::getPeriodicPhysicalDomain(grid_geometry, IntVector<NDIM>(ratio));
            const Box<NDIM> domain_box = Box<NDIM>::refine(grid_geometry->getPhysicalDomain()[0], ratio);
            const Box<NDIM> search_box = Box<NDIM>::grow(domain_box, domain_box.numberCells());
            int num_domain_cells = 0;
            for (Box<NDIM>::Iterator b(search_box); b; b++)
            {
                if (domain.contains(b()))
                {
                    ++num_domain_cells;
                }
            }
            plog << "ratio " << ratio << ": cells in the periodic physical domain = " << num_domain_cells << '\n';
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                for (int upperlower = 0; upperlower < 2; ++upperlower)
                {
                    const bool is_lower = upperlower == 0;
                    Box<NDIM> face_plane = search_box;
                    face_plane.lower(axis) = face_plane.upper(axis) =
                        is_lower ? domain_box.lower(axis) : domain_box.upper(axis) + 1;
                    int num_boundary_faces = 0;
                    for (Box<NDIM>::Iterator b(face_plane); b; b++)
                    {
                        if (PhysicalBoundaryUtilities::isBoundaryFace(b(), axis, is_lower, domain))
                        {
                            ++num_boundary_faces;
                        }
                    }
                    plog << "ratio " << ratio << ": faces on the " << (is_lower ? "lower" : "upper")
                         << " boundary normal to axis " << axis << " = " << num_boundary_faces << '\n';
                }
            }
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
