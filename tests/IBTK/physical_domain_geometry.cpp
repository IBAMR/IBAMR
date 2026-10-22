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
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PhysicalBoundaryUtilities.h>
#include <ibtk/ShiftedPatchGeometry.h>

#include <tbox/Pointer.h>
#include <tbox/Utilities.h>

#include <BergerRigoutsos.h>
#include <Box.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <StandardTagAndInitialize.h>

#include <array>
#include <fstream>

#include <ibtk/app_namespaces.h>

// Test of the shifted patch geometry and of the periodic physical domain and boundary-face test of the physical
// boundary utilities.

int
main(int argc, char* argv[])
{
    // Initialize IBAMR and libraries. Deinitialization is handled by this object as well.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input file, and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "physical_domain_geometry.log");

        // Create a single-level hierarchy of several patches.
        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
            "StandardTagAndInitialize", nullptr, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM>> load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);

        // Shift the geometry of each patch along each axis by half a cell, and check that the geometry is otherwise
        // the same, and that it is restored at the end of the scope.
        std::array<double, NDIM> total_shift;
        total_shift.fill(0.0);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                std::array<double, NDIM> shift;
                shift.fill(0.0);
                shift[axis] = -0.5 * pgeom->getDx()[axis];
                const ShiftedPatchGeometry shifted_pgeom(*patch, shift);
                {
                    const ShiftedPatchGeometry::Scope scope(*patch, shifted_pgeom);
                    Pointer<CartesianPatchGeometry<NDIM>> pgeom_shifted = patch->getPatchGeometry();
                    if (pgeom_shifted == pgeom || pgeom_shifted->getRatio() != pgeom->getRatio())
                    {
                        TBOX_ERROR("the shifted patch geometry is the original geometry or has another ratio\n");
                    }
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        if (pgeom_shifted->getDx()[d] != pgeom->getDx()[d] ||
                            pgeom_shifted->getXLower()[d] != pgeom->getXLower()[d] + shift[d] ||
                            pgeom_shifted->getXUpper()[d] != pgeom->getXUpper()[d] + shift[d])
                        {
                            TBOX_ERROR("the shifted patch geometry has the wrong grid spacing or corners\n");
                        }
                        for (int upperlower = 0; upperlower < 2; ++upperlower)
                        {
                            if (pgeom_shifted->getTouchesRegularBoundary(d, upperlower) !=
                                    pgeom->getTouchesRegularBoundary(d, upperlower) ||
                                pgeom_shifted->getTouchesPeriodicBoundary(d, upperlower) !=
                                    pgeom->getTouchesPeriodicBoundary(d, upperlower))
                            {
                                TBOX_ERROR("the shifted patch geometry has the wrong boundary flags\n");
                            }
                        }
                    }
                    total_shift[axis] +=
                        (pgeom_shifted->getXLower()[axis] - pgeom->getXLower()[axis]) / pgeom->getDx()[axis];
                }
                if (patch->getPatchGeometry() != pgeom)
                {
                    TBOX_ERROR("the patch geometry is not restored at the end of the scope\n");
                }
            }
        }

        // Count the cells and the faces on the boundaries of the periodic physical domain.
        IBTK_MPI::sumReduction(total_shift.data(), NDIM);
        const BoxArray<NDIM> domain =
            PhysicalBoundaryUtilities::getPeriodicPhysicalDomain(grid_geometry, IntVector<NDIM>(1));
        const Box<NDIM> domain_box = grid_geometry->getPhysicalDomain()[0];
        const Box<NDIM> search_box = Box<NDIM>::grow(domain_box, domain_box.numberCells());
        int num_domain_cells = 0;
        for (Box<NDIM>::Iterator b(search_box); b; b++)
        {
            if (domain.contains(b()))
            {
                ++num_domain_cells;
            }
        }
        std::array<std::array<int, 2>, NDIM> num_boundary_faces;
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            for (int upperlower = 0; upperlower < 2; ++upperlower)
            {
                const bool is_lower = upperlower == 0;
                Box<NDIM> face_plane = search_box;
                face_plane.lower(axis) = face_plane.upper(axis) =
                    is_lower ? domain_box.lower(axis) : domain_box.upper(axis) + 1;
                num_boundary_faces[axis][upperlower] = 0;
                for (Box<NDIM>::Iterator b(face_plane); b; b++)
                {
                    if (PhysicalBoundaryUtilities::isBoundaryFace(b(), axis, is_lower, domain))
                    {
                        ++num_boundary_faces[axis][upperlower];
                    }
                }
            }
        }

        if (IBTK_MPI::getRank() == 0)
        {
            std::ofstream out("output");
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                out << "sum over patches of the shift of the lower corner along axis " << axis
                    << ", in cells = " << total_shift[axis] << '\n';
            }
            out << "number of cells in the periodic physical domain = " << num_domain_cells << '\n';
            for (unsigned int axis = 0; axis < NDIM; ++axis)
            {
                out << "faces on the lower boundary normal to axis " << axis << " = " << num_boundary_faces[axis][0]
                    << '\n';
                out << "faces on the upper boundary normal to axis " << axis << " = " << num_boundary_faces[axis][1]
                    << '\n';
            }
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
