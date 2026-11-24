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

#include <BoundaryBox.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <CoarseFineBoundary.h>
#include <PatchHierarchy.h>
#include <ProcessorMapping.h>
#include <VariableDatabase.h>

#include <string>

#include <ibtk/app_namespaces.h>

namespace
{
// Print each boundary box with its location index.
void
print_boxes(const std::string& label, int codim, const Array<BoundaryBox<NDIM>>& boxes)
{
    for (int k = 0; k < boxes.size(); ++k)
    {
        plog << label << " codim " << codim << " location " << boxes[k].getLocationIndex() << " box "
             << boxes[k].getBox() << '\n';
    }
}

// Return a grid of number(d) tiles of the given size, separated by gap.
BoxArray<NDIM>
make_tiles(Pointer<Database> db)
{
    int lower[NDIM] = {}, size[NDIM] = {}, gap[NDIM] = {}, number[NDIM] = {};
    db->getIntegerArray("lower", lower, NDIM);
    db->getIntegerArray("size", size, NDIM);
    db->getIntegerArray("gap", gap, NDIM);
    db->getIntegerArray("number", number, NDIM);
    int num_tiles = 1;
    for (int d = 0; d < NDIM; ++d)
    {
        num_tiles *= number[d];
    }
    BoxArray<NDIM> tiles(num_tiles);
    for (int t = 0; t < num_tiles; ++t)
    {
        int rest = t;
        for (int d = 0; d < NDIM; ++d)
        {
            const int i = rest % number[d];
            rest /= number[d];
            tiles[t].lower(d) = lower[d] + i * (size[d] + gap[d]);
            tiles[t].upper(d) = tiles[t].lower(d) + size[d] - 1;
        }
    }
    return tiles;
}
} // namespace

int
main(int argc, char** argv)
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "boundary_boxes.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        const IntVector<NDIM> ghost_width(input_db->getInteger("GHOST_WIDTH"));
        const IntVector<NDIM> ratio(input_db->getInteger("REF_RATIO"));

        // The physical boundary boxes extend through the largest ghost width of the patch descriptor.
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<CellVariable<NDIM, double>> u_var = new CellVariable<NDIM, double>("u");
        var_db->registerVariableAndContext(u_var, var_db->getContext("context"), ghost_width);

        Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        for (int ln = 0; ln < 2; ++ln)
        {
            const BoxArray<NDIM> boxes = make_tiles(input_db->getDatabase("Level" + std::to_string(ln)));
            // Put every patch on rank 0 so that the output does not depend on the number of processes.
            ProcessorMapping mapping(boxes.size());
            for (int k = 0; k < boxes.size(); ++k)
            {
                mapping.setProcessorAssignment(k, 0);
            }
            hierarchy->makeNewPatchLevel(ln, ln == 0 ? IntVector<NDIM>(1) : ratio, boxes, mapping);
        }

        for (int ln = 0; ln < 2; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            plog << "level " << ln << " patches: " << level->getNumberOfPatches() << '\n';
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                const std::string label = "level " + std::to_string(ln) + " patch " + std::to_string(p()) + " physical";
                for (int codim = 1; codim <= NDIM; ++codim)
                {
                    print_boxes(
                        label, codim, level->getPatch(p())->getPatchGeometry()->getCodimensionBoundaries(codim));
                }
            }
        }

        CoarseFineBoundary<NDIM> coarse_fine_boundary(*hierarchy, 1, ghost_width);
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(1);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const std::string label = "level 1 patch " + std::to_string(p()) + " coarse-fine";
            for (int codim = 1; codim <= NDIM; ++codim)
            {
                print_boxes(label, codim, coarse_fine_boundary.getBoundaries(p(), codim));
            }
        }
    }
}
