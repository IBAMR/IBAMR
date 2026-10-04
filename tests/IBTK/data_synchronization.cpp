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
#include <ibtk/EdgeDataSynchronization.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/NodeDataSynchronization.h>

#include <tbox/Database.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>

#include <ArrayData.h>
#include <Box.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <EdgeData.h>
#include <EdgeGeometry.h>
#include <EdgeVariable.h>
#include <Index.h>
#include <IntVector.h>
#include <NodeData.h>
#include <NodeGeometry.h>
#include <NodeVariable.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <VariableDatabase.h>

#include <array>
#include <cstddef>
#include <functional>
#include <set>
#include <string>
#include <utility>

#include <ibtk/app_namespaces.h>

// Node and edge data have values on patch borders, and such a value is stored
// by every patch that touches it. Synchronizing the data gives all of these
// copies the value of one of them, the owner: the copy of the patch that
// contains the highest of the cells touching the value, comparing the last
// coordinate first. A periodic image of a patch counts as another patch.
//
// This test makes a hierarchy from the patches listed in the input file, gives
// every copy a different value, synchronizes the data on each level, and
// counts the copies that hold the original value of the owner.

namespace
{
// The refinement ratio between successive levels.
const int RATIO = 2;
// Larger than any index in the hierarchy.
const double MARK_BASE = 64.0;

// The index box of array n for data on the given box of cells.
using ArrayBoxFcn = std::function<Box<NDIM>(const Box<NDIM>&, int)>;

ArrayData<NDIM, double>&
get_array(NodeData<NDIM, double>& data, const int /*n*/)
{
    return data.getArrayData();
}

ArrayData<NDIM, double>&
get_array(EdgeData<NDIM, double>& data, const int n)
{
    return data.getArrayData(n);
}

// A number that identifies the copy that patch k holds of the value at index i
// of array n.
double
mark(const int k, const int n, const hier::Index<NDIM>& i)
{
    double value = k * NDIM + n;
    for (int d = 0; d < NDIM; ++d)
    {
        value = MARK_BASE * value + i(d);
    }
    return value;
}

// Find the patches of a level that touch the value at index i of array n,
// where a periodic image of a patch counts as another patch. Return their
// number and the mark of the copy that owns the value.
std::pair<int, double>
find_owner(Pointer<PatchLevel<NDIM>> level, const ArrayBoxFcn& array_box, const int n, const hier::Index<NDIM>& i)
{
    const BoxArray<NDIM>& boxes = level->getBoxes();
    const Box<NDIM>& domain = level->getPhysicalDomain()[0];
    const IntVector<NDIM> periodic_shift = level->getGridGeometry()->getPeriodicShift(level->getRatio());
    std::set<std::array<int, NDIM + 1>> touching_patches;
    double owner_mark = 0.0;
    // The iterator visits the cells in increasing order, comparing the last
    // coordinate first, so the last cell that is in a patch is the highest.
    for (Box<NDIM>::Iterator c(Box<NDIM>(i - 1, i)); c; c++)
    {
        if (!array_box(Box<NDIM>(c(), c()), n).contains(i))
        {
            continue;
        }
        // A cell outside a periodic domain is the image of the cell that the
        // periodic shift takes to it.
        IntVector<NDIM> shift(0);
        for (int d = 0; d < NDIM; ++d)
        {
            shift(d) = c()(d) < domain.lower(d) ? -periodic_shift(d) : c()(d) > domain.upper(d) ? periodic_shift(d) : 0;
        }
        for (int k = 0; k < boxes.size(); ++k)
        {
            if (boxes[k].contains(c() - shift))
            {
                std::array<int, NDIM + 1> patch_image;
                patch_image[0] = k;
                for (int d = 0; d < NDIM; ++d)
                {
                    patch_image[d + 1] = shift(d);
                }
                touching_patches.insert(patch_image);
                owner_mark = mark(k, n, i - shift);
            }
        }
    }
    return { static_cast<int>(touching_patches.size()), owner_mark };
}

// Synchronize data of one centering on each level of the hierarchy and print
// the number of copies of the values touched by a given number of patches,
// followed by the number of those copies that hold the value of the owner.
template <class VariableType, class DataType, class SynchronizationType>
void
test_synchronization(const std::string& name,
                     const int n_arrays,
                     const ArrayBoxFcn& array_box,
                     Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableType> var = new VariableType(name);
    const int idx = var_db->registerVariableAndContext(var, var_db->getContext("context"), IntVector<NDIM>(0));
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->allocatePatchData(idx);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<DataType> data = patch->getPatchData(idx);
            for (int n = 0; n < n_arrays; ++n)
            {
                for (Box<NDIM>::Iterator i(array_box(patch->getBox(), n)); i; i++)
                {
                    get_array(*data, n)(i(), 0) = mark(p(), n, i());
                }
            }
        }
    }

    // The levels are synchronized separately because there is no coarsening.
    SynchronizationType synch_op;
    synch_op.initializeOperatorState(typename SynchronizationType::SynchronizationTransactionComponent(idx, "NONE"),
                                     hierarchy);
    synch_op.synchronizeData(0.0);

    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        // Entry m - 1 is for the values touched by m patches.
        std::array<int, 1 << NDIM> n_copies{}, n_owner_copies{};
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<DataType> data = patch->getPatchData(idx);
            for (int n = 0; n < n_arrays; ++n)
            {
                for (Box<NDIM>::Iterator i(array_box(patch->getBox(), n)); i; i++)
                {
                    const std::pair<int, double> owner = find_owner(level, array_box, n, i());
                    ++n_copies[owner.first - 1];
                    if (get_array(*data, n)(i(), 0) == owner.second)
                    {
                        ++n_owner_copies[owner.first - 1];
                    }
                }
            }
        }
        level->deallocatePatchData(idx);
        IBTK_MPI::sumReduction(n_copies.data(), static_cast<int>(n_copies.size()));
        IBTK_MPI::sumReduction(n_owner_copies.data(), static_cast<int>(n_owner_copies.size()));

        plog << name << " data, level " << ln << '\n';
        for (std::size_t m = 1; m <= n_copies.size(); ++m)
        {
            if (n_copies[m - 1] > 0)
            {
                plog << "  values touched by " << m << (m == 1 ? " patch: " : " patches: ") << n_copies[m - 1]
                     << " copies, " << n_owner_copies[m - 1] << " with the value of the owner\n";
            }
        }
    }
}
} // namespace

int
main(int argc, char** argv)
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

    // Patch k of a level is on processor k mod the number of processors, so
    // that the patches are the same on any number of processors.
    Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<Database> patches_db = app_initializer->getComponentDatabase("Patches");
    IntVector<NDIM> ratio(1);
    for (int ln = 0; patches_db->keyExists("level_" + std::to_string(ln)); ++ln)
    {
        const BoxArray<NDIM> boxes = patches_db->getDatabaseBoxArray("level_" + std::to_string(ln));
        ProcessorMapping mapping(boxes.size());
        for (int k = 0; k < boxes.size(); ++k)
        {
            mapping.setProcessorAssignment(k, k % IBTK_MPI::getNodes());
        }
        hierarchy->makeNewPatchLevel(ln, ratio, boxes, mapping);
        ratio *= RATIO;
    }

    test_synchronization<NodeVariable<NDIM, double>, NodeData<NDIM, double>, NodeDataSynchronization>(
        "Node", 1, [](const Box<NDIM>& box, const int /*n*/) { return NodeGeometry<NDIM>::toNodeBox(box); }, hierarchy);
    test_synchronization<EdgeVariable<NDIM, double>, EdgeData<NDIM, double>, EdgeDataSynchronization>(
        "Edge", NDIM, EdgeGeometry<NDIM>::toEdgeBox, hierarchy);
}
