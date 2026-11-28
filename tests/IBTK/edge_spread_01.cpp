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

// Check that the LEInteractor::spread() overload for edge-centered data that takes the positions only gives the same
// grid values as the overload that takes index data, for points in the interior of one patch.

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/LNode.h>
#include <ibtk/LNodeSet.h>
#include <ibtk/LNodeSetData.h>

#include <tbox/Pointer.h>

#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CellIndex.h>
#include <EdgeData.h>
#include <EdgeGeometry.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
// The largest difference between the two edge-centered fields and the largest magnitude of the first, over the ghost
// box of every component.
std::pair<double, double>
max_difference(const EdgeData<NDIM, double>& a, const EdgeData<NDIM, double>& b)
{
    double difference = 0.0, largest = 0.0;
    for (unsigned int axis = 0; axis < NDIM; ++axis)
    {
        const Box<NDIM> edge_box = EdgeGeometry<NDIM>::toEdgeBox(a.getGhostBox(), axis);
        for (Box<NDIM>::Iterator it(edge_box); it; it++)
        {
            const double value = a.getArrayData(axis)(it(), 0);
            difference = std::max(difference, std::abs(value - b.getArrayData(axis)(it(), 0)));
            largest = std::max(largest, std::abs(value));
        }
    }
    return { difference, largest };
} // max_difference
} // namespace

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
    {
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        Pointer<CartesianGridGeometry<NDIM>> grid_geom = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        const int N = input_db->getInteger("N");
        const int num_points = input_db->getInteger("num_points");
        const std::vector<std::string> kernels = { "PIECEWISE_LINEAR", "IB_4", "IB_6", "BSPLINE_3" };
        int ghost_width = 0;
        for (const std::string& kernel : kernels)
        {
            ghost_width = std::max(ghost_width, LEInteractor::getMinimumGhostWidth(kernel));
        }

        // One patch covering the domain.
        BoxArray<NDIM> boxes(1);
        boxes[0] = Box<NDIM>(hier::Index<NDIM>(0), hier::Index<NDIM>(N - 1));
        ProcessorMapping mapping(1);
        mapping.setProcessorAssignment(0, 0);
        Pointer<PatchLevel<NDIM>> level = new PatchLevel<NDIM>(
            boxes, mapping, IntVector<NDIM>(1), grid_geom, VariableDatabase<NDIM>::getDatabase()->getPatchDescriptor());
        Pointer<Patch<NDIM>> patch = level->getPatch(0);
        const Box<NDIM>& patch_box = patch->getBox();

        // Points at golden-ratio positions at least one cell from every face of the patch, so that both overloads
        // select every point, with three values each.
        const double golden = 0.5 * (std::sqrt(5.0) - 1.0);
        std::vector<double> X(NDIM * num_points), F(NDIM * num_points);
        for (int k = 0; k < num_points; ++k)
        {
            for (int d = 0; d < NDIM; ++d)
            {
                const double h = (grid_geom->getXUpper()[d] - grid_geom->getXLower()[d]) / N;
                const double fraction = std::fmod((k + 1) * golden * (d + 1) + 0.37 * d, 1.0);
                X[NDIM * k + d] = grid_geom->getXLower()[d] + h * (1.0 + fraction * (N - 2));
                F[NDIM * k + d] = 1.0 + 0.1 * d + std::sin(1.3 * k + 0.7 * d);
            }
        }

        // Index data with each point in the cell of its position. The local PETSc index of a point is its number.
        Pointer<LNodeSetData> idx_data = new LNodeSetData(patch_box, IntVector<NDIM>(ghost_width));
        for (int k = 0; k < num_points; ++k)
        {
            CellIndex<NDIM> i(hier::Index<NDIM>(0));
            for (int d = 0; d < NDIM; ++d)
            {
                const double h = (grid_geom->getXUpper()[d] - grid_geom->getXLower()[d]) / N;
                i(d) = static_cast<int>(std::floor((X[NDIM * k + d] - grid_geom->getXLower()[d]) / h));
            }
            if (!idx_data->isElement(i))
            {
                idx_data->appendItemPointer(i, new LNodeSet());
            }
            idx_data->getItem(i)->push_back(new LNode(k, k, k));
        }
        idx_data->cacheLocalIndices(patch, grid_geom, IntVector<NDIM>(1));

        for (const std::string& kernel : kernels)
        {
            Pointer<EdgeData<NDIM, double>> q_index =
                new EdgeData<NDIM, double>(patch_box, 1, IntVector<NDIM>(ghost_width));
            Pointer<EdgeData<NDIM, double>> q_positions =
                new EdgeData<NDIM, double>(patch_box, 1, IntVector<NDIM>(ghost_width));
            q_index->fillAll(0.0);
            q_positions->fillAll(0.0);
            LEInteractor::spread(q_index, F.data(), NDIM, X.data(), NDIM, idx_data, patch, patch_box, kernel);
            LEInteractor::spread(q_positions, F, NDIM, X, NDIM, patch, patch_box, kernel);
            const std::pair<double, double> result = max_difference(*q_index, *q_positions);
            if (!(result.first <= 1.0e-12 * result.second))
            {
                TBOX_ERROR(kernel << ": the two spread overloads differ by " << result.first
                                  << ", relative to the largest value " << result.second << "\n");
            }
            plog << kernel << ": largest value = " << result.second << "\n";
        }
    }
    return 0;
} // main
