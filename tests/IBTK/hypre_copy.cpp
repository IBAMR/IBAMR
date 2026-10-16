// ---------------------------------------------------------------------
//
// Copyright (c) 2026 - 2026 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Copy cell- and side-centered data to hypre vectors and back.

#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/solver_utilities.h>

#include <Box.h>
#include <CellData.h>
#include <CellIterator.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIterator.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <vector>

#include <ibtk/app_namespaces.h>

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    initialize_hypre();
    std::ofstream out("output");

    // The source data have two layers of ghost cells and the destination data
    // have one, so neither array has the layout of the box that is copied.
    const Box<NDIM> box(hier::Index<NDIM>(2), hier::Index<NDIM>(9));
    const int src_ghosts = 2, dst_ghosts = 1;
    const MPI_Comm communicator = IBTK_MPI::getCommunicator();
    std::array<HYPRE_Int, NDIM> lower = hypre_array(box.lower()), upper = hypre_array(box.upper());
    const auto value = [](const hier::Index<NDIM>& i, const int component)
    {
        double v = 1.0 + 1000.0 * component;
        for (int d = 0; d < NDIM; ++d)
        {
            v += std::pow(10.0, d) * i(d);
        }
        return v;
    };

    {
        const int depth = 2;
        HYPRE_StructGrid grid;
        HYPRE_StructGridCreate(communicator, NDIM, &grid);
        HYPRE_StructGridSetExtents(grid, lower.data(), upper.data());
        HYPRE_StructGridAssemble(grid);
        std::vector<HYPRE_StructVector> vectors(depth);
        for (auto& vector : vectors)
        {
            HYPRE_StructVectorCreate(communicator, grid, &vector);
            HYPRE_StructVectorInitialize(vector);
        }

        CellData<NDIM, double> src(box, depth, src_ghosts), dst(box, depth, dst_ghosts);
        for (int k = 0; k < depth; ++k)
        {
            for (CellIterator<NDIM> ic(src.getGhostBox()); ic; ic++)
            {
                src(ic(), k) = value(ic(), k);
            }
        }
        dst.fillAll(-1.0);
        copyToHypre(vectors, src, box);
        for (auto& vector : vectors)
        {
            HYPRE_StructVectorAssemble(vector);
        }
        copyFromHypre(dst, vectors, box);

        double max_difference = 0.0;
        bool difference_is_finite = true;
        for (int k = 0; k < depth; ++k)
        {
            for (CellIterator<NDIM> ic(box); ic; ic++)
            {
                const double difference = std::abs(dst(ic(), k) - src(ic(), k));
                difference_is_finite = difference_is_finite && std::isfinite(difference);
                max_difference = std::max(max_difference, difference);
            }
        }
        if (!difference_is_finite)
        {
            TBOX_ERROR("The copied cell-centered data are not finite.\n");
        }
        out << "cell-centered maximum difference = " << max_difference << "\n";

        for (auto& vector : vectors)
        {
            HYPRE_StructVectorDestroy(vector);
        }
        HYPRE_StructGridDestroy(grid);
    }

    {
        const int part = 0;
        HYPRE_SStructGrid grid;
        HYPRE_SStructGridCreate(communicator, NDIM, 1, &grid);
        HYPRE_SStructGridSetExtents(grid, part, lower.data(), upper.data());
#if (NDIM == 2)
        HYPRE_SStructVariable variables[NDIM] = { HYPRE_SSTRUCT_VARIABLE_XFACE, HYPRE_SSTRUCT_VARIABLE_YFACE };
#endif
#if (NDIM == 3)
        HYPRE_SStructVariable variables[NDIM] = { HYPRE_SSTRUCT_VARIABLE_XFACE,
                                                  HYPRE_SSTRUCT_VARIABLE_YFACE,
                                                  HYPRE_SSTRUCT_VARIABLE_ZFACE };
#endif
        HYPRE_SStructGridSetVariables(grid, part, NDIM, variables);
        HYPRE_SStructGridAssemble(grid);
        HYPRE_SStructVector vector;
        HYPRE_SStructVectorCreate(communicator, grid, &vector);
        HYPRE_SStructVectorInitialize(vector);

        SideData<NDIM, double> src(box, 1, src_ghosts), dst(box, 1, dst_ghosts);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (SideIterator<NDIM> is(src.getGhostBox(), axis); is; is++)
            {
                src(is()) = value(is(), axis);
            }
        }
        dst.fillAll(-1.0);
        copyToHypre(vector, src, box);
        HYPRE_SStructVectorAssemble(vector);
        HYPRE_SStructVectorGather(vector);
        copyFromHypre(dst, vector, box);

        double max_difference = 0.0;
        bool difference_is_finite = true;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (SideIterator<NDIM> is(box, axis); is; is++)
            {
                const double difference = std::abs(dst(is()) - src(is()));
                difference_is_finite = difference_is_finite && std::isfinite(difference);
                max_difference = std::max(max_difference, difference);
            }
        }
        if (!difference_is_finite)
        {
            TBOX_ERROR("The copied side-centered data are not finite.\n");
        }
        out << "side-centered maximum difference = " << max_difference << "\n";

        HYPRE_SStructVectorDestroy(vector);
        HYPRE_SStructGridDestroy(grid);
    }
}
