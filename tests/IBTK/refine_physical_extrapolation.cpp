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

// Fill the ghost cells of the finest level with linear refinement and extrapolation at the physical boundary, without
// coarse-fine interpolation, and check the ghost cells between the lower boundary and a patch within the ghost cell
// width of it against a linear function, which both reproduce. The temporary coarse patch of the refinement is one
// cell wide there, so the extrapolation must declare the stencil width that it uses.

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/ibtk_utilities.h>

#include <petscsys.h>

#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellIterator.h>

#include <string>
#include <vector>

#include "../coarse_fine_ghost_values.h"

#include <ibtk/app_namespaces.h>

using namespace coarse_fine_ghost_values;

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "ghost_cells.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        const int ghost_width = input_db->getIntegerWithDefault("GHOST_WIDTH", 1);
        const std::string extrap_type = input_db->getString("EXTRAP_TYPE");

        Pointer<hier::Variable<NDIM>> Q_var = new CellVariable<NDIM, double>("Q");
        const auto Q_fcn = [](const VectorNd& x)
        {
            double val = 1.0;
            for (int d = 0; d < NDIM; ++d)
            {
                val += x[d];
            }
            return val;
        };

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int Q_idx =
            var_db->registerVariableAndContext(Q_var, var_db->getContext("CTX"), IntVector<NDIM>(ghost_width));

        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = make_hierarchy(app_initializer, grid_geometry);
        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(Q_idx);
        }
        fill_quadratic(patch_hierarchy, { Q_idx }, Q_fcn);

        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_cell_comps(1);
        ghost_cell_comps[0] = InterpolationTransactionComponent(Q_idx,
                                                                "CONSERVATIVE_LINEAR_REFINE",
                                                                false,
                                                                "CONSERVATIVE_COARSEN",
                                                                extrap_type,
                                                                false,
                                                                nullptr /*bdry_conds*/);
        HierarchyGhostCellInterpolation hier_ghost_cell;
        hier_ghost_cell.initializeOperatorState(ghost_cell_comps, patch_hierarchy);
        hier_ghost_cell.fillData(0.0);

        // Check the ghost cells below each patch, away from its corners and from the other boundaries.
        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(finest_ln);
        const FineLevel fine = make_fine_level(patch_hierarchy, grid_geometry);
        int num_checked = 0;
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> Q_data = patch->getPatchData(Q_idx);
            Box<NDIM> below = patch->getBox();
            below.upper(1) = below.lower(1) - 1;
            below.lower(1) -= ghost_width;
            for (int d = 0; d < NDIM; ++d)
            {
                if (d != 1)
                {
                    below.lower(d) += 2;
                    below.upper(d) -= 2;
                }
            }
            for (CellIterator<NDIM> ci(below); ci; ci++)
            {
                const CellIndex<NDIM>& idx = ci();
                if (!in_boxes(fine.domain, -1, idx))
                {
                    continue;
                }
                const double exact = exact_value(patch, -1, idx, Q_fcn);
                ++num_checked;
                if (!IBTK::abs_equal_eps((*Q_data)(idx)-exact, 0.0))
                {
                    TBOX_ERROR("refine_physical_extrapolation:\n"
                               << "  incorrect ghost value on cell " << idx << ", patch " << patch->getBox()
                               << ": computed " << (*Q_data)(idx) << ", exact " << exact << "\n");
                }
            }
        }
        pout << "Checked " << num_checked << " values below patches.\n";

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->deallocatePatchData(Q_idx);
        }
    } // cleanup dynamically allocated objects prior to shutdown
} // main
