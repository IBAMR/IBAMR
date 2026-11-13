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

// Fill the ghost values of the finest level with the refine operator alone, without coarse-fine interpolation. Check
// that all patches that contain a face of the fine level, as an interior or a ghost value, agree on its value, and that
// the fill does not change the interior values.

#include <ibtk/AppInitializer.h>
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/SideNoCornersFillPattern.h>

#include <tbox/Logger.h>

#include <CartesianGridGeometry.h>

#include <sstream>
#include <string>
#include <vector>

#include "../coarse_fine_ghost_values.h"
#include "../tests.h"

#include <ibtk/app_namespaces.h>

using namespace coarse_fine_ghost_values;

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    Pointer<Logger::Appender> abort_append(new TestAppender());
    Logger::getInstance()->setAbortAppender(abort_append);

    { // cleanup dynamically allocated objects prior to shutdown
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "ghost_cells.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        const int ghost_width = input_db->getIntegerWithDefault("GHOST_WIDTH", 1);
        const std::string refine_op_name = input_db->getStringWithDefault("REFINE_OP", "NONE");
        const bool is_cell = input_db->getString("CENTERING").compare("CELL") == 0;

        // Create cell-centered or side-centered quantities, and initialize them with a quadratic function. The second
        // quantity is filled with a no-corners fill pattern.
        Pointer<hier::Variable<NDIM>> Q_var, Q_nc_var;
        if (is_cell)
        {
            Q_var = new CellVariable<NDIM, double>("Q");
            Q_nc_var = new CellVariable<NDIM, double>("Q_no_corners");
        }
        else
        {
            Q_var = new SideVariable<NDIM, double>("Q");
            Q_nc_var = new SideVariable<NDIM, double>("Q_no_corners");
        }
        const auto Q_fcn = make_quadratic(grid_geometry);

        auto var_db = VariableDatabase<NDIM>::getDatabase();
        const int Q_idx =
            var_db->registerVariableAndContext(Q_var, var_db->getContext("CTX"), IntVector<NDIM>(ghost_width));
        const int Q_nc_idx =
            var_db->registerVariableAndContext(Q_nc_var, var_db->getContext("CTX"), IntVector<NDIM>(ghost_width));

        Pointer<PatchHierarchy<NDIM>> patch_hierarchy = make_hierarchy(app_initializer, grid_geometry);
        const int coarsest_ln = 0, finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(Q_idx);
            level->allocatePatchData(Q_nc_idx);
        }
        fill_quadratic(patch_hierarchy, { Q_idx, Q_nc_idx }, Q_fcn);

        // Fill the ghost cells with the refine operator only, with the default fill pattern and with a no-corners fill
        // pattern.
        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> ghost_cell_comps(1);
        ghost_cell_comps[0] = InterpolationTransactionComponent(
            Q_idx, refine_op_name, false, "NONE", "NONE", false, nullptr /*bdry_conds*/);
        HierarchyGhostCellInterpolation hier_ghost_cell;
        hier_ghost_cell.initializeOperatorState(ghost_cell_comps, patch_hierarchy);
        hier_ghost_cell.fillData(0.0);

        Pointer<VariableFillPattern<NDIM>> no_corners_pattern;
        if (is_cell)
        {
            no_corners_pattern = new CellNoCornersFillPattern(ghost_width, /*overwrite_interior*/ false);
        }
        else
        {
            no_corners_pattern = new SideNoCornersFillPattern(ghost_width, /*overwrite_interior*/ false);
        }
        std::vector<InterpolationTransactionComponent> no_corners_comps(1);
        no_corners_comps[0] = InterpolationTransactionComponent(Q_nc_idx,
                                                                Q_nc_idx,
                                                                refine_op_name,
                                                                false,
                                                                "NONE",
                                                                "NONE",
                                                                false,
                                                                nullptr /*bdry_conds*/,
                                                                no_corners_pattern);
        HierarchyGhostCellInterpolation hier_ghost_cell_no_corners;
        hier_ghost_cell_no_corners.initializeOperatorState(no_corners_comps, patch_hierarchy);
        hier_ghost_cell_no_corners.fillData(0.0);

        Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(finest_ln);
        const FineLevel fine = make_fine_level(patch_hierarchy, grid_geometry);
        const std::vector<int> axes = data_axes(is_cell);

        std::ostringstream failures;
        check_fine_level_agreement(level, fine, grid_geometry, Q_idx, axes, is_cell, failures);
        check_interior_unchanged(level, { Q_idx, Q_nc_idx }, axes, Q_fcn, failures);
        if (!failures.str().empty())
        {
            TBOX_ERROR("interface_face_ghost_copy:\n" << failures.str());
        }

        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->deallocatePatchData(Q_idx);
            level->deallocatePatchData(Q_nc_idx);
        }
    } // cleanup dynamically allocated objects prior to shutdown
}
