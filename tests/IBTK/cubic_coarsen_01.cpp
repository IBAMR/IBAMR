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

#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/muParserCartGridFunction.h>

#include <tbox/Logger.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>

#include <cmath>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

#include <ibtk/app_namespaces.h>

// Coarsen a function from the finest level of a two-level hierarchy through
// HierarchyGhostCellInterpolation with cubic coarsening and print the largest
// error of the covered coarse values. With garbage_ghosts = TRUE, the ghost
// cells of the fine level hold NaN before the coarsening. A second fill follows,
// and the program prints the largest change that it makes to the covered coarse
// values.

namespace
{
// Set the ghost values of the fine level (the values outside the patch interiors) to the given value.
void
set_ghosts(Pointer<PatchLevel<NDIM>> level, const int cc_idx, const int sc_idx, const double value)
{
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        Pointer<CellData<NDIM, double>> cc_data = patch->getPatchData(cc_idx);
        for (CellIterator<NDIM> ic(cc_data->getGhostBox()); ic; ic++)
        {
            if (!patch_box.contains(ic()))
            {
                (*cc_data)(ic()) = value;
            }
        }
        Pointer<SideData<NDIM, double>> sc_data = patch->getPatchData(sc_idx);
        for (unsigned int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> side_box = SideGeometry<NDIM>::toSideBox(patch_box, axis);
            for (SideIterator<NDIM> is(sc_data->getGhostBox(), axis); is; is++)
            {
                if (!side_box.contains(is()))
                {
                    (*sc_data)(is()) = value;
                }
            }
        }
    }
}

// Return the largest difference between the coarse values and the reference values in the cells covered by the fine
// level.
double
max_covered_error(Pointer<PatchHierarchy<NDIM>> hierarchy, const int idx, const int exact_idx, const bool side)
{
    Pointer<PatchLevel<NDIM>> coarse_level = hierarchy->getPatchLevel(0);
    Pointer<PatchLevel<NDIM>> fine_level = hierarchy->getPatchLevel(1);
    BoxArray<NDIM> covered_boxes = fine_level->getBoxes();
    covered_boxes.coarsen(fine_level->getRatioToCoarserLevel());
    double max_error = 0.0;
    for (PatchLevel<NDIM>::Iterator p(coarse_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = coarse_level->getPatch(p());
        for (int k = 0; k < covered_boxes.size(); ++k)
        {
            const Box<NDIM> box = covered_boxes[k] * patch->getBox();
            if (box.empty())
            {
                continue;
            }
            if (side)
            {
                Pointer<SideData<NDIM, double>> data = patch->getPatchData(idx);
                Pointer<SideData<NDIM, double>> exact = patch->getPatchData(exact_idx);
                for (unsigned int axis = 0; axis < NDIM; ++axis)
                {
                    for (SideIterator<NDIM> is(box, axis); is; is++)
                    {
                        const double error = std::abs((*data)(is()) - (*exact)(is()));
                        if (!(error <= max_error))
                        {
                            max_error = error;
                        }
                    }
                }
            }
            else
            {
                Pointer<CellData<NDIM, double>> data = patch->getPatchData(idx);
                Pointer<CellData<NDIM, double>> exact = patch->getPatchData(exact_idx);
                for (CellIterator<NDIM> ic(box); ic; ic++)
                {
                    const double error = std::abs((*data)(ic()) - (*exact)(ic()));
                    if (!(error <= max_error))
                    {
                        max_error = error;
                    }
                }
            }
        }
    }
    return IBTK_MPI::maxReduction(max_error);
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    // Suppress the weighted-averaging warning, which prints the path of the source file.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
    TimerManager::createManager(nullptr);
    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cubic_coarsen.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
        const bool garbage_ghosts = input_db->getBoolWithDefault("garbage_ghosts", false);
        const int ghost_width = input_db->getIntegerWithDefault("ghost_width", 1);
        const std::string refine_op_name = input_db->getStringWithDefault("refine_op_name", "NONE");
        const tbox::Array<std::string> function_names = input_db->getStringArray("function_names");

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

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");
        Pointer<CellVariable<NDIM, double>> u_cc_var = new CellVariable<NDIM, double>("u_cc");
        Pointer<CellVariable<NDIM, double>> exact_cc_var = new CellVariable<NDIM, double>("exact_cc");
        Pointer<SideVariable<NDIM, double>> u_sc_var = new SideVariable<NDIM, double>("u_sc");
        Pointer<SideVariable<NDIM, double>> exact_sc_var = new SideVariable<NDIM, double>("exact_sc");
        const int u_cc_idx = var_db->registerVariableAndContext(u_cc_var, ctx, IntVector<NDIM>(ghost_width));
        const int exact_cc_idx = var_db->registerVariableAndContext(exact_cc_var, ctx, IntVector<NDIM>(0));
        const int u_sc_idx = var_db->registerVariableAndContext(u_sc_var, ctx, IntVector<NDIM>(ghost_width));
        const int exact_sc_idx = var_db->registerVariableAndContext(exact_sc_var, ctx, IntVector<NDIM>(0));

        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = std::numeric_limits<int>::max();
        int level_number = 0;
        while ((gridding_algorithm->levelCanBeRefined(level_number)))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }
        if (patch_hierarchy->getFinestLevelNumber() != 1)
        {
            TBOX_ERROR("cubic_coarsen_01:\n"
                       << "  the hierarchy must have exactly two levels\n");
        }
        for (int ln = 0; ln <= 1; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(u_cc_idx, 0.0);
            level->allocatePatchData(exact_cc_idx, 0.0);
            level->allocatePatchData(u_sc_idx, 0.0);
            level->allocatePatchData(exact_sc_idx, 0.0);
        }

        using ITC = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<ITC> transaction_comps = {
            ITC(u_cc_idx, refine_op_name, false, "CUBIC_COARSEN", "NONE", false, nullptr),
            ITC(u_sc_idx, refine_op_name, false, "CUBIC_COARSEN", "NONE", false, nullptr)
        };
        HierarchyGhostCellInterpolation ghost_fill;
        ghost_fill.initializeOperatorState(transaction_comps, patch_hierarchy);

        plog << std::scientific << std::setprecision(3);
        for (int k = 0; k < function_names.size(); ++k)
        {
            const std::string& name = function_names[k];
            muParserCartGridFunction cc_fcn(
                name + "_cc", app_initializer->getComponentDatabase(name + "_cc"), grid_geometry);
            muParserCartGridFunction sc_fcn(
                name + "_sc", app_initializer->getComponentDatabase(name + "_sc"), grid_geometry);
            cc_fcn.setDataOnPatchHierarchy(u_cc_idx, u_cc_var, patch_hierarchy, 0.0);
            cc_fcn.setDataOnPatchHierarchy(exact_cc_idx, exact_cc_var, patch_hierarchy, 0.0);
            sc_fcn.setDataOnPatchHierarchy(u_sc_idx, u_sc_var, patch_hierarchy, 0.0);
            sc_fcn.setDataOnPatchHierarchy(exact_sc_idx, exact_sc_var, patch_hierarchy, 0.0);

            // Reset the coarse values; the coarsening must supply the covered ones.
            Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(0);
            for (PatchLevel<NDIM>::Iterator p(coarse_level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = coarse_level->getPatch(p());
                Pointer<CellData<NDIM, double>> cc_data = patch->getPatchData(u_cc_idx);
                cc_data->fillAll(0.0);
                Pointer<SideData<NDIM, double>> sc_data = patch->getPatchData(u_sc_idx);
                sc_data->fillAll(0.0);
            }
            set_ghosts(patch_hierarchy->getPatchLevel(1),
                       u_cc_idx,
                       u_sc_idx,
                       garbage_ghosts ? std::numeric_limits<double>::quiet_NaN() : 0.0);

            ghost_fill.fillData(0.0);

            const double cc_error = max_covered_error(patch_hierarchy, u_cc_idx, exact_cc_idx, false);
            const double sc_error = max_covered_error(patch_hierarchy, u_sc_idx, exact_sc_idx, true);
            if (!std::isfinite(cc_error) || !std::isfinite(sc_error))
            {
                TBOX_ERROR("cubic_coarsen_01:\n"
                           << "  coarse values of function " << name << " are not finite\n");
            }
            // Print roundoff-level errors as zero.
            plog << name << " cell-centered max error = " << (cc_error < 1.0e-11 ? 0.0 : cc_error) << "\n";
            plog << name << " side-centered max error = " << (sc_error < 1.0e-11 ? 0.0 : sc_error) << "\n";

            // Take the filled coarse values as the reference and fill again.
            for (PatchLevel<NDIM>::Iterator p(coarse_level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = coarse_level->getPatch(p());
                Pointer<CellData<NDIM, double>> cc_data = patch->getPatchData(u_cc_idx);
                Pointer<CellData<NDIM, double>> cc_exact = patch->getPatchData(exact_cc_idx);
                cc_exact->copy(*cc_data);
                Pointer<SideData<NDIM, double>> sc_data = patch->getPatchData(u_sc_idx);
                Pointer<SideData<NDIM, double>> sc_exact = patch->getPatchData(exact_sc_idx);
                sc_exact->copy(*sc_data);
            }
            ghost_fill.fillData(0.0);
            const double cc_change = max_covered_error(patch_hierarchy, u_cc_idx, exact_cc_idx, false);
            const double sc_change = max_covered_error(patch_hierarchy, u_sc_idx, exact_sc_idx, true);
            if (!std::isfinite(cc_change) || !std::isfinite(sc_change))
            {
                TBOX_ERROR("cubic_coarsen_01:\n"
                           << "  coarse values of function " << name << " are not finite after the second fill\n");
            }
            plog << name << " cell-centered max change by second fill = " << (cc_change < 1.0e-11 ? 0.0 : cc_change)
                 << "\n";
            plog << name << " side-centered max change by second fill = " << (sc_change < 1.0e-11 ? 0.0 : sc_change)
                 << "\n";
        }
    }
    return 0;
}
