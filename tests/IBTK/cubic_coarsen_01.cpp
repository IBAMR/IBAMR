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
#include <sstream>
#include <string>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Coarsen functions from the finest level of a two-level hierarchy through HierarchyGhostCellInterpolation with cubic
// coarsening, with NaN in the ghost cells of the fine level. Print the largest error of the covered coarse values, and
// after a second fill the largest change. Optional inputs: corner_cells (a box of coarse cells with its own error),
// second_function_names (a second context, coarsened by the same fills), fine_patch_boxes (the required fine patches),
// and nan_cell (a fine cell set to NaN, with the number of covered coarse cells that are NaN printed in place of the
// errors).

namespace
{
// The patch data indices of the data that the fills coarsen and of their exact values.
struct Fields
{
    int u_cc_idx;
    int exact_cc_idx;
    int u_sc_idx;
    int exact_sc_idx;
};

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

// Set the value of one cell of the fine level.
void
set_cell(Pointer<PatchLevel<NDIM>> level, const int cc_idx, const CellIndex<NDIM>& i, const double value)
{
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        if (patch->getBox().contains(i))
        {
            Pointer<CellData<NDIM, double>> cc_data = patch->getPatchData(cc_idx);
            (*cc_data)(i) = value;
        }
    }
}

// Return the number of cells covered by the fine level whose coarse value is NaN.
int
count_covered_nans(Pointer<PatchHierarchy<NDIM>> hierarchy, const int cc_idx)
{
    Pointer<PatchLevel<NDIM>> coarse_level = hierarchy->getPatchLevel(0);
    Pointer<PatchLevel<NDIM>> fine_level = hierarchy->getPatchLevel(1);
    BoxArray<NDIM> covered_boxes = fine_level->getBoxes();
    covered_boxes.coarsen(fine_level->getRatioToCoarserLevel());
    int count = 0;
    for (PatchLevel<NDIM>::Iterator p(coarse_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = coarse_level->getPatch(p());
        Pointer<CellData<NDIM, double>> data = patch->getPatchData(cc_idx);
        for (int k = 0; k < covered_boxes.size(); ++k)
        {
            for (CellIterator<NDIM> ic(covered_boxes[k] * patch->getBox()); ic; ic++)
            {
                if (std::isnan((*data)(ic())))
                {
                    ++count;
                }
            }
        }
    }
    return IBTK_MPI::sumReduction(count);
}

// Return the largest difference between the coarse values and the reference values in the cells covered by the fine
// level.
double
max_covered_error(Pointer<PatchHierarchy<NDIM>> hierarchy,
                  const int idx,
                  const int exact_idx,
                  const bool side,
                  const Box<NDIM>* const region = nullptr)
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
            const Box<NDIM> box =
                region ? covered_boxes[k] * patch->getBox() * *region : covered_boxes[k] * patch->getBox();
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
    // Suppress the weighted-averaging warning, which contains the path of the source file.
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
    // Print the message of an error without file names or line numbers.
    Logger::getInstance()->setAbortAppender(Pointer<Logger::Appender>(new TestAppender()));
    TimerManager::createManager(nullptr);
    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "cubic_coarsen.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();
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
        Pointer<CellVariable<NDIM, double>> u_cc_var = new CellVariable<NDIM, double>("u_cc");
        Pointer<CellVariable<NDIM, double>> exact_cc_var = new CellVariable<NDIM, double>("exact_cc");
        Pointer<SideVariable<NDIM, double>> u_sc_var = new SideVariable<NDIM, double>("u_sc");
        Pointer<SideVariable<NDIM, double>> exact_sc_var = new SideVariable<NDIM, double>("exact_sc");
        std::vector<std::string> context_names = { "context" };
        tbox::Array<std::string> second_function_names;
        if (input_db->keyExists("second_function_names"))
        {
            context_names.push_back("second_context");
            second_function_names = input_db->getStringArray("second_function_names");
            TBOX_ASSERT(second_function_names.size() == function_names.size());
        }
        std::vector<Fields> fields;
        for (const std::string& context_name : context_names)
        {
            Pointer<VariableContext> ctx = var_db->getContext(context_name);
            fields.push_back({ var_db->registerVariableAndContext(u_cc_var, ctx, IntVector<NDIM>(1)),
                               var_db->registerVariableAndContext(exact_cc_var, ctx, IntVector<NDIM>(0)),
                               var_db->registerVariableAndContext(u_sc_var, ctx, IntVector<NDIM>(1)),
                               var_db->registerVariableAndContext(exact_sc_var, ctx, IntVector<NDIM>(0)) });
        }

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
        if (input_db->keyExists("fine_patch_boxes"))
        {
            // fine_patch_boxes contains the lower and upper corners of the required fine patches.
            const tbox::Array<int> bounds = input_db->getIntegerArray("fine_patch_boxes");
            Pointer<PatchLevel<NDIM>> fine_level = patch_hierarchy->getPatchLevel(1);
            const BoxArray<NDIM>& fine_boxes = fine_level->getBoxes();
            bool matches = bounds.size() == 2 * NDIM * fine_boxes.size();
            for (int k = 0; matches && k < fine_boxes.size(); ++k)
            {
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    matches = matches && bounds[2 * NDIM * k + d] == fine_boxes[k].lower(d) &&
                              bounds[2 * NDIM * k + NDIM + d] == fine_boxes[k].upper(d);
                }
            }
            if (!matches)
            {
                std::ostringstream actual;
                for (int k = 0; k < fine_boxes.size(); ++k)
                {
                    actual << " " << fine_boxes[k];
                }
                TBOX_ERROR("cubic_coarsen_01:\n"
                           << "  the fine patches differ from fine_patch_boxes\n"
                           << "  actual patches:" << actual.str() << "\n");
            }
        }
        for (int ln = 0; ln <= 1; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const Fields& f : fields)
            {
                level->allocatePatchData(f.u_cc_idx, 0.0);
                level->allocatePatchData(f.exact_cc_idx, 0.0);
                level->allocatePatchData(f.u_sc_idx, 0.0);
                level->allocatePatchData(f.exact_sc_idx, 0.0);
            }
        }

        using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
        std::vector<InterpolationTransactionComponent> transaction_comps;
        for (const Fields& f : fields)
        {
            transaction_comps.emplace_back(f.u_cc_idx, refine_op_name, false, "CUBIC_COARSEN", "NONE", false, nullptr);
            transaction_comps.emplace_back(f.u_sc_idx, refine_op_name, false, "CUBIC_COARSEN", "NONE", false, nullptr);
        }
        HierarchyGhostCellInterpolation ghost_fill;
        ghost_fill.initializeOperatorState(transaction_comps, patch_hierarchy);

        plog << std::scientific << std::setprecision(3);
        const bool nan_in_source = input_db->keyExists("nan_cell");
        Pointer<PatchLevel<NDIM>> coarse_level = patch_hierarchy->getPatchLevel(0);
        for (int k = 0; k < function_names.size(); ++k)
        {
            std::vector<std::string> names = { function_names[k] };
            if (fields.size() > 1)
            {
                names.push_back(second_function_names[k]);
            }
            for (unsigned int i = 0; i < fields.size(); ++i)
            {
                const Fields& f = fields[i];
                muParserCartGridFunction cc_fcn(
                    names[i] + "_cc", app_initializer->getComponentDatabase(names[i] + "_cc"), grid_geometry);
                muParserCartGridFunction sc_fcn(
                    names[i] + "_sc", app_initializer->getComponentDatabase(names[i] + "_sc"), grid_geometry);
                cc_fcn.setDataOnPatchHierarchy(f.u_cc_idx, u_cc_var, patch_hierarchy, 0.0);
                cc_fcn.setDataOnPatchHierarchy(f.exact_cc_idx, exact_cc_var, patch_hierarchy, 0.0);
                sc_fcn.setDataOnPatchHierarchy(f.u_sc_idx, u_sc_var, patch_hierarchy, 0.0);
                sc_fcn.setDataOnPatchHierarchy(f.exact_sc_idx, exact_sc_var, patch_hierarchy, 0.0);

                // Reset the coarse values: the coarsening must supply the covered ones.
                for (PatchLevel<NDIM>::Iterator p(coarse_level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = coarse_level->getPatch(p());
                    Pointer<CellData<NDIM, double>> cc_data = patch->getPatchData(f.u_cc_idx);
                    cc_data->fillAll(0.0);
                    Pointer<SideData<NDIM, double>> sc_data = patch->getPatchData(f.u_sc_idx);
                    sc_data->fillAll(0.0);
                }
                set_ghosts(patch_hierarchy->getPatchLevel(1),
                           f.u_cc_idx,
                           f.u_sc_idx,
                           std::numeric_limits<double>::quiet_NaN());
            }

            if (nan_in_source)
            {
                const tbox::Array<int> nan_cell = input_db->getIntegerArray("nan_cell");
                CellIndex<NDIM> index;
                for (unsigned int d = 0; d < NDIM; ++d)
                {
                    index(d) = nan_cell[d];
                }
                set_cell(patch_hierarchy->getPatchLevel(1),
                         fields[0].u_cc_idx,
                         index,
                         std::numeric_limits<double>::quiet_NaN());
            }

            ghost_fill.fillData(0.0);

            if (nan_in_source)
            {
                plog << names[0] << " covered coarse cells that are NaN = "
                     << count_covered_nans(patch_hierarchy, fields[0].u_cc_idx) << "\n";
                continue;
            }

            for (unsigned int i = 0; i < fields.size(); ++i)
            {
                const Fields& f = fields[i];
                const double cc_error = max_covered_error(patch_hierarchy, f.u_cc_idx, f.exact_cc_idx, false);
                const double sc_error = max_covered_error(patch_hierarchy, f.u_sc_idx, f.exact_sc_idx, true);
                if (!std::isfinite(cc_error) || !std::isfinite(sc_error))
                {
                    TBOX_ERROR("cubic_coarsen_01:\n"
                               << "  coarse values of function " << names[i] << " are not finite\n");
                }
                // Print roundoff-level errors as zero.
                plog << names[i] << " cell-centered max error = " << (cc_error < 1.0e-11 ? 0.0 : cc_error) << "\n";
                plog << names[i] << " side-centered max error = " << (sc_error < 1.0e-11 ? 0.0 : sc_error) << "\n";
                if (input_db->keyExists("corner_cells"))
                {
                    // The coarse cells at a reentrant corner or edge of the refined region.
                    const tbox::Array<int> bounds = input_db->getIntegerArray("corner_cells");
                    Box<NDIM> corner_cells;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        corner_cells.lower(d) = bounds[d];
                        corner_cells.upper(d) = bounds[NDIM + d];
                    }
                    const double cc_corner_error =
                        max_covered_error(patch_hierarchy, f.u_cc_idx, f.exact_cc_idx, false, &corner_cells);
                    const double sc_corner_error =
                        max_covered_error(patch_hierarchy, f.u_sc_idx, f.exact_sc_idx, true, &corner_cells);
                    plog << names[i] << " cell-centered max error in the corner cells = "
                         << (cc_corner_error < 1.0e-11 ? 0.0 : cc_corner_error) << "\n";
                    plog << names[i] << " side-centered max error in the corner cells = "
                         << (sc_corner_error < 1.0e-11 ? 0.0 : sc_corner_error) << "\n";
                }

                // Take the filled coarse values as the reference and fill again.
                for (PatchLevel<NDIM>::Iterator p(coarse_level); p; p++)
                {
                    Pointer<Patch<NDIM>> patch = coarse_level->getPatch(p());
                    Pointer<CellData<NDIM, double>> cc_data = patch->getPatchData(f.u_cc_idx);
                    Pointer<CellData<NDIM, double>> cc_exact = patch->getPatchData(f.exact_cc_idx);
                    cc_exact->copy(*cc_data);
                    Pointer<SideData<NDIM, double>> sc_data = patch->getPatchData(f.u_sc_idx);
                    Pointer<SideData<NDIM, double>> sc_exact = patch->getPatchData(f.exact_sc_idx);
                    sc_exact->copy(*sc_data);
                }
            }

            ghost_fill.fillData(0.0);

            for (unsigned int i = 0; i < fields.size(); ++i)
            {
                const Fields& f = fields[i];
                const double cc_change = max_covered_error(patch_hierarchy, f.u_cc_idx, f.exact_cc_idx, false);
                const double sc_change = max_covered_error(patch_hierarchy, f.u_sc_idx, f.exact_sc_idx, true);
                if (!std::isfinite(cc_change) || !std::isfinite(sc_change))
                {
                    TBOX_ERROR("cubic_coarsen_01:\n"
                               << "  coarse values of function " << names[i]
                               << " are not finite after the second fill\n");
                }
                plog << names[i]
                     << " cell-centered max change by second fill = " << (cc_change < 1.0e-11 ? 0.0 : cc_change)
                     << "\n";
                plog << names[i]
                     << " side-centered max change by second fill = " << (sc_change < 1.0e-11 ? 0.0 : sc_change)
                     << "\n";
            }
        }
    }
    return 0;
}
