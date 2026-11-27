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

#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/Logger.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>

#include <BergerRigoutsos.h>
#include <BoundaryBox.h>
#include <Box.h>
#include <BoxArray.h>
#include <BoxList.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <IntVector.h>
#include <LoadBalancer.h>
#include <Patch.h>
#include <PatchGeometry.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <RefineAlgorithm.h>
#include <RefineOperator.h>
#include <RefineSchedule.h>
#include <StandardTagAndInitStrategy.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Generate patches from tagged cells, regrid once, and check the patch layout and a ghost fill.

namespace
{
bool
in_domain(const hier::Index<NDIM>& cell, const BoxArray<NDIM>& domain)
{
    for (int k = 0; k < domain.size(); ++k)
    {
        if (domain[k].contains(cell))
        {
            return true;
        }
    }
    return false;
}

// Tag no cells. Inputs that list REFINE_BOXES also list GRADIENT_DETECTOR, so the gridding algorithm clusters the
// refine-box tags instead of using the boxes as patches.
class NoTags : public StandardTagAndInitStrategy<NDIM>
{
public:
    void initializeLevelData(Pointer<BasePatchHierarchy<NDIM>> /*hierarchy*/,
                             int /*level_number*/,
                             double /*init_data_time*/,
                             bool /*can_be_refined*/,
                             bool /*initial_time*/,
                             Pointer<BasePatchLevel<NDIM>> /*old_level*/,
                             bool /*allocate_data*/) override
    {
    }

    void resetHierarchyConfiguration(Pointer<BasePatchHierarchy<NDIM>> /*hierarchy*/,
                                     int /*coarsest_level*/,
                                     int /*finest_level*/) override
    {
    }

    void applyGradientDetector(Pointer<BasePatchHierarchy<NDIM>> /*hierarchy*/,
                               int /*level_number*/,
                               double /*error_data_time*/,
                               int /*tag_index*/,
                               bool /*initial_time*/,
                               bool /*uses_richardson_extrapolation_too*/) override
    {
    }
};
} // namespace

int
main(int argc, char** argv)
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<Logger::Appender> abort_append(new TestAppender());
    Logger::getInstance()->setAbortAppender(abort_append);
    Logger::getInstance()->setWarning(false);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    const IntVector<NDIM> ghosts(input_db->getInteger("ghost_width"));
    Pointer<CellVariable<NDIM, double>> var = new CellVariable<NDIM, double>("q");
    const int src_idx = var_db->registerVariableAndContext(var, var_db->getContext("src"), IntVector<NDIM>(0));
    const int dst_idx = var_db->registerVariableAndContext(var, var_db->getContext("dst"), ghosts);

    NoTags no_tags;
    Pointer<StandardTagAndInitialize<NDIM>> tagger = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", &no_tags, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    tagger,
                                    box_generator,
                                    load_balancer);

    gridding->makeCoarsestLevel(hierarchy, 0.0);
    for (int ln = 0; gridding->levelCanBeRefined(ln); ++ln)
    {
        gridding->makeFinerLevel(hierarchy, 0.0, true, 0);
        if (!hierarchy->finerLevelExists(ln))
        {
            break;
        }
    }
    // Check the layout after the initial gridding and again after one regrid.
    auto check_layout = [&](const char* stage)
    {
        plog << stage << ":\n";
        const double src_value = 1.0, not_set = -1.0e3;
        int too_narrow = 0, outside = 0, uncovered = 0, overlapping_pairs = 0;
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            BoxList<NDIM> domain(level->getPhysicalDomain());
            const BoxArray<NDIM>& boxes = level->getBoxes();
            plog << "level " << ln << ":";
            for (int k = 0; k < boxes.size(); ++k)
            {
                plog << ' ' << boxes[k];
                for (int d = 0; d < NDIM; ++d)
                {
                    if (boxes[k].numberCells(d) < ghosts(d))
                    {
                        ++too_narrow;
                    }
                }
                BoxList<NDIM> outside_cells(boxes[k]);
                outside_cells.removeIntersections(domain);
                outside += outside_cells.getTotalSizeOfBoxes();
                for (int j = k + 1; j < boxes.size(); ++j)
                {
                    if (boxes[k].intersects(boxes[j]))
                    {
                        ++overlapping_pairs;
                    }
                }
            }
            plog << '\n';

            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<PatchGeometry<NDIM>> patch_geometry = patch->getPatchGeometry();
                BoxList<NDIM> exterior(Box<NDIM>::grow(patch->getBox(), ghosts));
                exterior.removeIntersections(domain);
                for (int codim = 1; codim <= NDIM; ++codim)
                {
                    const tbox::Array<BoundaryBox<NDIM>>& bdry_boxes = patch_geometry->getCodimensionBoundaries(codim);
                    for (int k = 0; k < bdry_boxes.size(); ++k)
                    {
                        exterior.removeIntersections(
                            patch_geometry->getBoundaryFillBox(bdry_boxes[k], patch->getBox(), ghosts));
                    }
                }
                uncovered += exterior.getTotalSizeOfBoxes();
            }

            level->allocatePatchData(src_idx);
            level->allocatePatchData(dst_idx);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> src_data = patch->getPatchData(src_idx);
                Pointer<CellData<NDIM, double>> dst_data = patch->getPatchData(dst_idx);
                src_data->fillAll(src_value);
                dst_data->fillAll(not_set);
            }
        }
        plog << "patches narrower than the ghost width, by direction: " << too_narrow << '\n';
        plog << "patch cells outside the domain: " << outside << '\n';
        plog << "overlapping patch pairs: " << overlapping_pairs << '\n';
        plog << "ghost cells outside the domain that no boundary fill box covers: " << IBTK_MPI::sumReduction(uncovered)
             << '\n';

        Pointer<RefineOperator<NDIM>> refine_op = grid_geometry->lookupRefineOperator(var, "CONSTANT_REFINE");
        int cells_not_set = 0;
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            RefineAlgorithm<NDIM> refine_alg;
            refine_alg.registerRefine(dst_idx, src_idx, dst_idx, ln > 0 ? refine_op : Pointer<RefineOperator<NDIM>>());
            refine_alg.createSchedule(level, level, ln - 1, hierarchy)->fillData(0.0);

            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> dst_data = patch->getPatchData(dst_idx);
                for (CellIterator<NDIM> c(dst_data->getGhostBox()); c; c++)
                {
                    if (in_domain(c(), level->getPhysicalDomain()) && (*dst_data)(c()) == not_set)
                    {
                        ++cells_not_set;
                    }
                }
            }
        }
        plog << "cells in the domain the fill did not set: " << IBTK_MPI::sumReduction(cells_not_set) << '\n';
        for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
        {
            hierarchy->getPatchLevel(ln)->deallocatePatchData(src_idx);
            hierarchy->getPatchLevel(ln)->deallocatePatchData(dst_idx);
        }
    };
    check_layout("initial gridding");
    // Buffer the tags by one cell, so the regrid grows the layout.
    tbox::Array<int> tag_buffer(gridding->getMaxLevels());
    for (int ln = 0; ln < tag_buffer.size(); ++ln)
    {
        tag_buffer[ln] = 1;
    }
    gridding->regridAllFinerLevels(hierarchy, 0, 1.0, tag_buffer);
    check_layout("after one regrid");
}
