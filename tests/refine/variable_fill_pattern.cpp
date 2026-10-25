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
#include <ibtk/CellNoCornersFillPattern.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/PIO.h>
#include <tbox/Pointer.h>

#include <Box.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <FirstLayerCellNoCornersFillPattern.h>
#include <Index.h>
#include <IntVector.h>
#include <LocallyActiveDataPatchLevelManager.h>
#include <LocallyActiveDataRefineAlgorithm.h>
#include <LocallyActiveVariableDatabase.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <VariableDatabase.h>
#include <VariableFillPattern.h>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <ibtk/app_namespaces.h>

// Report which cells refine schedules fill for refine items with different
// fill patterns. Cell values are set if they hold the exact value of a
// function that is linear in space and time.
// Ghost cells inside the domain are counted as face ghosts (outside the patch
// in one direction) or corner ghosts (outside in two or more).

namespace
{
const double NOT_SET = -1.0e3;
const int RATIO = 2;

double
exact(const hier::Patch<NDIM>& patch, const hier::Index<NDIM>& i, const double time)
{
    tbox::Pointer<geom::CartesianPatchGeometry<NDIM>> patch_geom = patch.getPatchGeometry();
    const double* const dx = patch_geom->getDx();
    const double* const x_lower = patch_geom->getXLower();
    double value = 1.0 + time;
    for (int d = 0; d < NDIM; ++d)
    {
        value += (d + 1) * (x_lower[d] + (i(d) - patch.getBox().lower(d) + 0.5) * dx[d]);
    }
    return value;
}

// Set the interior to the exact values and the ghost cells to NOT_SET.
void
set_data(tbox::Pointer<hier::PatchLevel<NDIM>> level, const int idx, const double time)
{
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        tbox::Pointer<pdat::CellData<NDIM, double>> data = patch->getPatchData(idx);
        data->fillAll(NOT_SET);
        for (pdat::CellIterator<NDIM> i(patch->getBox()); i; i++)
        {
            (*data)(i()) = exact(*patch, i(), time);
        }
        data->setTime(time);
    }
}

void
report(const std::string& label, tbox::Pointer<hier::PatchLevel<NDIM>> level, const int idx, const double time)
{
    const hier::Box<NDIM> domain = level->getPhysicalDomain()[0];
    int face = 0, face_set = 0, corner = 0, corner_set = 0;
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        tbox::Pointer<pdat::CellData<NDIM, double>> data = patch->getPatchData(idx);
        const hier::Box<NDIM>& box = patch->getBox();
        for (pdat::CellIterator<NDIM> i(data->getGhostBox()); i; i++)
        {
            if (box.contains(i()) || !domain.contains(i()))
            {
                continue;
            }
            int n_outside = 0;
            for (int d = 0; d < NDIM; ++d)
            {
                if (i()(d) < box.lower(d) || i()(d) > box.upper(d))
                {
                    ++n_outside;
                }
            }
            const bool set = std::abs((*data)(i()) - exact(*patch, i(), time)) < 1.0e-10;
            if (n_outside == 1)
            {
                ++face;
                face_set += set;
            }
            else
            {
                ++corner;
                corner_set += set;
            }
        }
    }
    face = IBTK::IBTK_MPI::sumReduction(face);
    face_set = IBTK::IBTK_MPI::sumReduction(face_set);
    corner = IBTK::IBTK_MPI::sumReduction(corner);
    corner_set = IBTK::IBTK_MPI::sumReduction(corner_set);
    tbox::plog << label << ": face ghosts set " << face_set << " of " << face << ", corner ghosts set " << corner_set
               << " of " << corner << '\n';
}

// On level 0, the lower half in the first direction has 2^(NDIM-1) patches and
// the upper half is one patch, so that face neighbors supply ghost cells at
// corners, and in 3D diagonal neighbors do too. Level 1, if present, covers
// the middle of the domain away from the physical boundary, so its ghost cells
// lie on the coarse-fine interface.
tbox::Pointer<hier::PatchHierarchy<NDIM>>
make_hierarchy(tbox::Pointer<geom::CartesianGridGeometry<NDIM>> grid_geometry,
               const int n_levels,
               const std::vector<int>& idxs)
{
    tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy =
        new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    const hier::Box<NDIM> domain = grid_geometry->getPhysicalDomain()[0];
    const int n = domain.numberCells(0) / 2;
    const int n_lower = 1 << (NDIM - 1);
    hier::BoxArray<NDIM> coarse_boxes(n_lower + 1);
    for (int k = 0; k < n_lower; ++k)
    {
        hier::Index<NDIM> lower(0);
        for (int d = 1; d < NDIM; ++d)
        {
            lower(d) = ((k >> (d - 1)) & 1) * n;
        }
        coarse_boxes[k] = hier::Box<NDIM>(lower, lower + (n - 1));
    }
    hier::Index<NDIM> upper_lower(0);
    upper_lower(0) = n;
    coarse_boxes[n_lower] = hier::Box<NDIM>(upper_lower, domain.upper());
    hier::BoxArray<NDIM> fine_boxes(1);
    fine_boxes[0] = hier::Box<NDIM>(hier::Index<NDIM>(RATIO * n / 2), hier::Index<NDIM>(RATIO * 3 * n / 2 - 1));
    const std::vector<hier::BoxArray<NDIM>> boxes = { coarse_boxes, fine_boxes };
    for (int ln = 0; ln < n_levels; ++ln)
    {
        hier::ProcessorMapping mapping(boxes[ln].size());
        for (int k = 0; k < boxes[ln].size(); ++k)
        {
            mapping.setProcessorAssignment(k, k % IBTK::IBTK_MPI::getNodes());
        }
        hierarchy->makeNewPatchLevel(ln, hier::IntVector<NDIM>(ln == 0 ? 1 : RATIO), boxes[ln], mapping);
        for (const int idx : idxs)
        {
            hierarchy->getPatchLevel(ln)->allocatePatchData(idx, 0.0);
        }
    }
    return hierarchy;
}
} // namespace

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    tbox::Pointer<IBTK::AppInitializer> app_initializer = new IBTK::AppInitializer(argc, argv);
    tbox::Pointer<geom::CartesianGridGeometry<NDIM>> grid_geometry = new geom::CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

    auto* var_db = hier::VariableDatabase<NDIM>::getDatabase();
    tbox::Pointer<pdat::CellVariable<NDIM, double>> var = new pdat::CellVariable<NDIM, double>("u");
    const int u_idx = var_db->registerVariableAndContext(var, var_db->getContext("u"), hier::IntVector<NDIM>(1));
    const int v_idx = var_db->registerVariableAndContext(var, var_db->getContext("v"), hier::IntVector<NDIM>(1));
    const int narrow_idx =
        var_db->registerVariableAndContext(var, var_db->getContext("narrow"), hier::IntVector<NDIM>(2));
    const int wide_idx = var_db->registerVariableAndContext(var, var_db->getContext("wide"), hier::IntVector<NDIM>(2));
    const int old_idx = var_db->registerVariableAndContext(var, var_db->getContext("old"), hier::IntVector<NDIM>(0));
    const int new_idx = var_db->registerVariableAndContext(var, var_db->getContext("new"), hier::IntVector<NDIM>(0));
    tbox::Pointer<xfer::RefineOperator<NDIM>> refine_op =
        grid_geometry->lookupRefineOperator(var, "CONSERVATIVE_LINEAR_REFINE");
    tbox::Pointer<xfer::TimeInterpolateOperator<NDIM>> time_op =
        grid_geometry->lookupTimeInterpolateOperator(var, "STD_LINEAR_TIME_INTERPOLATE");
    auto no_corners = []() -> tbox::Pointer<xfer::VariableFillPattern<NDIM>>
    { return new IBTK::CellNoCornersFillPattern(1, /*overwrite_interior*/ false); };

    // Two items of one patch data type and ghost width on one level, one with
    // the default pattern and one without corners, registered in both orders.
    // Without refine operators, the item registered first heads a shared
    // class.
    {
        tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy = make_hierarchy(grid_geometry, 1, { u_idx, v_idx });
        tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        for (const bool default_first : { true, false })
        {
            xfer::RefineAlgorithm<NDIM> refine_alg;
            if (default_first)
            {
                refine_alg.registerRefine(u_idx, u_idx, u_idx, nullptr);
            }
            refine_alg.registerRefine(v_idx, v_idx, v_idx, nullptr, no_corners());
            if (!default_first)
            {
                refine_alg.registerRefine(u_idx, u_idx, u_idx, nullptr);
            }
            set_data(level, u_idx, 0.0);
            set_data(level, v_idx, 0.0);
            refine_alg.createSchedule(level)->fillData(0.0);
            const std::string order = default_first ? "default registered first" : "no corners registered first";
            report("one level, " + order + ", default", level, u_idx, 0.0);
            report("one level, " + order + ", no corners", level, v_idx, 0.0);
        }
    }

    // Two items on one level without corners, with stencil widths 1 and 2.
    // Without refine operators, the item registered first heads a shared
    // class.
    {
        tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy =
            make_hierarchy(grid_geometry, 1, { narrow_idx, wide_idx });
        tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        xfer::RefineAlgorithm<NDIM> refine_alg;
        refine_alg.registerRefine(
            narrow_idx, narrow_idx, narrow_idx, nullptr, new IBTK::CellNoCornersFillPattern(1, false));
        refine_alg.registerRefine(wide_idx, wide_idx, wide_idx, nullptr, new IBTK::CellNoCornersFillPattern(2, false));
        set_data(level, narrow_idx, 0.0);
        set_data(level, wide_idx, 0.0);
        refine_alg.createSchedule(level)->fillData(0.0);
        report("one level, no corners, stencil width 1", level, narrow_idx, 0.0);
        report("one level, no corners, stencil width 2", level, wide_idx, 0.0);
    }

    // A time-interpolated item without corners.
    {
        tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy =
            make_hierarchy(grid_geometry, 1, { u_idx, old_idx, new_idx });
        tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        xfer::RefineAlgorithm<NDIM> refine_alg;
        refine_alg.registerRefine(u_idx, new_idx, old_idx, new_idx, u_idx, nullptr, time_op, no_corners());
        set_data(level, u_idx, 0.5);
        set_data(level, old_idx, 0.0);
        set_data(level, new_idx, 1.0);
        refine_alg.createSchedule(level, level, nullptr, /*use_time_interpolation*/ true)->fillData(0.5);
        report("one level, time interpolated, no corners", level, u_idx, 0.5);
    }

    // Fills of a refined level, each item on its own. No-corners patterns
    // still fill corner ghost cells on the coarse-fine interface.
    {
        tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy = make_hierarchy(grid_geometry, 2, { u_idx });
        const std::vector<std::pair<std::string, tbox::Pointer<xfer::VariableFillPattern<NDIM>>>> patterns = {
            { "default", nullptr },
            { "IBTK no corners", no_corners() },
            { "SAMRAI first layer no corners", new pdat::FirstLayerCellNoCornersFillPattern<NDIM>() }
        };
        for (const auto& pattern : patterns)
        {
            xfer::RefineAlgorithm<NDIM> refine_alg;
            refine_alg.registerRefine(u_idx, u_idx, u_idx, refine_op, pattern.second);
            set_data(hierarchy->getPatchLevel(0), u_idx, 0.0);
            set_data(hierarchy->getPatchLevel(1), u_idx, 0.0);
            refine_alg.createSchedule(hierarchy->getPatchLevel(1), 0, hierarchy)->fillData(0.0);
            report("two levels, " + pattern.first, hierarchy->getPatchLevel(1), u_idx, 0.0);
        }
    }

    // A time-interpolated fill of two equivalent items with the locally-active
    // classes, which do not set fill patterns.
    {
        auto* active_var_db = hier::LocallyActiveVariableDatabase<NDIM>::getDatabase();
        const int active_idx =
            active_var_db->registerVariable(new pdat::CellVariable<NDIM, double>("active"), hier::IntVector<NDIM>(1));
        const int active_2_idx =
            active_var_db->registerVariable(new pdat::CellVariable<NDIM, double>("active 2"), hier::IntVector<NDIM>(1));
        const int active_old_idx = active_var_db->registerVariable(new pdat::CellVariable<NDIM, double>("active::old"),
                                                                   hier::IntVector<NDIM>(0));
        const int active_new_idx = active_var_db->registerVariable(new pdat::CellVariable<NDIM, double>("active::new"),
                                                                   hier::IntVector<NDIM>(0));
        tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy = make_hierarchy(grid_geometry, 1, {});
        tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        tbox::Pointer<hier::LocallyActiveDataPatchLevelManager<NDIM>> level_mgr =
            active_var_db->getLocallyActiveDataPatchLevelManager(level);
        const std::vector<std::pair<int, double>> idxs_and_times = {
            { active_idx, 0.5 }, { active_2_idx, 0.5 }, { active_old_idx, 0.0 }, { active_new_idx, 1.0 }
        };
        for (const auto& idx_and_time : idxs_and_times)
        {
            for (int k = 0; k < level->getNumberOfPatches(); ++k)
            {
                level_mgr->setPatchDataActive(hier::PatchDataId(idx_and_time.first), hier::PatchNumber(k));
            }
            level_mgr->allocatePatchData(hier::PatchDataId(idx_and_time.first), idx_and_time.second);
            set_data(level, idx_and_time.first, idx_and_time.second);
        }
        xfer::LocallyActiveDataRefineAlgorithm<NDIM> refine_alg;
        refine_alg.registerRefine(
            active_idx, active_new_idx, active_old_idx, active_new_idx, active_idx, nullptr, time_op);
        refine_alg.registerRefine(
            active_2_idx, active_new_idx, active_old_idx, active_new_idx, active_2_idx, nullptr, time_op);
        refine_alg.createSchedule(level, level_mgr, level, level_mgr, nullptr, /*use_time_interpolation*/ true)
            ->fillData(0.5, /*do_physical_boundary_fill*/ false);
        report("one level, locally active, time interpolated", level, active_idx, 0.5);
        report("one level, locally active, time interpolated, second item", level, active_2_idx, 0.5);
    }
}
