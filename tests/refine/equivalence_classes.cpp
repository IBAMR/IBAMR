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

#include <tbox/PIO.h>
#include <tbox/Pointer.h>

#include <Box.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <CoarsenAlgorithm.h>
#include <CoarsenPatchStrategy.h>
#include <CoarsenSchedule.h>
#include <Index.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIterator.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <cmath>
#include <string>
#include <vector>

#include <ibtk/app_namespaces.h>

// Report which values refine and coarsen schedules set for items that could
// wrongly share an equivalence class: refine items whose scratch data differ in
// ghost width, side data items with different directions, and coarsen items
// that coarsen different ghost cell widths. Cell values are set if they hold
// the exact value of a function that is linear in space; side values are set
// if they hold 1 plus their axis.

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
report_interior(const std::string& label, tbox::Pointer<hier::PatchLevel<NDIM>> level, const int idx, const double time)
{
    int cells = 0, cells_set = 0;
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        tbox::Pointer<pdat::CellData<NDIM, double>> data = patch->getPatchData(idx);
        for (pdat::CellIterator<NDIM> i(patch->getBox()); i; i++)
        {
            ++cells;
            cells_set += std::abs((*data)(i()) - exact(*patch, i(), time)) < 1.0e-10;
        }
    }
    cells = IBTK::IBTK_MPI::sumReduction(cells);
    cells_set = IBTK::IBTK_MPI::sumReduction(cells_set);
    tbox::plog << label << ": interior cells set " << cells_set << " of " << cells << '\n';
}

// Side values hold 1 plus their axis.
void
set_side_data(tbox::Pointer<hier::PatchLevel<NDIM>> level, const int idx)
{
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        tbox::Pointer<pdat::SideData<NDIM, double>> data = patch->getPatchData(idx);
        data->fillAll(NOT_SET);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            if (!data->getDirectionVector()(axis))
            {
                continue;
            }
            for (pdat::SideIterator<NDIM> i(patch->getBox(), axis); i; i++)
            {
                (*data)(i()) = 1.0 + axis;
            }
        }
    }
}

void
report_side(const std::string& label, tbox::Pointer<hier::PatchLevel<NDIM>> level, const int idx)
{
    const hier::Box<NDIM> domain = level->getPhysicalDomain()[0];
    tbox::plog << label << ": ghost values set";
    for (int axis = 0; axis < NDIM; ++axis)
    {
        int values = 0, values_set = 0;
        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
            tbox::Pointer<pdat::SideData<NDIM, double>> data = patch->getPatchData(idx);
            const hier::Box<NDIM> side_box = pdat::SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
            const hier::Box<NDIM> domain_side_box = pdat::SideGeometry<NDIM>::toSideBox(domain, axis);
            for (pdat::SideIterator<NDIM> i(data->getGhostBox(), axis); i; i++)
            {
                if (side_box.contains(i()) || !domain_side_box.contains(i()))
                {
                    continue;
                }
                ++values;
                values_set += (*data)(i()) == 1.0 + axis;
            }
        }
        values = IBTK::IBTK_MPI::sumReduction(values);
        values_set = IBTK::IBTK_MPI::sumReduction(values_set);
        tbox::plog << (axis == 0 ? " " : ", ") << "axis " << axis << ' ' << values_set << " of " << values;
    }
    tbox::plog << '\n';
}

// Coarsen operators set only the interior of the coarsened data. Set the ring
// of width 1 around it to a value.
class SetRing : public xfer::CoarsenPatchStrategy<NDIM>
{
public:
    SetRing(const int idx, const double value) : d_idx(idx), d_value(value)
    {
    }

    hier::IntVector<NDIM> getCoarsenOpStencilWidth() const override
    {
        return hier::IntVector<NDIM>(0);
    }

    void preprocessCoarsen(hier::Patch<NDIM>& /*coarse*/,
                           const hier::Patch<NDIM>& /*fine*/,
                           const hier::Box<NDIM>& /*coarse_box*/,
                           const hier::IntVector<NDIM>& /*ratio*/) override
    {
    }

    void postprocessCoarsen(hier::Patch<NDIM>& coarse,
                            const hier::Patch<NDIM>& /*fine*/,
                            const hier::Box<NDIM>& coarse_box,
                            const hier::IntVector<NDIM>& /*ratio*/) override
    {
        tbox::Pointer<pdat::CellData<NDIM, double>> data = coarse.getPatchData(d_idx);
        for (pdat::CellIterator<NDIM> i(hier::Box<NDIM>::grow(coarse_box, hier::IntVector<NDIM>(1))); i; i++)
        {
            if (!coarse_box.contains(i()))
            {
                (*data)(i()) = d_value;
            }
        }
    }

private:
    const int d_idx;
    const double d_value;
};

// Count the cells of level 0 in the ring of width 1 around level 1 that hold
// value, and the cells covered by level 1 that hold covered_value.
void
report_ring(const std::string& label,
            tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy,
            const int idx,
            const double value,
            const double covered_value)
{
    tbox::Pointer<hier::PatchLevel<NDIM>> fine_level = hierarchy->getPatchLevel(1);
    hier::Box<NDIM> covered = fine_level->getBoxes()[0];
    covered.coarsen(hier::IntVector<NDIM>(RATIO));
    const hier::Box<NDIM> grown = hier::Box<NDIM>::grow(covered, hier::IntVector<NDIM>(1));
    tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    int ring = 0, ring_set = 0, inside = 0, inside_set = 0;
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        tbox::Pointer<pdat::CellData<NDIM, double>> data = patch->getPatchData(idx);
        for (pdat::CellIterator<NDIM> i(patch->getBox()); i; i++)
        {
            if (covered.contains(i()))
            {
                ++inside;
                inside_set += (*data)(i()) == covered_value;
            }
            else if (grown.contains(i()))
            {
                ++ring;
                ring_set += (*data)(i()) == value;
            }
        }
    }
    ring = IBTK::IBTK_MPI::sumReduction(ring);
    ring_set = IBTK::IBTK_MPI::sumReduction(ring_set);
    inside = IBTK::IBTK_MPI::sumReduction(inside);
    inside_set = IBTK::IBTK_MPI::sumReduction(inside_set);
    tbox::plog << label << ": ring cells holding " << value << ' ' << ring_set << " of " << ring
               << ", covered cells holding " << covered_value << ' ' << inside_set << " of " << inside << '\n';
}

// Count, for each axis, the side values of level 0 strictly inside the region
// covered by level 1 that hold 1 plus their axis.
void
report_coarse_side(const std::string& label, tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy, const int idx)
{
    tbox::Pointer<hier::PatchLevel<NDIM>> fine_level = hierarchy->getPatchLevel(1);
    hier::Box<NDIM> covered = fine_level->getBoxes()[0];
    covered.coarsen(hier::IntVector<NDIM>(RATIO));
    tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    tbox::plog << label << ": covered values set";
    for (int axis = 0; axis < NDIM; ++axis)
    {
        int values = 0, values_set = 0;
        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
            tbox::Pointer<pdat::SideData<NDIM, double>> data = patch->getPatchData(idx);
            for (pdat::SideIterator<NDIM> i(patch->getBox() * covered, axis); i; i++)
            {
                if (i()(axis) == covered.lower(axis) || i()(axis) == covered.upper(axis) + 1)
                {
                    continue;
                }
                ++values;
                values_set += (*data)(i()) == 1.0 + axis;
            }
        }
        values = IBTK::IBTK_MPI::sumReduction(values);
        values_set = IBTK::IBTK_MPI::sumReduction(values_set);
        tbox::plog << (axis == 0 ? " " : ", ") << "axis " << axis << ' ' << values_set << " of " << values;
    }
    tbox::plog << '\n';
}

// Count the cells of the ghost boxes inside the domain that hold the exact
// value, including the patch interiors.
void
report_cells(const std::string& label, tbox::Pointer<hier::PatchLevel<NDIM>> level, const int idx)
{
    const hier::Box<NDIM> domain = level->getPhysicalDomain()[0];
    int cells = 0, cells_set = 0;
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        tbox::Pointer<pdat::CellData<NDIM, double>> data = patch->getPatchData(idx);
        for (pdat::CellIterator<NDIM> i(data->getGhostBox() * domain); i; i++)
        {
            ++cells;
            cells_set += std::abs((*data)(i()) - exact(*patch, i(), 0.0)) < 1.0e-10;
        }
    }
    cells = IBTK::IBTK_MPI::sumReduction(cells);
    cells_set = IBTK::IBTK_MPI::sumReduction(cells_set);
    tbox::plog << label << ": cells set " << cells_set << " of " << cells << '\n';
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
    const int x_idx = var_db->registerVariableAndContext(var, var_db->getContext("x"), hier::IntVector<NDIM>(0));
    const int y_idx = var_db->registerVariableAndContext(var, var_db->getContext("y"), hier::IntVector<NDIM>(0));
    const int x_scratch_idx =
        var_db->registerVariableAndContext(var, var_db->getContext("x::scratch"), hier::IntVector<NDIM>(0));
    const int y_scratch_idx =
        var_db->registerVariableAndContext(var, var_db->getContext("y::scratch"), hier::IntVector<NDIM>(1));
    tbox::Pointer<pdat::SideVariable<NDIM, double>> side_var = new pdat::SideVariable<NDIM, double>("s");
    tbox::Pointer<pdat::SideVariable<NDIM, double>> axis_0_side_var =
        new pdat::SideVariable<NDIM, double>("s::axis 0", 1, true, /*direction*/ 0);
    const int side_idx =
        var_db->registerVariableAndContext(side_var, var_db->getContext("s"), hier::IntVector<NDIM>(1));
    const int axis_0_side_idx =
        var_db->registerVariableAndContext(axis_0_side_var, var_db->getContext("s"), hier::IntVector<NDIM>(1));
    const int ring_idx =
        var_db->registerVariableAndContext(var, var_db->getContext("ring"), hier::IntVector<NDIM>(RATIO));
    const int no_ring_idx =
        var_db->registerVariableAndContext(var, var_db->getContext("no ring"), hier::IntVector<NDIM>(RATIO));
    const int copy_idx = var_db->registerVariableAndContext(var, var_db->getContext("copy"), hier::IntVector<NDIM>(1));
    const int copy_source_idx =
        var_db->registerVariableAndContext(var, var_db->getContext("copy source"), hier::IntVector<NDIM>(1));
    tbox::Pointer<xfer::RefineOperator<NDIM>> refine_op =
        grid_geometry->lookupRefineOperator(var, "CONSERVATIVE_LINEAR_REFINE");

    // Two items whose scratch data differ in ghost width and whose refine
    // operators differ, filling a refined level from level 0 alone, registered
    // in both orders. The item registered last heads a shared class, since
    // both have refine operators. Report the interior of the item with linear
    // refinement.
    {
        tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy = make_hierarchy(grid_geometry, 2, { x_idx, y_idx });
        tbox::Pointer<hier::PatchLevel<NDIM>> fine_level = hierarchy->getPatchLevel(1);
        tbox::Pointer<xfer::RefineOperator<NDIM>> constant_refine_op =
            grid_geometry->lookupRefineOperator(var, "CONSTANT_REFINE");
        for (const bool linear_first : { true, false })
        {
            xfer::RefineAlgorithm<NDIM> refine_alg;
            if (linear_first)
            {
                refine_alg.registerRefine(y_idx, y_idx, y_scratch_idx, refine_op);
            }
            refine_alg.registerRefine(x_idx, x_idx, x_scratch_idx, constant_refine_op);
            if (!linear_first)
            {
                refine_alg.registerRefine(y_idx, y_idx, y_scratch_idx, refine_op);
            }
            for (const int idx : { x_idx, y_idx })
            {
                set_data(hierarchy->getPatchLevel(0), idx, 0.0);
                for (hier::PatchLevel<NDIM>::Iterator p(fine_level); p; p++)
                {
                    tbox::Pointer<pdat::CellData<NDIM, double>> data = fine_level->getPatch(p())->getPatchData(idx);
                    data->fillAll(NOT_SET);
                }
            }
            refine_alg.createSchedule(fine_level, tbox::Pointer<hier::PatchLevel<NDIM>>(), 0, hierarchy)->fillData(0.0);
            const std::string order =
                linear_first ? "linear refine registered first" : "constant refine registered first";
            report_interior("two levels, from level 0 only, " + order, fine_level, y_idx, 0.0);
        }
    }

    // Two side data items on one level, one with all directions and one with
    // axis 0 only, registered in both orders. The item registered first heads
    // a shared class, since neither has a refine operator. Report the item with
    // all directions.
    {
        tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy =
            make_hierarchy(grid_geometry, 1, { side_idx, axis_0_side_idx });
        tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        for (const bool all_first : { true, false })
        {
            xfer::RefineAlgorithm<NDIM> refine_alg;
            if (all_first)
            {
                refine_alg.registerRefine(side_idx, side_idx, side_idx, nullptr);
            }
            refine_alg.registerRefine(axis_0_side_idx, axis_0_side_idx, axis_0_side_idx, nullptr);
            if (!all_first)
            {
                refine_alg.registerRefine(side_idx, side_idx, side_idx, nullptr);
            }
            set_side_data(level, side_idx);
            set_side_data(level, axis_0_side_idx);
            refine_alg.createSchedule(level)->fillData(0.0);
            const std::string order = all_first ? "all directions registered first" : "axis 0 registered first";
            report_side("one level, side data, " + order, level, side_idx);
        }
    }

    // Two cell items coarsened from level 1 onto level 0, one of them also in
    // a ring of one ghost cell around level 1, registered in both orders. The
    // item registered last heads a shared class, since both have coarsen
    // operators. Level 1 holds 5, level 0 holds 3, and the patch strategy sets
    // the ring of the item that has one to 7.
    {
        tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy =
            make_hierarchy(grid_geometry, 2, { ring_idx, no_ring_idx });
        tbox::Pointer<xfer::CoarsenOperator<NDIM>> coarsen_op =
            grid_geometry->lookupCoarsenOperator(var, "CONSERVATIVE_COARSEN");
        for (const bool ring_first : { true, false })
        {
            xfer::CoarsenAlgorithm<NDIM> coarsen_alg;
            if (ring_first)
            {
                coarsen_alg.registerCoarsen(ring_idx, ring_idx, coarsen_op, hier::IntVector<NDIM>(1));
            }
            coarsen_alg.registerCoarsen(no_ring_idx, no_ring_idx, coarsen_op);
            if (!ring_first)
            {
                coarsen_alg.registerCoarsen(ring_idx, ring_idx, coarsen_op, hier::IntVector<NDIM>(1));
            }
            for (int ln = 0; ln < 2; ++ln)
            {
                tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
                for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
                {
                    for (const int idx : { ring_idx, no_ring_idx })
                    {
                        tbox::Pointer<pdat::CellData<NDIM, double>> data = level->getPatch(p())->getPatchData(idx);
                        data->fillAll(ln == 0 ? 3.0 : 5.0);
                    }
                }
            }
            SetRing set_ring(ring_idx, 7.0);
            coarsen_alg.createSchedule(hierarchy->getPatchLevel(0), hierarchy->getPatchLevel(1), &set_ring)
                ->coarsenData();
            const std::string order = ring_first ? "ring registered first" : "no ring registered first";
            report_ring("coarsen, " + order + ", item with ring", hierarchy, ring_idx, 7.0, 5.0);
            report_ring("coarsen, " + order + ", item without ring", hierarchy, no_ring_idx, 3.0, 5.0);
        }
    }

    // Two side data items coarsened from level 1 onto level 0, one with all
    // directions and one with axis 0 only, registered in both orders. The item
    // registered last heads a shared class, since both have coarsen operators.
    // Level 1 holds 1 plus the axis and level 0 holds 0. Report the item with
    // all directions.
    {
        tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy =
            make_hierarchy(grid_geometry, 2, { side_idx, axis_0_side_idx });
        tbox::Pointer<xfer::CoarsenOperator<NDIM>> side_coarsen_op =
            grid_geometry->lookupCoarsenOperator(side_var, "CONSERVATIVE_COARSEN");
        for (const bool all_first : { true, false })
        {
            xfer::CoarsenAlgorithm<NDIM> coarsen_alg;
            if (all_first)
            {
                coarsen_alg.registerCoarsen(side_idx, side_idx, side_coarsen_op);
            }
            coarsen_alg.registerCoarsen(axis_0_side_idx, axis_0_side_idx, side_coarsen_op);
            if (!all_first)
            {
                coarsen_alg.registerCoarsen(side_idx, side_idx, side_coarsen_op);
            }
            for (const int idx : { side_idx, axis_0_side_idx })
            {
                set_side_data(hierarchy->getPatchLevel(1), idx);
                tbox::Pointer<hier::PatchLevel<NDIM>> coarse_level = hierarchy->getPatchLevel(0);
                for (hier::PatchLevel<NDIM>::Iterator p(coarse_level); p; p++)
                {
                    tbox::Pointer<pdat::SideData<NDIM, double>> data = coarse_level->getPatch(p())->getPatchData(idx);
                    data->fillAll(0.0);
                }
            }
            coarsen_alg.createSchedule(hierarchy->getPatchLevel(0), hierarchy->getPatchLevel(1))->coarsenData();
            const std::string order = all_first ? "all directions registered first" : "axis 0 registered first";
            report_coarse_side("coarsen, side data, " + order, hierarchy, side_idx);
        }
    }

    // A schedule built for an item whose source is other data, which copies
    // from each patch to itself, reset to an item whose source is its scratch
    // data. Before the reset, the ghost values and the old source are unset
    // again, so only an in-place fill after the reset sets the ghost values.
    {
        tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy =
            make_hierarchy(grid_geometry, 1, { copy_idx, copy_source_idx });
        tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
        xfer::RefineAlgorithm<NDIM> copy_alg;
        copy_alg.registerRefine(copy_idx, copy_source_idx, copy_idx, nullptr);
        set_data(level, copy_idx, 0.0);
        set_data(level, copy_source_idx, 0.0);
        tbox::Pointer<xfer::RefineSchedule<NDIM>> schedule = copy_alg.createSchedule(level);
        schedule->fillData(0.0);
        set_data(level, copy_idx, 0.0);
        for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            tbox::Pointer<pdat::CellData<NDIM, double>> data = level->getPatch(p())->getPatchData(copy_source_idx);
            data->fillAll(NOT_SET);
        }
        xfer::RefineAlgorithm<NDIM> in_place_alg;
        in_place_alg.registerRefine(copy_idx, copy_idx, copy_idx, nullptr);
        in_place_alg.resetSchedule(schedule);
        schedule->fillData(0.0);
        report_cells("one level, copying schedule reset to in place", level, copy_idx);
    }
}
