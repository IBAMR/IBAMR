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

#include <ibamr/PhaseChangeHierarchyIntegrator.h>
#include <ibamr/PhaseChangeUtilities.h>

#include "PhaseChangeTestUtilities.h"

#include <ibamr/app_namespaces.h>

RefinementRegion::RefinementRegion(const double dt) : d_dt(dt)
{
}

double
RefinementRegion::getCenter(const double time) const
{
    return (static_cast<int>(std::round(time / d_dt)) + static_cast<int>(d_swap_regions)) % 2 == 0 ? 0.25 : 0.75;
}

void
RefinementRegion::swapRegions()
{
    d_swap_regions = !d_swap_regions;
}

void
tag_moving_refinement_region(Pointer<BasePatchHierarchy<NDIM>> hierarchy,
                             const int level_number,
                             const double time,
                             int tag_idx,
                             bool,
                             bool,
                             void* ctx)
{
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(level_number);
    const int axis = NDIM - 1;
    const RefinementRegion& region = *static_cast<RefinementRegion*>(ctx);
    const double center = region.getCenter(time);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, int>> tags = patch->getPatchData(tag_idx);
        Pointer<CartesianPatchGeometry<NDIM>> geom = patch->getPatchGeometry();
        tags->fillAll(0);
        for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
        {
            const CellIndex<NDIM> ci(it());
            const double x =
                geom->getXLower()[axis] + (ci(axis) - patch->getBox().lower(axis) + 0.5) * geom->getDx()[axis];
            if (std::abs(x - center) < 0.125)
            {
                (*tags)(ci) = 1;
            }
        }
    }
}

void
check_restart_fields(Pointer<PatchHierarchy<NDIM>> hierarchy,
                     const std::vector<int>& cell_indices,
                     const std::vector<int>& side_indices,
                     const int step,
                     const double time,
                     const bool from_restart,
                     std::ostream& results)
{
    std::map<std::string, double> fields;
    fields["time"] = time;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            for (unsigned k = 0; k < cell_indices.size(); ++k)
            {
                Pointer<CellData<NDIM, double>> data = patch->getPatchData(cell_indices[k]);
                for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
                {
                    for (int depth = 0; depth < data->getDepth(); ++depth)
                    {
                        std::ostringstream key;
                        key << "cell/" << ln << '/' << k << '/' << depth << '/' << it();
                        fields[key.str()] = (*data)(CellIndex<NDIM>(it()), depth);
                    }
                }
            }
            for (unsigned k = 0; k < side_indices.size(); ++k)
            {
                Pointer<SideData<NDIM, double>> data = patch->getPatchData(side_indices[k]);
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    for (Box<NDIM>::Iterator it(SideGeometry<NDIM>::toSideBox(patch->getBox(), axis)); it; it++)
                    {
                        std::ostringstream key;
                        key << "side/" << ln << '/' << k << '/' << axis << '/' << it();
                        fields[key.str()] = (*data)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower));
                    }
                }
            }
        }
    }
    const std::string filename = "trajectory." + std::to_string(step) + "." + std::to_string(IBTK_MPI::getRank());
    std::ifstream reference;
    std::ofstream output;
    if (from_restart)
    {
        reference.open(filename);
        std::size_t size = 0;
        reference >> size;
        if (!reference || size != fields.size())
        {
            TBOX_ERROR("Restart trajectory has missing data or a different local mesh.\n");
        }
    }
    else
    {
        output.open(filename);
        output << fields.size() << '\n' << std::setprecision(17);
    }
    double max_error = 0.0;
    for (const std::pair<const std::string, double>& field : fields)
    {
        if (!std::isfinite(field.second))
        {
            TBOX_ERROR("Nonfinite phase-change field: " << field.first << "\n");
        }
        if (from_restart)
        {
            std::string key;
            double expected = 0.0;
            reference >> std::quoted(key) >> expected;
            if (!reference || key != field.first || !std::isfinite(expected))
            {
                TBOX_ERROR("Restart trajectory has invalid data or a different local mesh.\n");
            }
            max_error = std::max(max_error, std::abs(field.second - expected) / std::max(1.0, std::abs(expected)));
        }
        else
        {
            output << std::quoted(field.first) << ' ' << field.second << '\n';
        }
    }
    if (from_restart)
    {
        results << std::setprecision(13) << "Step " << step
                << " maximum scaled restart error = " << IBTK_MPI::maxReduction(max_error) << '\n';
    }
    if (!from_restart && !output)
    {
        TBOX_ERROR("Could not write uninterrupted trajectory.\n");
    }
}

void
check_liquid_fraction_tags(Pointer<AdvDiffHierarchyIntegrator> integrator, Pointer<PatchHierarchy<NDIM>> hierarchy)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellVariable<NDIM, double>> fraction_var = new CellVariable<NDIM, double>("tag_test_fraction");
    Pointer<CellVariable<NDIM, double>> gradient_var = new CellVariable<NDIM, double>("tag_test_gradient", NDIM);
    Pointer<CellVariable<NDIM, int>> tags_var = new CellVariable<NDIM, int>("tag_test_tags");
    const int fraction_idx = var_db->registerVariableAndContext(fraction_var, integrator->getCurrentContext());
    const int gradient_idx = var_db->registerVariableAndContext(gradient_var, integrator->getCurrentContext());
    const int tags_idx = var_db->registerVariableAndContext(tags_var, integrator->getCurrentContext());
    PhaseChangeUtilities::TagLiquidFractionRefinementCells tagger(integrator, fraction_var, gradient_var, 0.25, 0.75);
    const double fractions[] = { 0.0, 0.25, 0.5, 0.75, 1.0, 0.0 };
    int checked = 0;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (const int idx : { fraction_idx, gradient_idx, tags_idx })
        {
            level->allocatePatchData(idx);
        }
        for (const bool initial_time : { true, false })
        {
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, double>> fraction = patch->getPatchData(fraction_idx);
                Pointer<CellData<NDIM, double>> gradient = patch->getPatchData(gradient_idx);
                Pointer<CellData<NDIM, int>> tags = patch->getPatchData(tags_idx);
                gradient->fillAll(0.0);
                for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
                {
                    const CellIndex<NDIM> ci(it());
                    const int sample = ci(0) % 6;
                    (*fraction)(ci) = fractions[sample];
                    (*tags)(ci) = sample == 5 ? 1 : 0;
                    if (sample > 0 && sample <= NDIM)
                    {
                        (*gradient)(ci, sample - 1) = sample % 2 ? 2.0 : -2.0;
                    }
                }
            }
            PhaseChangeUtilities::call_tag_liquid_fraction_cells_callback(
                hierarchy, ln, 0.0, tags_idx, initial_time, false, &tagger);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<CellData<NDIM, int>> tags = patch->getPatchData(tags_idx);
                for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
                {
                    const CellIndex<NDIM> ci(it());
                    const int sample = ci(0) % 6;
                    const bool expected =
                        sample == 5 || (initial_time ? (sample >= 1 && sample <= 3) : (sample >= 1 && sample <= NDIM));
                    if ((*tags)(ci) != static_cast<int>(expected))
                    {
                        TBOX_ERROR("Incorrect liquid-fraction refinement tag at " << ci << "\n");
                    }
                    ++checked;
                }
            }
        }
        for (const int idx : { fraction_idx, gradient_idx, tags_idx })
        {
            level->deallocatePatchData(idx);
        }
    }
    if (IBTK_MPI::sumReduction(checked) == 0)
    {
        TBOX_ERROR("Liquid-fraction tagging test did not examine any cells.\n");
    }
    for (const int idx : { fraction_idx, gradient_idx, tags_idx })
    {
        var_db->removePatchDataIndex(idx);
    }
}

// On the finest level, interior patch cells have the same centered stencil as
// HierarchyMathOps. Check every gradient component against the live fraction.
void
check_liquid_fraction_gradient(Pointer<AdvDiffHierarchyIntegrator> integrator,
                               Pointer<PatchHierarchy<NDIM>> hierarchy,
                               Pointer<CellVariable<NDIM, double>> fraction_var,
                               Pointer<CellVariable<NDIM, double>> gradient_var,
                               std::ostream& results)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    const int fraction_idx = var_db->mapVariableAndContextToIndex(fraction_var, integrator->getCurrentContext());
    const int gradient_idx = var_db->mapVariableAndContextToIndex(gradient_var, integrator->getCurrentContext());
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(hierarchy->getFinestLevelNumber());
    int checked = 0;
    double max_error = 0.0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CartesianPatchGeometry<NDIM>> geom = patch->getPatchGeometry();
        Pointer<CellData<NDIM, double>> fraction = patch->getPatchData(fraction_idx);
        Pointer<CellData<NDIM, double>> gradient = patch->getPatchData(gradient_idx);
        Box<NDIM> interior = patch->getBox();
        interior.grow(-1);
        for (Box<NDIM>::Iterator it(interior); it; it++)
        {
            const CellIndex<NDIM> ci(it());
            for (int axis = 0; axis < NDIM; ++axis)
            {
                CellIndex<NDIM> left(ci), right(ci);
                --left(axis);
                ++right(axis);
                const double expected = ((*fraction)(right) - (*fraction)(left)) / (2.0 * geom->getDx()[axis]);
                const double actual = (*gradient)(ci, axis);
                if (!std::isfinite(actual) || !std::isfinite(expected))
                {
                    TBOX_ERROR("Incorrect stored liquid-fraction gradient in direction " << axis << "\n");
                }
                max_error = std::max(max_error, std::abs(actual - expected) / std::max(1.0, std::abs(expected)));
                ++checked;
            }
        }
    }
    results << std::setprecision(13) << "Maximum scaled fraction gradient error = " << IBTK_MPI::maxReduction(max_error)
            << '\n';
    if (IBTK_MPI::sumReduction(checked) == 0)
    {
        TBOX_ERROR("No interior liquid-fraction gradients were checked.\n");
    }
}

void
check_divergence_source_transfer(Pointer<HierarchyIntegrator> root_integrator,
                                 Pointer<PhaseChangeHierarchyIntegrator> phase_integrator,
                                 Pointer<PatchHierarchy<NDIM>> hierarchy,
                                 RefinementRegion& region,
                                 std::ostream& results)
{
    const int source_idx = phase_integrator->getVelocityDivergencePatchDataIndex();
    HierarchyCellDataOpsReal<NDIM, double> ops(hierarchy, 0, hierarchy->getFinestLevelNumber());
    // A constant is an exact transfer oracle on both the old and new mesh.
    ops.setToScalar(source_idx, 2.5);
    region.swapRegions();
    root_integrator->regridHierarchy();
    int checked = 0;
    double max_error = 0.0;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CellData<NDIM, double>> source = patch->getPatchData(source_idx);
            for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
            {
                const double value = (*source)(CellIndex<NDIM>(it()));
                if (!std::isfinite(value))
                {
                    TBOX_ERROR("Divergence source was not transferred to the new mesh.\n");
                }
                max_error = std::max(max_error, std::abs(value - 2.5));
                ++checked;
            }
        }
    }
    if (IBTK_MPI::sumReduction(checked) == 0)
    {
        TBOX_ERROR("No divergence-source cells were checked.\n");
    }
    results << std::setprecision(13)
            << "Maximum divergence source transfer error = " << IBTK_MPI::maxReduction(max_error) << '\n';
}
