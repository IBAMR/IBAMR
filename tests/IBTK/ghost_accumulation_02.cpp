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
#include <ibtk/CartCellRobinPhysBdryOp.h>
#include <ibtk/CartSideRobinPhysBdryOp.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/SAMRAIGhostDataAccumulator.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ArrayData.h>
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <SAMRAI_config.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <memory>
#include <vector>

#include <ibtk/app_namespaces.h>

// test stuff
#include "../tests.h"

// Test of SAMRAIGhostDataAccumulator::accumulateGhostData() with a physical boundary operator. The data hold
// pseudorandom values in the ghost regions. The test compares the result with the sum, over all patches of a level and
// all periodic images, of the copies of each value after the transpose of the boundary fill. It also compares a call
// that is restricted to active patches with a call that is not, on data that are zero on the other patches.

namespace
{
using PatchArrays = std::vector<std::unique_ptr<ArrayData<NDIM, double>>>;

// A value that depends only on the level, the patch box, the component, the index and the depth.
double
random_value(const int ln, const Box<NDIM>& patch_box, const int component, const hier::Index<NDIM>& i, const int depth)
{
    std::uint64_t h = 0x9E3779B97F4A7C15ull;
    auto mix = [&h](const std::int64_t v)
    {
        h ^= static_cast<std::uint64_t>(v + (1 << 20)) + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
        h *= 0xBF58476D1CE4E5B9ull;
        h ^= h >> 31;
    };
    mix(ln);
    for (int d = 0; d < NDIM; ++d)
    {
        mix(patch_box.lower(d));
        mix(patch_box.upper(d));
    }
    mix(component);
    for (int d = 0; d < NDIM; ++d)
    {
        mix(i(d));
    }
    mix(depth);
    return 2.0 * static_cast<double>(h >> 11) / 9007199254740992.0 - 1.0;
}

Box<NDIM>
data_box(const bool cell, const Box<NDIM>& cell_box, const int component)
{
    return cell ? cell_box : SideGeometry<NDIM>::toSideBox(cell_box, component);
}

int
num_components(const bool cell)
{
    return cell ? 1 : NDIM;
}

ArrayData<NDIM, double>&
array_data(const bool cell, PatchData<NDIM>& data, const int component)
{
    if (cell)
    {
        return dynamic_cast<CellData<NDIM, double>&>(data).getArrayData();
    }
    return dynamic_cast<SideData<NDIM, double>&>(data).getArrayData(component);
}

// All shifts that map a location to its periodic images.
std::vector<IntVector<NDIM>>
periodic_shifts(Pointer<PatchLevel<NDIM>> level)
{
    const IntVector<NDIM> period = level->getGridGeometry()->getPeriodicShift(level->getRatio());
    std::vector<IntVector<NDIM>> shifts(1, IntVector<NDIM>(0));
    for (int d = 0; d < NDIM; ++d)
    {
        if (period(d) == 0)
        {
            continue;
        }
        const std::size_t n = shifts.size();
        for (std::size_t k = 0; k < n; ++k)
        {
            for (const int sign : { -1, 1 })
            {
                IntVector<NDIM> shift = shifts[k];
                shift(d) += sign * period(d);
                shifts.push_back(shift);
            }
        }
    }
    return shifts;
}

// Copies of the data of every patch of the level on every process.
std::vector<PatchArrays>
gather_level(Pointer<PatchLevel<NDIM>> level,
             const int idx,
             const bool cell,
             const int depth,
             const IntVector<NDIM>& gcw)
{
    const BoxArray<NDIM>& boxes = level->getBoxes();
    std::vector<PatchArrays> result(boxes.getNumberOfBoxes());
    for (int p = 0; p < boxes.getNumberOfBoxes(); ++p)
    {
        const Box<NDIM> ghost_box = Box<NDIM>::grow(boxes[p], gcw);
        for (int component = 0; component < num_components(cell); ++component)
        {
            result[p].emplace_back(new ArrayData<NDIM, double>(data_box(cell, ghost_box, component), depth));
            ArrayData<NDIM, double>& array = *result[p].back();
            array.fillAll(0.0);
            if (level->getMappingForPatch(p) == IBTK_MPI::getRank())
            {
                Pointer<PatchData<NDIM>> patch_data = level->getPatch(p)->getPatchData(idx);
                array.copy(array_data(cell, *patch_data, component), array.getBox());
            }
            IBTK_MPI::sumReduction(array.getPointer(), array.getBox().size() * depth);
        }
    }
    return result;
}

void
fill_levels(Pointer<PatchHierarchy<NDIM>> hierarchy, const int idx, const bool cell, const int depth)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<PatchData<NDIM>> patch_data = patch->getPatchData(idx);
            for (int component = 0; component < num_components(cell); ++component)
            {
                ArrayData<NDIM, double>& array = array_data(cell, *patch_data, component);
                for (int d = 0; d < depth; ++d)
                {
                    for (Box<NDIM>::Iterator i(array.getBox()); i; i++)
                    {
                        array(i(), d) = random_value(ln, patch->getBox(), component, i(), d);
                    }
                }
            }
        }
    }
}

// Set the values at sides that are in the ghost region of a patch p and on the boundary of another patch q (or of a
// periodic image of q), where p and q share no cell, to the given value. Returns the number of such sides on this
// process.
int
set_undeliverable_sides(Pointer<PatchLevel<NDIM>> level, const int idx, const int depth, const double value)
{
    const BoxArray<NDIM>& boxes = level->getBoxes();
    const std::vector<IntVector<NDIM>> shifts = periodic_shifts(level);
    int count = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, double>> side_data = patch->getPatchData(idx);
        const Box<NDIM> ghost_box = side_data->getGhostBox();
        for (int axis = 0; axis < NDIM; ++axis)
        {
            ArrayData<NDIM, double>& array = side_data->getArrayData(axis);
            for (Box<NDIM>::Iterator i(array.getBox()); i; i++)
            {
                bool undeliverable = false;
                for (int q = 0; q < boxes.getNumberOfBoxes() && !undeliverable; ++q)
                {
                    for (const IntVector<NDIM>& shift : shifts)
                    {
                        if (q == p() && shift == IntVector<NDIM>(0))
                        {
                            continue;
                        }
                        const Box<NDIM> q_box = Box<NDIM>::shift(boxes[q], shift);
                        const Box<NDIM> q_side_box = SideGeometry<NDIM>::toSideBox(q_box, axis);
                        if (q_side_box.contains(i()) && (ghost_box * q_box).empty())
                        {
                            undeliverable = true;
                            break;
                        }
                    }
                }
                if (undeliverable)
                {
                    ++count;
                    for (int d = 0; d < depth; ++d)
                    {
                        array(i(), d) = value;
                    }
                }
            }
        }
    }
    return count;
}

void
copy_levels(Pointer<PatchHierarchy<NDIM>> hierarchy, const int dst_idx, const int src_idx)
{
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            patch->getPatchData(dst_idx)->copy(*patch->getPatchData(src_idx));
        }
    }
}

bool
is_active(const Pointer<Patch<NDIM>> patch)
{
    return patch->getPatchNumber() % 2 == 0;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    // prevent a warning about timer initializations
    TimerManager::createManager(nullptr);
    {
        Pointer<Logger::Appender> abort_append(new TestAppender());
        Logger::getInstance()->setAbortAppender(abort_append);

        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "ghost_accumulation_02.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

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

        const bool cell = input_db->getString("var_type") == "CELL";
        const int depth = input_db->getInteger("depth");
        const IntVector<NDIM> gcw(input_db->getInteger("ghost_width"));
        const double gap_value = input_db->getDouble("undeliverable_side_value");

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("context");
        Pointer<hier::Variable<NDIM>> u_var =
            cell ? Pointer<hier::Variable<NDIM>>(new CellVariable<NDIM, double>("u", depth)) :
                   Pointer<hier::Variable<NDIM>>(new SideVariable<NDIM, double>("u", depth));
        const int original_idx = var_db->registerVariableAndContext(u_var, ctx, gcw);
        const int transposed_idx = var_db->registerClonedPatchDataIndex(u_var, original_idx);
        const int result_idx = var_db->registerClonedPatchDataIndex(u_var, original_idx);
        const int unpruned_idx = var_db->registerClonedPatchDataIndex(u_var, original_idx);
        const int pruned_idx = var_db->registerClonedPatchDataIndex(u_var, original_idx);

        // set up grid
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);
        const int tag_buffer = std::numeric_limits<int>::max();
        int level_number = 0;
        while (gridding_algorithm->levelCanBeRefined(level_number))
        {
            gridding_algorithm->makeFinerLevel(patch_hierarchy, 0.0, 0.0, tag_buffer);
            ++level_number;
        }
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (const int idx : { original_idx, transposed_idx, result_idx, unpruned_idx, pruned_idx })
            {
                level->allocatePatchData(idx, 0.0);
            }
        }

        // Boundary operator with Dirichlet and Neumann conditions on different parts of the boundary. The data have the
        // given depth, and the conditions are the same for all depths.
        std::vector<std::unique_ptr<muParserRobinBcCoefs>> bc_coefs;
        std::vector<RobinBcCoefStrategy<NDIM>*> bc_coef_ptrs;
        for (int d = 0; d < num_components(cell) * depth; ++d)
        {
            // The coefficients of all depths are those of the first depth.
            bc_coefs.emplace_back(
                new muParserRobinBcCoefs("bc_coefs_" + std::to_string(d),
                                         input_db->getDatabase("BcCoefs_" + std::to_string(d % num_components(cell))),
                                         grid_geometry));
            bc_coef_ptrs.push_back(bc_coefs.back().get());
        }
        std::unique_ptr<RobinPhysBdryPatchStrategy> bdry_op;
        if (cell)
        {
            bdry_op = std::make_unique<CartCellRobinPhysBdryOp>(result_idx, bc_coef_ptrs, /*homogeneous_bc*/ false);
        }
        else
        {
            bdry_op = std::make_unique<CartSideRobinPhysBdryOp>(result_idx, bc_coef_ptrs, /*homogeneous_bc*/ false);
        }

        // Random values in the interiors and ghost regions of all levels, and a chosen value at sides that the
        // accumulation cannot deliver.
        fill_levels(patch_hierarchy, original_idx, cell, depth);
        int num_undeliverable = 0;
        if (!cell)
        {
            for (int ln = 0; ln <= finest_ln; ++ln)
            {
                num_undeliverable +=
                    set_undeliverable_sides(patch_hierarchy->getPatchLevel(ln), original_idx, depth, gap_value);
            }
        }
        num_undeliverable = IBTK_MPI::sumReduction(num_undeliverable);
        for (const int idx : { transposed_idx, result_idx })
        {
            copy_levels(patch_hierarchy, idx, original_idx);
        }

        // Accumulate with the new call.
        SAMRAIGhostDataAccumulator accumulator(patch_hierarchy, u_var, gcw, 0, finest_ln);
        accumulator.accumulateGhostData(result_idx, bdry_op.get(), 0.0);

        // The reference: apply the transpose of the boundary fill to the patches of the copy and add up the copies of
        // each value over all patches and periodic images.
        bdry_op->setPatchDataIndex(transposed_idx);
        int num_values = 0;
        double sum_values = 0.0;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                bdry_op->accumulateFromPhysicalBoundaryData(*level->getPatch(p()), 0.0, gcw);
            }
            const std::vector<PatchArrays> copies = gather_level(level, transposed_idx, cell, depth, gcw);
            const BoxArray<NDIM>& boxes = level->getBoxes();
            const std::vector<IntVector<NDIM>> shifts = periodic_shifts(level);
            double max_error = 0.0, max_value = 0.0;
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<PatchData<NDIM>> result = patch->getPatchData(result_idx);
                for (int component = 0; component < num_components(cell); ++component)
                {
                    ArrayData<NDIM, double>& array = array_data(cell, *result, component);
                    for (Box<NDIM>::Iterator i(data_box(cell, patch->getBox(), component)); i; i++)
                    {
                        for (int d = 0; d < depth; ++d)
                        {
                            double total = 0.0;
                            for (int q = 0; q < boxes.getNumberOfBoxes(); ++q)
                            {
                                for (const IntVector<NDIM>& shift : shifts)
                                {
                                    const hier::Index<NDIM> j = i() + shift;
                                    const ArrayData<NDIM, double>& copy = *copies[q][component];
                                    if (copy.getBox().contains(j))
                                    {
                                        total += copy(j, d);
                                    }
                                }
                            }
                            max_error = std::max(max_error, std::abs(array(i(), d) - total));
                            max_value = std::max(max_value, std::abs(total));
                            sum_values += array(i(), d);
                            ++num_values;
                        }
                    }
                }
            }
            max_error = IBTK_MPI::maxReduction(max_error);
            max_value = IBTK_MPI::maxReduction(max_value);
            if (!(max_error <= 1.0e-12 * std::max(1.0, max_value)))
            {
                TBOX_ERROR("accumulated values on level " << ln << " differ from the sum of the copies by " << max_error
                                                          << "\n");
            }
        }
        num_values = IBTK_MPI::sumReduction(num_values);
        sum_values = IBTK_MPI::sumReduction(sum_values);

        // Accumulate with and without the restriction to active patches, on data that are zero on the other patches.
        std::vector<std::vector<int>> active_patch_nums(finest_ln + 1);
        int num_active = 0;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                for (const int idx : { unpruned_idx, pruned_idx })
                {
                    Pointer<PatchData<NDIM>> patch_data = patch->getPatchData(idx);
                    patch_data->copy(*patch->getPatchData(original_idx));
                    if (is_active(patch))
                    {
                        continue;
                    }
                    for (int component = 0; component < num_components(cell); ++component)
                    {
                        array_data(cell, *patch_data, component).fillAll(0.0);
                    }
                }
                if (is_active(patch))
                {
                    active_patch_nums[ln].push_back(patch->getPatchNumber());
                    ++num_active;
                }
            }
        }
        num_active = IBTK_MPI::sumReduction(num_active);
        accumulator.accumulateGhostData(unpruned_idx, bdry_op.get(), 0.0);
        accumulator.accumulateGhostData(pruned_idx, bdry_op.get(), 0.0, &active_patch_nums);
        double max_difference = 0.0, max_value = 0.0;
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = level->getPatch(p());
                Pointer<PatchData<NDIM>> unpruned = patch->getPatchData(unpruned_idx);
                Pointer<PatchData<NDIM>> pruned = patch->getPatchData(pruned_idx);
                for (int component = 0; component < num_components(cell); ++component)
                {
                    const ArrayData<NDIM, double>& a = array_data(cell, *unpruned, component);
                    const ArrayData<NDIM, double>& b = array_data(cell, *pruned, component);
                    for (Box<NDIM>::Iterator i(data_box(cell, patch->getBox(), component)); i; i++)
                    {
                        for (int d = 0; d < depth; ++d)
                        {
                            max_difference = std::max(max_difference, std::abs(a(i(), d) - b(i(), d)));
                            max_value = std::max(max_value, std::abs(a(i(), d)));
                        }
                    }
                }
            }
        }
        max_difference = IBTK_MPI::maxReduction(max_difference);
        max_value = IBTK_MPI::maxReduction(max_value);
        if (!(max_difference <= 1.0e-14 * std::max(1.0, max_value)))
        {
            TBOX_ERROR("accumulation restricted to active patches differs by " << max_difference << "\n");
        }

        plog << "number of levels: " << finest_ln + 1 << '\n';
        for (int ln = 0; ln <= finest_ln; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(ln);
            plog << "number of patches on level " << ln << ": " << level->getNumberOfPatches() << '\n';
        }
        plog << "number of undeliverable sides set: " << num_undeliverable << '\n';
        plog << "number of values compared: " << num_values << '\n';
        plog << "sum of accumulated values: " << std::setprecision(10) << sum_values << '\n';
        plog << "number of active patches: " << num_active << '\n';
    }
}
