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
#include <ibtk/CartGridPatchwiseFunction.h>
#include <ibtk/CartesianCentering.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/Logger.h>

#include <BergerRigoutsos.h>
#include <CartesianPatchGeometry.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <utility>

#include <ibtk/app_namespaces.h>

namespace
{
using Cell = CartesianCentering<DataCentering::CELL>;
using Side = CartesianCentering<DataCentering::SIDE>;

VectorNd
position(const Patch<NDIM>& patch, const SAMRAI::hier::Index<NDIM>& index, const VectorNd& offset)
{
    const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch.getPatchGeometry();
    VectorNd x;
    for (int d = 0; d < NDIM; ++d)
    {
        x[d] = geometry->getXLower()[d] + geometry->getDx()[d] * (index(d) - patch.getBox().lower()(d) + offset[d]);
    }
    return x;
}

std::pair<double, double>
source_values(const VectorNd& x)
{
    double h = 0.25, liquid = 0.5;
    for (int d = 0; d < NDIM; ++d)
    {
        h += (d + 1) * x[d] * x[d] / 1024.0;
        liquid += (d + 2) * x[d] * x[d] / 2048.0;
    }
    return { h, liquid };
}

template <Centering Layout>
double
result_error(PatchHierarchy<NDIM>& hierarchy, const int data_idx, const double time, const bool initial_time)
{
    double error = 0.0;
    for (int ln = 0; ln <= hierarchy.getFinestLevelNumber(); ++ln)
    {
        const Pointer<PatchLevel<NDIM>> level = hierarchy.getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Pointer<typename Layout::template Data<double>> data = patch->getPatchData(data_idx);
            const Pointer<CartesianPatchGeometry<NDIM>> geometry = patch->getPatchGeometry();
            const double* const dx = geometry->getDx();
            for (int axis = 0; axis < (Layout::is_oriented() ? NDIM : 1); ++axis)
            {
                const ArrayData<NDIM, double>& array = [&]() -> const ArrayData<NDIM, double>&
                {
                    if constexpr (Layout::is_oriented())
                    {
                        return data->getArrayData(axis);
                    }
                    else
                    {
                        return data->getArrayData();
                    }
                }();
                Box<NDIM> interior = patch->getBox();
                double scale = 1.0;
                if constexpr (Layout::is_oriented())
                {
                    ++interior.upper()(axis);
                    scale = dx[axis];
                }
                else
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        scale *= dx[d];
                    }
                }
                for (Box<NDIM>::Iterator it(array.getBox()); it; it++)
                {
                    const double actual = array(it(), 0);
                    if (!interior.contains(it()))
                    {
                        TBOX_ASSERT(std::isnan(actual));
                        continue;
                    }
                    const VectorNd x = position(*patch, it(), Layout::offset(axis));
                    std::pair<double, double> values = source_values(x);
                    if constexpr (Layout::is_oriented())
                    {
                        // The average of a quadratic at x +/- dx/2 includes this curvature term.
                        values.first += (axis + 1) * dx[axis] * dx[axis] / 4096.0;
                        values.second += (axis + 2) * dx[axis] * dx[axis] / 8192.0;
                    }
                    const double expected = 1.0 + 2.0 * values.first + 4.0 * values.first * values.second +
                                            time * scale + (initial_time ? 1.0 : 0.0);
                    TBOX_ASSERT(std::isfinite(actual));
                    error = std::max(error, std::abs(actual - expected));
                }
            }
        }
    }
    return IBTK_MPI::maxReduction(error);
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    const Pointer<AppInitializer> app = new AppInitializer(argc, argv);
    Logger::getInstance()->setWarning(false);
    const Pointer<CartesianGridGeometry<NDIM>> geometry =
        new CartesianGridGeometry<NDIM>("geometry", app->getComponentDatabase("CartesianGeometry"));
    const Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("hierarchy", geometry);
    const Pointer<StandardTagAndInitialize<NDIM>> error_detector =
        new StandardTagAndInitialize<NDIM>("tagging", nullptr, app->getComponentDatabase("StandardTagAndInitialize"));
    const Pointer<BergerRigoutsos<NDIM>> boxes = new BergerRigoutsos<NDIM>();
    const Pointer<LoadBalancer<NDIM>> load = new LoadBalancer<NDIM>("load", app->getComponentDatabase("LoadBalancer"));
    const Pointer<GriddingAlgorithm<NDIM>> gridding = new GriddingAlgorithm<NDIM>(
        "gridding", app->getComponentDatabase("GriddingAlgorithm"), error_detector, boxes, load);

    VariableDatabase<NDIM>* const var_db = VariableDatabase<NDIM>::getDatabase();
    const Pointer<VariableContext> context = var_db->getContext("data");
    const Pointer<CellVariable<NDIM, double>> h_var = new CellVariable<NDIM, double>("H");
    const Pointer<CellVariable<NDIM, double>> liquid_var = new CellVariable<NDIM, double>("liquid fraction");
    const Pointer<CellVariable<NDIM, double>> cell_var = new CellVariable<NDIM, double>("cell blend");
    const Pointer<SideVariable<NDIM, double>> side_var = new SideVariable<NDIM, double>("side blend");
    IntVector<NDIM> output_ghosts(1);
    output_ghosts(1) = 2;
    const int h_idx = var_db->registerVariableAndContext(h_var, context, IntVector<NDIM>(1));
    const int liquid_idx = var_db->registerVariableAndContext(liquid_var, context, IntVector<NDIM>(1));
    const int cell_idx = var_db->registerVariableAndContext(cell_var, context, output_ghosts);
    const int side_idx = var_db->registerVariableAndContext(side_var, context, output_ghosts);
    gridding->makeCoarsestLevel(hierarchy, 0.0);
    gridding->makeFinerLevel(hierarchy, 0.0, 0.0, 1);
    TBOX_ASSERT(hierarchy->getFinestLevelNumber() == 1);

    int local_patches = 0;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        const Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (const int idx : { h_idx, liquid_idx, cell_idx, side_idx })
        {
            level->allocatePatchData(idx);
        }
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            ++local_patches;
            const Pointer<Patch<NDIM>> patch = level->getPatch(p());
            const Pointer<Cell::Data<double>> h = patch->getPatchData(h_idx);
            const Pointer<Cell::Data<double>> liquid = patch->getPatchData(liquid_idx);
            const Pointer<Cell::Data<double>> cell = patch->getPatchData(cell_idx);
            const Pointer<Side::Data<double>> side = patch->getPatchData(side_idx);
            cell->fillAll(std::numeric_limits<double>::quiet_NaN());
            side->fillAll(std::numeric_limits<double>::quiet_NaN());
            // Source ghosts are supplied by the caller, including physical and coarse-fine boundaries.
            for (auto it = Cell::begin(h->getGhostBox(), 0); it; it++)
            {
                const std::pair<double, double> values = source_values(position(*patch, it(), Cell::offset(0)));
                (*h)(it()) = values.first;
                (*liquid)(it()) = values.second;
            }
        }
    }
    TBOX_ASSERT(IBTK_MPI::sumReduction(local_patches) > 2);

    double expected_time = 0.25;
    bool expected_initial = true, expected_level = true;
    std::map<const Patch<NDIM>*, int> cell_calls, side_calls;
    const auto check_context =
        [&](Pointer<Patch<NDIM>> patch, const double time, const bool initial_time, Pointer<PatchLevel<NDIM>> level)
    {
        TBOX_ASSERT(time == expected_time && initial_time == expected_initial);
        TBOX_ASSERT(static_cast<bool>(level) == expected_level);
        if (level)
        {
            TBOX_ASSERT(level == hierarchy->getPatchLevel(patch->getPatchLevelNumber()));
            TBOX_ASSERT(level->getPatch(patch->getPatchNumber()) == patch);
        }
    };
    auto cell_callback =
        [&, densities = std::make_unique<std::array<double, 3>>(std::array<double, 3>{ 1.0, 3.0, 7.0 })](
            const int data_idx,
            Pointer<Variable<NDIM>> var,
            Pointer<Patch<NDIM>> patch,
            const double time,
            const bool initial_time,
            Pointer<PatchLevel<NDIM>> level)
    {
        TBOX_ASSERT(data_idx == cell_idx && var == cell_var);
        check_context(patch, time, initial_time, level);
        ++cell_calls[patch.getPointer()];
        const Pointer<Cell::Data<double>> h_data = patch->getPatchData(h_idx);
        const Pointer<Cell::Data<double>> liquid_data = patch->getPatchData(liquid_idx);
        const Cell::Data<double>& h = *h_data;
        const Cell::Data<double>& liquid = *liquid_data;
        const Pointer<Cell::Data<double>> dst = patch->getPatchData(data_idx);
        const Pointer<CartesianPatchGeometry<NDIM>> patch_geometry = patch->getPatchGeometry();
        const double* const dx = patch_geometry->getDx();
        double volume = 1.0;
        for (int d = 0; d < NDIM; ++d)
        {
            volume *= dx[d];
        }
        for (auto it = Cell::begin(patch->getBox(), 0); it; it++)
        {
            const double material = (*densities)[1] + ((*densities)[2] - (*densities)[1]) * liquid(it());
            (*dst)(it()) =
                (*densities)[0] * (1.0 - h(it())) + material * h(it()) + time * volume + (initial_time ? 1.0 : 0.0);
        }
    };
    static_assert(PatchwiseCallback<decltype(cell_callback)>);
    static_assert(!std::copy_constructible<decltype(cell_callback)>);
    const Pointer<CartGridFunction> cell_function =
        make_cart_grid_patchwise_function("cell blend", std::move(cell_callback));

    Pointer<CartGridFunction> side_function;
    int observed_side_calls = 0;
    {
        // The factory must copy this const lvalue and invoke its mutable copy after this scope ends.
        const auto side_callback = [&, calls = 0](const int data_idx,
                                                  Pointer<Variable<NDIM>> var,
                                                  Pointer<Patch<NDIM>> patch,
                                                  const double time,
                                                  const bool initial_time,
                                                  Pointer<PatchLevel<NDIM>> level) mutable
        {
            TBOX_ASSERT(data_idx == side_idx && var == side_var);
            check_context(patch, time, initial_time, level);
            ++side_calls[patch.getPointer()];
            observed_side_calls = ++calls;
            const Pointer<Cell::Data<double>> h_data = patch->getPatchData(h_idx);
            const Pointer<Cell::Data<double>> liquid_data = patch->getPatchData(liquid_idx);
            const Cell::Data<double>& h = *h_data;
            const Cell::Data<double>& liquid = *liquid_data;
            const Pointer<Side::Data<double>> dst = patch->getPatchData(data_idx);
            const Pointer<CartesianPatchGeometry<NDIM>> patch_geometry = patch->getPatchGeometry();
            const double* const dx = patch_geometry->getDx();
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (auto it = Side::begin(patch->getBox(), axis); it; it++)
                {
                    const Side::Index& side = it();
                    const double h_avg = 0.5 * (h(side.toCell(0)) + h(side.toCell(1)));
                    const double liquid_avg = 0.5 * (liquid(side.toCell(0)) + liquid(side.toCell(1)));
                    (*dst)(side) =
                        (1.0 - h_avg) + (3.0 + 4.0 * liquid_avg) * h_avg + time * dx[axis] + (initial_time ? 1.0 : 0.0);
                }
            }
        };
        static_assert(!PatchwiseCallback<decltype(side_callback)>);
        static_assert(PatchwiseCallback<std::decay_t<decltype(side_callback)>>);
        side_function = make_cart_grid_patchwise_function("side blend", side_callback);
    }
    TBOX_ASSERT(cell_function->isTimeDependent() && side_function->isTimeDependent());
    cell_function->setDataOnPatchHierarchy(cell_idx, cell_var, hierarchy, expected_time, expected_initial);
    const double cell_error = result_error<Cell>(*hierarchy, cell_idx, expected_time, expected_initial);
    expected_time = 0.5;
    expected_initial = false;
    side_function->setDataOnPatchHierarchy(side_idx, side_var, hierarchy, expected_time, expected_initial);
    const double side_error = result_error<Side>(*hierarchy, side_idx, expected_time, expected_initial);
    TBOX_ASSERT(observed_side_calls == local_patches);

    expected_time = 1.25;
    expected_level = false;
    for (int ln = 0; ln <= hierarchy->getFinestLevelNumber(); ++ln)
    {
        const Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            const Pointer<Patch<NDIM>> patch = level->getPatch(p());
            TBOX_ASSERT(cell_calls[patch.getPointer()] == 1 && side_calls[patch.getPointer()] == 1);
            side_function->setDataOnPatch(side_idx, side_var, patch, expected_time);
            TBOX_ASSERT(side_calls[patch.getPointer()] == 2);
            const Pointer<Cell::Data<double>> h_data = patch->getPatchData(h_idx);
            const Pointer<Cell::Data<double>> liquid_data = patch->getPatchData(liquid_idx);
            const Cell::Data<double>& h = *h_data;
            const Cell::Data<double>& liquid = *liquid_data;
            for (auto it = Cell::begin(h.getGhostBox(), 0); it; it++)
            {
                const std::pair<double, double> expected = source_values(position(*patch, it(), Cell::offset(0)));
                TBOX_ASSERT(h(it()) == expected.first && liquid(it()) == expected.second);
            }
        }
    }
    TBOX_ASSERT(observed_side_calls == 2 * local_patches);
    TBOX_ASSERT(cell_calls.size() == static_cast<std::size_t>(local_patches));
    TBOX_ASSERT(side_calls.size() == static_cast<std::size_t>(local_patches));
    const double direct_error = result_error<Side>(*hierarchy, side_idx, expected_time, expected_initial);
    plog << std::scientific << std::setprecision(12);
    plog << "cell blend error: " << cell_error << '\n';
    plog << "side blend errors (hierarchy, patch): " << side_error << ' ' << direct_error << '\n';
    return 0;
}
