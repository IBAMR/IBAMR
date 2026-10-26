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

#include <ArrayData.h>
#include <Box.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CellIterator.h>
#include <EdgeData.h>
#include <EdgeGeometry.h>
#include <EdgeVariable.h>
#include <FaceData.h>
#include <FaceGeometry.h>
#include <FaceVariable.h>
#include <Index.h>
#include <IntVector.h>
#include <NodeData.h>
#include <NodeGeometry.h>
#include <NodeVariable.h>
#include <OuteredgeData.h>
#include <OuteredgeVariable.h>
#include <OuterfaceData.h>
#include <OuterfaceVariable.h>
#include <OuternodeData.h>
#include <OuternodeVariable.h>
#include <OutersideData.h>
#include <OutersideVariable.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <string>
#include <vector>

#include <ibtk/app_namespaces.h>

// Report how many ghost values a refine schedule sets for items whose sources
// are outer node, side, face or edge data, each paired with an item of the same
// centering whose source is its own data. Each source value depends on its
// index and axis, and a ghost value inside the domain is set if it holds the
// value for its own index and axis.

namespace
{
const double NOT_SET = -1.0e3;

// On level 0, the lower half in the first direction has 2^(NDIM-1) patches and
// the upper half is one patch, so that sources include both face and diagonal
// neighbors in 3D.
tbox::Pointer<hier::PatchLevel<NDIM>>
make_level(tbox::Pointer<geom::CartesianGridGeometry<NDIM>> grid_geometry, const std::vector<int>& idxs)
{
    tbox::Pointer<hier::PatchHierarchy<NDIM>> hierarchy =
        new hier::PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    const hier::Box<NDIM> domain = grid_geometry->getPhysicalDomain()[0];
    const int n = domain.numberCells(0) / 2;
    const int n_lower = 1 << (NDIM - 1);
    hier::BoxArray<NDIM> boxes(n_lower + 1);
    for (int k = 0; k < n_lower; ++k)
    {
        hier::Index<NDIM> lower(0);
        for (int d = 1; d < NDIM; ++d)
        {
            lower(d) = ((k >> (d - 1)) & 1) * n;
        }
        boxes[k] = hier::Box<NDIM>(lower, lower + (n - 1));
    }
    hier::Index<NDIM> upper_lower(0);
    upper_lower(0) = n;
    boxes[n_lower] = hier::Box<NDIM>(upper_lower, domain.upper());
    hier::ProcessorMapping mapping(boxes.size());
    for (int k = 0; k < boxes.size(); ++k)
    {
        mapping.setProcessorAssignment(k, k % IBTK::IBTK_MPI::getNodes());
    }
    hierarchy->makeNewPatchLevel(0, hier::IntVector<NDIM>(1), boxes, mapping);
    tbox::Pointer<hier::PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    for (const int idx : idxs)
    {
        level->allocatePatchData(idx, 0.0);
    }
    return level;
}

pdat::ArrayData<NDIM, double>&
array(pdat::NodeData<NDIM, double>& data, int)
{
    return data.getArrayData();
}

template <class DataType>
pdat::ArrayData<NDIM, double>&
array(DataType& data, const int axis)
{
    return data.getArrayData(axis);
}

// The source value at index i of the given axis.
double
source_value(const hier::Index<NDIM>& i, const int axis)
{
    double value = axis + 1;
    for (int d = 0; d < NDIM; ++d)
    {
        value = 64.0 * value + i(d);
    }
    return value;
}

// Set the values on the patch to the source values and the ghost values to
// NOT_SET.
template <class DataType>
void
set_values(tbox::Pointer<hier::PatchLevel<NDIM>> level,
           const int idx,
           const int n_arrays,
           hier::Box<NDIM> (*index_box)(const hier::Box<NDIM>&, int))
{
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        tbox::Pointer<DataType> data = patch->getPatchData(idx);
        data->fillAll(NOT_SET);
        for (int n = 0; n < n_arrays; ++n)
        {
            for (pdat::CellIterator<NDIM> i(index_box(patch->getBox(), n)); i; i++)
            {
                array(*data, n)(i(), 0) = source_value(i(), n);
            }
        }
    }
}

// Count the ghost values inside the domain that hold their source values.
template <class DataType>
void
report_values(const std::string& label,
              tbox::Pointer<hier::PatchLevel<NDIM>> level,
              const int idx,
              const int n_arrays,
              hier::Box<NDIM> (*index_box)(const hier::Box<NDIM>&, int))
{
    const hier::Box<NDIM> domain = level->getPhysicalDomain()[0];
    int values = 0, values_set = 0;
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        tbox::Pointer<DataType> data = patch->getPatchData(idx);
        for (int n = 0; n < n_arrays; ++n)
        {
            pdat::ArrayData<NDIM, double>& values_array = array(*data, n);
            const hier::Box<NDIM> interior = index_box(patch->getBox(), n);
            const hier::Box<NDIM> domain_indices = index_box(domain, n);
            for (pdat::CellIterator<NDIM> i(values_array.getBox()); i; i++)
            {
                if (interior.contains(i()) || !domain_indices.contains(i()))
                {
                    continue;
                }
                ++values;
                values_set += values_array(i(), 0) == source_value(i(), n);
            }
        }
    }
    values = IBTK::IBTK_MPI::sumReduction(values);
    values_set = IBTK::IBTK_MPI::sumReduction(values_set);
    tbox::plog << label << ": ghost values set " << values_set << " of " << values << '\n';
}

hier::Box<NDIM>
node_box(const hier::Box<NDIM>& box, int)
{
    return pdat::NodeGeometry<NDIM>::toNodeBox(box);
}

template <template <int, class> class VariableType,
          template <int, class>
          class OuterVariableType,
          template <int, class>
          class DataType,
          template <int, class>
          class OuterDataType>
void
test_outer_source(const std::string& name,
                  tbox::Pointer<geom::CartesianGridGeometry<NDIM>> grid_geometry,
                  const int n_arrays,
                  hier::Box<NDIM> (*index_box)(const hier::Box<NDIM>&, int))
{
    auto* var_db = hier::VariableDatabase<NDIM>::getDatabase();
    tbox::Pointer<VariableType<NDIM, double>> var = new VariableType<NDIM, double>(name);
    tbox::Pointer<OuterVariableType<NDIM, double>> outer_var = new OuterVariableType<NDIM, double>(name + "::outer");
    const int from_outer_idx =
        var_db->registerVariableAndContext(var, var_db->getContext("from outer"), hier::IntVector<NDIM>(2));
    const int plain_idx =
        var_db->registerVariableAndContext(var, var_db->getContext("plain"), hier::IntVector<NDIM>(2));
    const int outer_idx =
        var_db->registerVariableAndContext(outer_var, var_db->getContext("outer"), hier::IntVector<NDIM>(0));
    tbox::Pointer<hier::PatchLevel<NDIM>> level = make_level(grid_geometry, { from_outer_idx, plain_idx, outer_idx });
    xfer::RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(plain_idx, plain_idx, plain_idx, nullptr);
    refine_alg.registerRefine(from_outer_idx, outer_idx, from_outer_idx, nullptr);
    set_values<DataType<NDIM, double>>(level, plain_idx, n_arrays, index_box);
    set_values<DataType<NDIM, double>>(level, from_outer_idx, n_arrays, index_box);
    for (hier::PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        tbox::Pointer<hier::Patch<NDIM>> patch = level->getPatch(p());
        tbox::Pointer<OuterDataType<NDIM, double>> outer_data = patch->getPatchData(outer_idx);
        outer_data->copy(*patch->getPatchData(plain_idx));
    }
    refine_alg.createSchedule(level)->fillData(0.0);
    report_values<DataType<NDIM, double>>(name + ", plain item", level, plain_idx, n_arrays, index_box);
    report_values<DataType<NDIM, double>>(name + ", item from outer data", level, from_outer_idx, n_arrays, index_box);
}
} // namespace

int
main(int argc, char** argv)
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    tbox::Pointer<IBTK::AppInitializer> app_initializer = new IBTK::AppInitializer(argc, argv);
    tbox::Pointer<geom::CartesianGridGeometry<NDIM>> grid_geometry = new geom::CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    test_outer_source<pdat::NodeVariable, pdat::OuternodeVariable, pdat::NodeData, pdat::OuternodeData>(
        "node", grid_geometry, 1, node_box);
    test_outer_source<pdat::SideVariable, pdat::OutersideVariable, pdat::SideData, pdat::OutersideData>(
        "side", grid_geometry, NDIM, pdat::SideGeometry<NDIM>::toSideBox);
    test_outer_source<pdat::FaceVariable, pdat::OuterfaceVariable, pdat::FaceData, pdat::OuterfaceData>(
        "face", grid_geometry, NDIM, pdat::FaceGeometry<NDIM>::toFaceBox);
    test_outer_source<pdat::EdgeVariable, pdat::OuteredgeVariable, pdat::EdgeData, pdat::OuteredgeData>(
        "edge", grid_geometry, NDIM, pdat::EdgeGeometry<NDIM>::toEdgeBox);
}
