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

#include <ibtk/HierarchyGhostCellInterpolation.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/SideDataSynchronization.h>

#include <tbox/MemoryDatabase.h>

#include <Box.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CartesianPatchGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <RefineOperator.h>
#include <SideData.h>
#include <SideIndex.h>
#include <SideIterator.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <fstream>
#include <iomanip>
#include <string>

#include <ibtk/app_namespaces.h>

// Create grid geometries that exist at the same time, each with its own two-level hierarchy, and fill ghost cells and
// synchronize side data on each. Then destroy one geometry and do the same on a new one.
namespace
{
struct Mesh
{
    Pointer<CartesianGridGeometry<NDIM>> geometry;
    Pointer<PatchHierarchy<NDIM>> hierarchy;
};

struct Indices
{
    int cell = invalid_index;
    int side = invalid_index;
    int side_sync = invalid_index;
};

// The data are linear, so that the transfers reproduce them exactly.
double
exact_value(const VectorNd& x)
{
    double value = 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        value += static_cast<double>(d + 1) * x[d];
    }
    return value;
}

// Eight coarse cells in each direction and one fine patch that covers coarse cells 2 through 5.
Mesh
make_mesh(const std::string& name, const Indices& indices)
{
    Pointer<Database> db = new MemoryDatabase(name + "_db");
    const Box<NDIM> coarse_box(hier::Index<NDIM>(0), hier::Index<NDIM>(7));
    tbox::Array<tbox::DatabaseBox> domain_boxes(1);
    domain_boxes[0] = tbox::DatabaseBox(coarse_box);
    db->putDatabaseBoxArray("domain_boxes", domain_boxes);
    double x_lo[NDIM], x_up[NDIM];
    for (int d = 0; d < NDIM; ++d)
    {
        x_lo[d] = 0.0;
        x_up[d] = 1.0;
    }
    db->putDoubleArray("x_lo", x_lo, NDIM);
    db->putDoubleArray("x_up", x_up, NDIM);

    Mesh mesh;
    mesh.geometry = new CartesianGridGeometry<NDIM>(name, db);
    mesh.hierarchy = new PatchHierarchy<NDIM>(name + "_hierarchy", mesh.geometry);
    ProcessorMapping mapping(1);
    mapping.setProcessorAssignment(0, 0);
    const Box<NDIM> fine_box(hier::Index<NDIM>(4), hier::Index<NDIM>(11));
    mesh.hierarchy->makeNewPatchLevel(0, IntVector<NDIM>(1), BoxArray<NDIM>(coarse_box), mapping);
    mesh.hierarchy->makeNewPatchLevel(1, IntVector<NDIM>(2), BoxArray<NDIM>(fine_box), mapping);

    // The side data to synchronize is zero on the coarse level.
    for (int ln = 0; ln <= 1; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = mesh.hierarchy->getPatchLevel(ln);
        level->allocatePatchData(indices.cell);
        level->allocatePatchData(indices.side);
        level->allocatePatchData(indices.side_sync);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM>> patch_geom = patch->getPatchGeometry();
            const double* const dx = patch_geom->getDx();
            const double* const x_lower = patch_geom->getXLower();
            const hier::Index<NDIM>& index_lower = patch->getBox().lower();
            Pointer<CellData<NDIM, double>> cell_data = patch->getPatchData(indices.cell);
            Pointer<SideData<NDIM, double>> side_data = patch->getPatchData(indices.side);
            Pointer<SideData<NDIM, double>> side_sync_data = patch->getPatchData(indices.side_sync);
            cell_data->fillAll(0.0);
            side_data->fillAll(0.0);
            side_sync_data->fillAll(0.0);
            for (CellIterator<NDIM> ci(patch->getBox()); ci; ci++)
            {
                VectorNd x;
                for (int d = 0; d < NDIM; ++d)
                {
                    x[d] = x_lower[d] + dx[d] * (static_cast<double>(ci()(d) - index_lower(d)) + 0.5);
                }
                (*cell_data)(ci()) = exact_value(x);
            }
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (SideIterator<NDIM> si(patch->getBox(), axis); si; si++)
                {
                    VectorNd x;
                    for (int d = 0; d < NDIM; ++d)
                    {
                        x[d] = x_lower[d] +
                               dx[d] * (static_cast<double>(si()(d) - index_lower(d)) + (d == axis ? 0.0 : 0.5));
                    }
                    (*side_data)(si()) = exact_value(x);
                    if (ln == 1) (*side_sync_data)(si()) = exact_value(x);
                }
            }
        }
    }
    return mesh;
}

Pointer<Patch<NDIM>>
get_first_patch(const Mesh& mesh, const int ln)
{
    Pointer<PatchLevel<NDIM>> level = mesh.hierarchy->getPatchLevel(ln);
    return level->getPatch(PatchLevel<NDIM>::Iterator(level)());
}

// Fill ghost cells with cubic coarsening on a cell and a side component, synchronize side data with cubic coarsening,
// and write one ghost value and one synchronized value.
void
transfer(const std::string& name, const Mesh& mesh, const Indices& indices, std::ostream& out)
{
    using InterpolationTransactionComponent = HierarchyGhostCellInterpolation::InterpolationTransactionComponent;
    std::vector<InterpolationTransactionComponent> components;
    components.emplace_back(indices.cell, "NONE", true, "CUBIC_COARSEN", "LINEAR", false, nullptr);
    components.emplace_back(indices.side, "NONE", true, "CUBIC_COARSEN", "LINEAR", false, nullptr);
    HierarchyGhostCellInterpolation fill;
    fill.initializeOperatorState(components, mesh.hierarchy);
    fill.fillData(0.0);

    // The first ghost layer to the left of the fine patch, away from its corners.
    Pointer<Patch<NDIM>> fine_patch = get_first_patch(mesh, 1);
    hier::Index<NDIM> ghost_index(7);
    ghost_index(0) = fine_patch->getBox().lower(0) - 1;
    Pointer<CellData<NDIM, double>> cell_data = fine_patch->getPatchData(indices.cell);
    Pointer<SideData<NDIM, double>> side_data = fine_patch->getPatchData(indices.side);
    out << name << " cell ghost value = " << (*cell_data)(CellIndex<NDIM>(ghost_index)) << "\n";
    out << name << " side ghost value = " << (*side_data)(SideIndex<NDIM>(ghost_index, 0, SideIndex<NDIM>::Lower))
        << "\n";

    SideDataSynchronization synchronization;
    synchronization.initializeOperatorState(
        SideDataSynchronization::SynchronizationTransactionComponent(indices.side_sync, "CUBIC_COARSEN"),
        mesh.hierarchy);
    synchronization.synchronizeData(0.0);

    // A coarse side under the fine patch.
    Pointer<Patch<NDIM>> coarse_patch = get_first_patch(mesh, 0);
    Pointer<SideData<NDIM, double>> side_sync_data = coarse_patch->getPatchData(indices.side_sync);
    out << name << " synchronized side value = "
        << (*side_sync_data)(SideIndex<NDIM>(hier::Index<NDIM>(4), 0, SideIndex<NDIM>::Lower)) << "\n";
}

// Look up the RT0 side refine operator, which the transfer objects above register.
void
look_up_refine_operator(const std::string& name, const Mesh& mesh, const Indices& indices, std::ostream& out)
{
    Pointer<Variable<NDIM>> var;
    VariableDatabase<NDIM>::getDatabase()->mapIndexToVariable(indices.side, var);
    Pointer<RefineOperator<NDIM>> refine_op = mesh.geometry->lookupRefineOperator(var, "RT0_REFINE");
    out << name << " refine operator = " << refine_op->getOperatorName() << "\n";
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    {
        std::ofstream out("output");
        out << std::setprecision(10);

        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> context = var_db->getContext("CTX");
        Pointer<CellVariable<NDIM, double>> cell_var = new CellVariable<NDIM, double>("cell");
        Pointer<SideVariable<NDIM, double>> side_var = new SideVariable<NDIM, double>("side");
        Pointer<SideVariable<NDIM, double>> side_sync_var = new SideVariable<NDIM, double>("side_sync");
        Indices indices;
        indices.cell = var_db->registerVariableAndContext(cell_var, context, IntVector<NDIM>(1));
        indices.side = var_db->registerVariableAndContext(side_var, context, IntVector<NDIM>(1));
        indices.side_sync = var_db->registerVariableAndContext(side_sync_var, context, IntVector<NDIM>(1));

        // Two geometries at the same time.
        Mesh mesh_0 = make_mesh("geometry_0", indices);
        {
            Mesh mesh_1 = make_mesh("geometry_1", indices);
            transfer("geometry_0", mesh_0, indices, out);
            transfer("geometry_1", mesh_1, indices, out);
            look_up_refine_operator("geometry_0", mesh_0, indices, out);
            look_up_refine_operator("geometry_1", mesh_1, indices, out);
        }

        // A new geometry after the destruction of the second one.
        Mesh mesh_2 = make_mesh("geometry_2", indices);
        transfer("geometry_2", mesh_2, indices, out);
        look_up_refine_operator("geometry_2", mesh_2, indices, out);
    }
    return 0;
}
