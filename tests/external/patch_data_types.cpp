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

// Check that cell-centered data of types bool, char, and dcomplex are copied
// between processes when ghost cell values are filled.

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/Complex.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <GriddingAlgorithm.h>
#include <LoadBalancer.h>
#include <RefineAlgorithm.h>
#include <StandardTagAndInitialize.h>

#include <string>

#include <ibtk/app_namespaces.h>

namespace
{
template <typename T>
T cell_value(const hier::Index<NDIM>& idx);

template <>
bool
cell_value<bool>(const hier::Index<NDIM>& idx)
{
    return (idx(0) + 2 * idx(1)) % 3 == 0;
}

template <>
char
cell_value<char>(const hier::Index<NDIM>& idx)
{
    return static_cast<char>('a' + idx(0) + 6 * idx(1));
}

template <>
dcomplex
cell_value<dcomplex>(const hier::Index<NDIM>& idx)
{
    return dcomplex(idx(0), idx(1));
}

// Fill the ghost cell values of a cell-centered variable of type T on the
// coarsest level, and log the number of values that were received from other
// processes and how many of them are incorrect.
template <typename T>
void
fill_and_check_ghost_values(const std::string& name, Pointer<PatchHierarchy<NDIM>> patch_hierarchy)
{
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellVariable<NDIM, T>> var = new CellVariable<NDIM, T>(name);
    const int idx = var_db->registerVariableAndContext(var, var_db->getContext("context"), IntVector<NDIM>(1));
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);
    level->allocatePatchData(idx);

    // Start with ghost cell values that differ from the values in the cells
    // that they overlap.
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, T>> data = patch->getPatchData(idx);
        for (CellIterator<NDIM> ci(data->getGhostBox()); ci; ci++)
        {
            const T value = cell_value<T>(ci());
            (*data)(ci()) = patch->getBox().contains(ci()) ? value : (value == T() ? T(1) : T());
        }
    }

    RefineAlgorithm<NDIM> refine_alg;
    refine_alg.registerRefine(idx, idx, idx, nullptr);
    refine_alg.createSchedule(level)->fillData(0.0);

    int n_received = 0;
    int n_incorrect = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, T>> data = patch->getPatchData(idx);
        for (int src_patch_num = 0; src_patch_num < level->getNumberOfPatches(); ++src_patch_num)
        {
            if (level->getMappingForPatch(src_patch_num) == IBTK_MPI::getRank())
            {
                continue;
            }
            const Box<NDIM> received_box = data->getGhostBox() * level->getBoxForPatch(src_patch_num);
            for (CellIterator<NDIM> ci(received_box); ci; ci++)
            {
                ++n_received;
                if ((*data)(ci()) != cell_value<T>(ci()))
                {
                    ++n_incorrect;
                }
            }
        }
    }
    level->deallocatePatchData(idx);

    plog << name << " values received from other processes: " << IBTK_MPI::sumReduction(n_received)
         << ", incorrect: " << IBTK_MPI::sumReduction(n_incorrect) << "\n";
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    { // cleanup dynamically allocated objects prior to shutdown
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "patch_data_types.log");

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
        gridding_algorithm->makeCoarsestLevel(patch_hierarchy, 0.0);

        fill_and_check_ghost_values<bool>("bool", patch_hierarchy);
        fill_and_check_ghost_values<char>("char", patch_hierarchy);
        fill_and_check_ghost_values<dcomplex>("dcomplex", patch_hierarchy);
    } // cleanup dynamically allocated objects prior to shutdown
} // main
