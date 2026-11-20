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

#include <tbox/Logger.h>
#include <tbox/PIO.h>
#include <tbox/Pointer.h>

#include <Box.h>
#include <BoxArray.h>
#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellIterator.h>
#include <CellVariable.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <ProcessorMapping.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <VariableDatabase.h>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// Reset a refine schedule on one level that was built for an item whose source
// is its scratch data to an item whose source is other data. The schedule left
// out the copy from each patch to itself, so the reset must be rejected.

int
main(int argc, char** argv)
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<Logger::Appender> abort_append(new TestAppender());
    Logger::getInstance()->setAbortAppender(abort_append);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv);
    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));

    auto* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<CellVariable<NDIM, double>> var = new CellVariable<NDIM, double>("q");
    const int u_idx = var_db->registerVariableAndContext(var, var_db->getContext("u"), IntVector<NDIM>(1));
    const int v_idx = var_db->registerVariableAndContext(var, var_db->getContext("v"), IntVector<NDIM>(1));

    // Two patches split in the first direction.
    Pointer<PatchHierarchy<NDIM>> hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    const Box<NDIM> domain = grid_geometry->getPhysicalDomain()[0];
    BoxArray<NDIM> boxes(2);
    boxes[0] = domain;
    boxes[0].upper(0) = (domain.lower(0) + domain.upper(0)) / 2;
    boxes[1] = domain;
    boxes[1].lower(0) = boxes[0].upper(0) + 1;
    ProcessorMapping mapping(2);
    for (int k = 0; k < 2; ++k)
    {
        mapping.setProcessorAssignment(k, k % IBTK_MPI::getNodes());
    }
    hierarchy->makeNewPatchLevel(0, IntVector<NDIM>(1), boxes, mapping);
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    level->allocatePatchData(u_idx, 0.0);
    level->allocatePatchData(v_idx, 0.0);
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<CellData<NDIM, double>> u_data = level->getPatch(p())->getPatchData(u_idx);
        Pointer<CellData<NDIM, double>> v_data = level->getPatch(p())->getPatchData(v_idx);
        u_data->fillAll(-1.0);
        v_data->fillAll(13.0);
    }

    RefineAlgorithm<NDIM> in_place_alg;
    in_place_alg.registerRefine(u_idx, u_idx, u_idx, nullptr);
    Pointer<RefineSchedule<NDIM>> schedule = in_place_alg.createSchedule(level);

    RefineAlgorithm<NDIM> copy_alg;
    copy_alg.registerRefine(u_idx, v_idx, u_idx, nullptr);
    copy_alg.resetSchedule(schedule);
    schedule->fillData(0.0);

    int values = 0, values_set = 0;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<CellData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        for (CellIterator<NDIM> i(u_data->getGhostBox() * domain); i; i++)
        {
            ++values;
            values_set += (*u_data)(i()) == 13.0;
        }
    }
    values = IBTK_MPI::sumReduction(values);
    values_set = IBTK_MPI::sumReduction(values_set);
    plog << "reset accepted; values set " << values_set << " of " << values << '\n';
}
