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
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScMatUtilities.h>
#include <ibtk/PETScVecUtilities.h>

#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

// The RT0 prolongation matrix must not mix axes or depths for a multi-depth side variable. Build
// its diagonal-block preallocation for a depth-2 field, and confirm both that assembly needs no
// extra allocation and that each axis's depth is prolonged independently.
int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app);
    Pointer<PatchHierarchy<NDIM>> hierarchy = std::get<0>(hierarchy_tuple);
    Pointer<PatchLevel<NDIM>> coarse = hierarchy->getPatchLevel(0);
    Pointer<PatchLevel<NDIM>> fine = hierarchy->getPatchLevel(1);
    VariableDatabase<NDIM>* db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = db->getContext("rt0_depth");
    Pointer<SideVariable<NDIM, int>> dof_var = new SideVariable<NDIM, int>("dof", 2);
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u", 2);
    const int dof = db->registerVariableAndContext(dof_var, context, IntVector<NDIM>(1));
    const int u = db->registerVariableAndContext(u_var, context, IntVector<NDIM>(1));
    std::vector<int> counts[2];
    for (int ln = 0; ln < 2; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        level->allocatePatchData(dof);
        level->allocatePatchData(u);
        PETScVecUtilities::constructPatchLevelDOFIndices(counts[ln], dof, level);
    }
    AO ordering = nullptr;
    PETScVecUtilities::constructPatchLevelAO(ordering, counts[0], dof, coarse, 0);
    Mat P = nullptr;
    PETScMatUtilities::constructProlongationOp(P, "RT0", dof, counts[1], counts[0], fine, coarse, ordering, 0);
    MatInfo info = {};
    int ierr = MatGetInfo(P, MAT_GLOBAL_SUM, &info);
    IBTK_CHKERRQ(ierr);
    // PETSc counts reallocations during AIJ insertion, before final assembly, but only with logging
    // enabled. Report the count on standard output, because pout is mirrored into the compared log
    // file, so that the compared output does not depend on the PETSc configuration.
#if defined(PETSC_USE_LOG)
    if (IBTK_MPI::getRank() == 0) std::cout << "assembly reallocations = " << info.mallocs << '\n';
    if (info.mallocs != 0.0)
    {
        TBOX_ERROR("RT0 preallocation required additional storage.\n");
    }
#else
    if (IBTK_MPI::getRank() == 0) std::cout << "assembly reallocations unavailable without PETSc logging\n";
#endif
    PetscInt first = 0, last = 0;
    ierr = MatGetOwnershipRange(P, &first, &last);
    IBTK_CHKERRQ(ierr);
    for (PetscInt row = first; row < last; ++row)
    {
        PetscInt n = 0;
        ierr = MatGetRow(P, row, &n, nullptr, nullptr);
        IBTK_CHKERRQ(ierr);
        if (n < 1 || n > 2)
        {
            TBOX_ERROR("Unexpected RT0 row sparsity.\n");
        }
        ierr = MatRestoreRow(P, row, &n, nullptr, nullptr);
        IBTK_CHKERRQ(ierr);
    }
    Vec x = nullptr, y = nullptr, expected = nullptr;
    ierr = MatCreateVecs(P, &x, &y);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(y, &expected);
    IBTK_CHKERRQ(ierr);
    double error = 0.0;
    for (int component = 0; component < 2 * NDIM; ++component)
    {
        for (int ln = 0; ln < 2; ++ln)
        {
            Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
            for (PatchLevel<NDIM>::Iterator patch_it(level); patch_it; patch_it++)
            {
                Pointer<SideData<NDIM, double>> data = level->getPatch(patch_it())->getPatchData(u);
                data->fillAll(0.0);
                const int axis = component / 2;
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(data->getGhostBox(), axis)); b; b++)
                {
                    (*data)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower), component % 2) = component + 1.0;
                }
            }
        }
        PETScVecUtilities::copyToPatchLevelVec(x, u, dof, coarse);
        PETScVecUtilities::copyToPatchLevelVec(expected, u, dof, fine);
        ierr = MatMult(P, x, y);
        IBTK_CHKERRQ(ierr);
        ierr = VecAXPY(y, -1.0, expected);
        IBTK_CHKERRQ(ierr);
        double norm = 0.0;
        ierr = VecNorm(y, NORM_INFINITY, &norm);
        IBTK_CHKERRQ(ierr);
        if (!std::isfinite(norm))
        {
            TBOX_ERROR("Nonfinite RT0 depth error.\n");
        }
        error = std::max(error, norm);
    }
    plog << "axis/depth isolation error = " << error << '\n';
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&y);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&expected);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&P);
    IBTK_CHKERRQ(ierr);
    ierr = AODestroy(&ordering);
    IBTK_CHKERRQ(ierr);
    if (!(error < 1.0e-12))
    {
        TBOX_ERROR("RT0 axis/depth isolation error = " << error << " exceeds 1.0e-12.\n");
    }
    return 0;
}
