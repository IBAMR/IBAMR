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

#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/CartExtrapPhysBdryOp.h>
#include <ibtk/CartSideDoubleRT0Refine.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScMatUtilities.h>

#include <CartesianGridGeometry.h>
#include <CellData.h>
#include <CellVariable.h>
#include <RefineAlgorithm.h>
#include <RefineSchedule.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideIndex.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
// A separable profile gives independent one-dimensional interpolation weights.
double
profile_factor(int index, int axis, int direction, const Box<NDIM>& domain, const IntVector<NDIM>& periodic)
{
    if (!periodic(direction))
    {
        index = std::max(domain.lower()(direction), std::min(index, domain.upper()(direction) + (direction == axis)));
    }
    const double coordinate = index - domain.lower()(direction) + (direction == axis ? 0.0 : 0.5);
    const int extent = domain.upper()(direction) - domain.lower()(direction) + 1;
    return 1.0 + 0.1 * (direction + 1) * std::cos(2.0 * std::acos(-1.0) * coordinate / extent);
}

double
profile_value(const hier::Index<NDIM>& index,
              int axis,
              const Box<NDIM>& domain,
              const IntVector<NDIM>& periodic,
              bool normal_only)
{
    double value = axis + 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        if (!normal_only || d == axis)
        {
            value *= profile_factor(index(d), axis, d, domain, periodic);
        }
    }
    return value;
}

double
interpolated_value(const hier::Index<NDIM>& fine_index,
                   int axis,
                   const Box<NDIM>& domain,
                   const IntVector<NDIM>& periodic,
                   const IntVector<NDIM>& ratio,
                   bool rt0)
{
    double value = axis + 1.0;
    for (int d = 0; d < NDIM; ++d)
    {
        const double shift = d == axis ? 0.0 : 0.5;
        const double coordinate = (fine_index(d) + shift) / ratio(d) - shift;
        if (rt0 && d != axis)
        {
            const int coarse_index = static_cast<int>(std::floor(static_cast<double>(fine_index(d)) / ratio(d)));
            value *= profile_factor(coarse_index, axis, d, domain, periodic);
        }
        else
        {
            const int left = static_cast<int>(std::floor(coordinate));
            const double right_weight = coordinate - left;
            value *= (1.0 - right_weight) * profile_factor(left, axis, d, domain, periodic) +
                     right_weight * profile_factor(left + 1, axis, d, domain, periodic);
        }
    }
    return value;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    Pointer<Database> test_db = app->getInputDatabase()->getDatabase("test");
    const bool rt0 = test_db->getStringWithDefault("operator", "LINEAR") == "RT0";
    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app);
    Pointer<PatchHierarchy<NDIM>> hierarchy = std::get<0>(hierarchy_tuple);
    if (hierarchy->getFinestLevelNumber() != 1)
    {
        TBOX_ERROR("Side prolongation test requires two levels.\n");
    }
    Pointer<PatchLevel<NDIM>> coarse = hierarchy->getPatchLevel(0);
    Pointer<PatchLevel<NDIM>> fine = hierarchy->getPatchLevel(1);
    const Box<NDIM>& domain = coarse->getPhysicalDomain()[0];
    const IntVector<NDIM> ratio = fine->getRatio() / coarse->getRatio();
    Pointer<CartesianGridGeometry<NDIM>> geometry = hierarchy->getGridGeometry();
    const IntVector<NDIM> periodic = geometry->getPeriodicShift(coarse->getRatio());
    VariableDatabase<NDIM>* db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = db->getContext("side_prolongation");
    Pointer<SideVariable<NDIM, int>> u_dof = new SideVariable<NDIM, int>("u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof = new CellVariable<NDIM, int>("p_dof");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("p");
    Pointer<SideVariable<NDIM, double>> expected_var = new SideVariable<NDIM, double>("expected_u");
    const int u_dof_idx = db->registerVariableAndContext(u_dof, context, IntVector<NDIM>(1));
    const int p_dof_idx = db->registerVariableAndContext(p_dof, context, IntVector<NDIM>(1));
    const int u_idx = db->registerVariableAndContext(u_var, context, IntVector<NDIM>(1));
    const int p_idx = db->registerVariableAndContext(p_var, context, IntVector<NDIM>(1));
    const int expected_idx = db->registerVariableAndContext(expected_var, context, IntVector<NDIM>(1));
    for (int ln = 0; ln < 2; ++ln)
    {
        Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(ln);
        for (const int idx : { u_dof_idx, p_dof_idx, u_idx, p_idx, expected_idx })
        {
            level->allocatePatchData(idx);
        }
    }
    std::vector<int> coarse_counts, fine_counts;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(coarse_counts, u_dof_idx, p_dof_idx, coarse);
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(fine_counts, u_dof_idx, p_dof_idx, fine);
    if (*std::min_element(coarse_counts.begin(), coarse_counts.end()) <= 0 ||
        *std::min_element(fine_counts.begin(), fine_counts.end()) <= 0)
    {
        TBOX_ERROR("Side prolongation test requires DOFs on every rank at both levels.\n");
    }
    AO ordering = nullptr;
    int u_offset = 0, p_offset = 0;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelAO(
        ordering, coarse_counts, u_dof_idx, p_dof_idx, coarse, u_offset, p_offset);
    Mat prolongation = nullptr;
    PETScMatUtilities::constructProlongationOp(
        prolongation, rt0 ? "RT0" : "LINEAR", u_dof_idx, fine_counts, coarse_counts, fine, coarse, ordering, u_offset);
    Vec x = nullptr, result = nullptr, expected = nullptr;
    int ierr = MatCreateVecs(prolongation, &x, &result);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(result, &expected);
    IBTK_CHKERRQ(ierr);
    double reference_error = 0.0, pressure_leak = 0.0, schedule_error = 0.0, result_norm = 0.0;
    // Exercise each component separately, then pressure, then the RT0 subspace.
    for (int component = 0; component < NDIM + 2; ++component)
    {
        const bool normal_only = component == NDIM + 1;
        for (PatchLevel<NDIM>::Iterator p(coarse); p; p++)
        {
            Pointer<Patch<NDIM>> patch = coarse->getPatch(p());
            Pointer<SideData<NDIM, double>> u = patch->getPatchData(u_idx);
            Pointer<CellData<NDIM, double>> pressure = patch->getPatchData(p_idx);
            u->fillAll(0.0);
            pressure->fillAll(component == NDIM ? 2.0 : 0.0);
            for (int axis = 0; axis < NDIM; ++axis)
            {
                if (!normal_only && component != axis)
                {
                    continue;
                }
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(u->getGhostBox(), axis)); b; b++)
                {
                    (*u)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower)) =
                        profile_value(b(), axis, domain, periodic, normal_only);
                }
            }
        }
        for (PatchLevel<NDIM>::Iterator p(fine); p; p++)
        {
            Pointer<Patch<NDIM>> patch = fine->getPatch(p());
            Pointer<SideData<NDIM, double>> u = patch->getPatchData(expected_idx);
            Pointer<CellData<NDIM, double>> pressure = patch->getPatchData(p_idx);
            u->fillAll(0.0);
            pressure->fillAll(0.0);
            if (component < NDIM)
            {
                for (Box<NDIM>::Iterator b(SideGeometry<NDIM>::toSideBox(patch->getBox(), component)); b; b++)
                {
                    (*u)(SideIndex<NDIM>(b(), component, SideIndex<NDIM>::Lower)) =
                        interpolated_value(b(), component, domain, periodic, ratio, rt0);
                }
            }
        }
        if (normal_only)
        {
            // Tangentially constant data make LINEAR and RT0 refinement identical.
            Pointer<RefineAlgorithm<NDIM>> algorithm = new RefineAlgorithm<NDIM>();
            algorithm->registerRefine(expected_idx, u_idx, expected_idx, new CartSideDoubleRT0Refine());
            CartExtrapPhysBdryOp boundary(expected_idx, "CONSTANT");
            Pointer<RefineSchedule<NDIM>> schedule =
                algorithm->createSchedule(fine, Pointer<PatchLevel<NDIM>>(), 0, hierarchy, &boundary);
            schedule->fillData(0.0);
        }
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(x, u_idx, u_dof_idx, p_idx, p_dof_idx, coarse);
        IBAMR::StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
            expected, expected_idx, u_dof_idx, p_idx, p_dof_idx, fine);
        ierr = MatMult(prolongation, x, result);
        IBTK_CHKERRQ(ierr);
        double norm = 0.0;
        ierr = VecNorm(result, NORM_INFINITY, &norm);
        IBTK_CHKERRQ(ierr);
        if (!std::isfinite(norm))
        {
            TBOX_ERROR("Side prolongation produced a nonfinite norm.\n");
        }
        result_norm = std::max(result_norm, norm);
        ierr = VecAXPY(result, -1.0, expected);
        IBTK_CHKERRQ(ierr);
        ierr = VecNorm(result, NORM_INFINITY, &norm);
        IBTK_CHKERRQ(ierr);
        if (!std::isfinite(norm))
        {
            TBOX_ERROR("Side prolongation error has a nonfinite norm.\n");
        }
        if (normal_only)
        {
            schedule_error = norm;
        }
        else if (component == NDIM)
        {
            pressure_leak = norm;
        }
        else
        {
            reference_error = std::max(reference_error, norm);
        }
    }
    plog << std::setprecision(12) << "result norm = " << result_norm << '\n'
         << "interpolation error = " << reference_error << '\n'
         << "pressure leakage = " << pressure_leak << '\n'
         << "schedule error = " << schedule_error << '\n';
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&result);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&expected);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&prolongation);
    IBTK_CHKERRQ(ierr);
    ierr = AODestroy(&ordering);
    IBTK_CHKERRQ(ierr);
    return std::isfinite(reference_error) && std::isfinite(schedule_error) && result_norm > 0.0 &&
                   reference_error < 1.0e-12 && schedule_error < 1.0e-12 && pressure_leak < 1.0e-12 ?
               0 :
               1;
}
