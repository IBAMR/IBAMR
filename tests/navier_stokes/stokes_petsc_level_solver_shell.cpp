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

#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>
#include <ibamr/private/StaggeredStokesEigenSchurComplementShellBackend.h>

#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/private/PETScLevelSolverShellBackend.h>

#include <tbox/MemoryDatabase.h>

#include <CellData.h>
#include <CellVariable.h>
#include <LocationIndexRobinBcCoefs.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <numeric>

#include "../tests.h"

#include <ibamr/app_namespaces.h>

namespace
{
using CACell = std::array<int, NDIM>;
using CAFields = std::map<CACell, std::array<int, NDIM + 1>>;

/*! \brief Wrap a logical cell into the periodic test domain. */
CACell
ca_wrap(CACell cell, const int n)
{
    for (int& index : cell)
    {
        index = (index + n) % n;
    }
    return cell;
}

/*! \brief Return the pressure and incident velocities of one logical cell. */
std::set<int>
ca_cell_stencil(const CAFields& fields, const CACell& cell, const int n)
{
    const std::array<int, NDIM + 1>& local = fields.at(ca_wrap(cell, n));
    std::set<int> result(local.begin(), local.end());
    for (int axis = 0; axis < NDIM; ++axis)
    {
        CACell upper = cell;
        ++upper[axis];
        result.insert(fields.at(ca_wrap(upper, n))[axis]);
    }
    return result;
}

/*! \brief Return the two complete cells incident to a logical face. */
std::set<int>
ca_face_stencil(const CAFields& fields, CACell cell, const int axis, const int n)
{
    std::set<int> result = ca_cell_stencil(fields, cell, n);
    --cell[axis];
    const std::set<int> lower = ca_cell_stencil(fields, cell, n);
    result.insert(lower.begin(), lower.end());
    return result;
}

/*! \brief Read the existing PETSc subdomain query into ordered sets. */
std::vector<std::set<int>>
ca_read_sets(const std::vector<IS>& indices)
{
    std::vector<std::set<int>> result;
    for (IS is : indices)
    {
        PetscInt count = 0;
        const PetscInt* dofs = nullptr;
        int ierr = ISGetLocalSize(is, &count);
        IBTK_CHKERRQ(ierr);
        ierr = ISGetIndices(is, &dofs);
        IBTK_CHKERRQ(ierr);
        result.emplace_back(dofs, dofs + count);
        ierr = ISRestoreIndices(is, &dofs);
        IBTK_CHKERRQ(ierr);
    }
    return result;
}

/*! \brief Check velocity-seeded construction against independent logical stencils. */
int
check_ca_construction(Pointer<AppInitializer> app, Pointer<Database> test)
{
    const auto hierarchy_data = setup_hierarchy<NDIM>(app);
    Pointer<PatchHierarchy<NDIM>> hierarchy = std::get<0>(hierarchy_data);
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = variables->getContext("ca_test");
    Pointer<SideVariable<NDIM, int>> u_dof = new SideVariable<NDIM, int>("ca_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof = new CellVariable<NDIM, int>("ca_p_dof");
    const int udi = variables->registerVariableAndContext(u_dof, context, IntVector<NDIM>(1));
    const int pdi = variables->registerVariableAndContext(p_dof, context, IntVector<NDIM>(1));
    level->allocatePatchData(udi);
    level->allocatePatchData(pdi);
    std::vector<int> counts;
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(counts, udi, pdi, level);
    TBOX_ASSERT(counts.size() == 1);
    const int total = counts.front();
    const int n = app->getInputDatabase()->getInteger("N");
    const std::string scenario = test->getString("ca_construction");
    const bool sparse = scenario == "sparse_map_closure";
    CAFields fields;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, int>> u = patch->getPatchData(udi);
        Pointer<CellData<NDIM, int>> pressure = patch->getPatchData(pdi);
        if (sparse)
        {
            // A bijection deliberately interleaves pressure and velocity IDs.
            for (int axis = 0; axis < NDIM; ++axis)
            {
                const Box<NDIM> sides = SideGeometry<NDIM>::toSideBox(u->getGhostBox(), axis);
                for (Box<NDIM>::Iterator b(sides); b; b++)
                {
                    int& dof = (*u)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower));
                    if (dof >= 0)
                    {
                        dof = (5 * dof + 1) % total;
                    }
                }
            }
            for (Box<NDIM>::Iterator b(pressure->getGhostBox()); b; b++)
            {
                int& dof = (*pressure)(b());
                if (dof >= 0)
                {
                    dof = (5 * dof + 1) % total;
                }
            }
        }
        for (Box<NDIM>::Iterator b(patch->getBox()); b; b++)
        {
            CACell cell{};
            for (int d = 0; d < NDIM; ++d)
            {
                cell[d] = b()(d);
            }
            for (int axis = 0; axis < NDIM; ++axis)
            {
                fields[cell][axis] = (*u)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower));
            }
            fields[cell][NDIM] = (*pressure)(b());
        }
    }
    int failures = 0;
    Mat matrix = nullptr;
    int ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, total, total, 8, nullptr, &matrix);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(matrix, MAT_IGNORE_ZERO_ENTRIES, PETSC_FALSE);
    IBTK_CHKERRQ(ierr);
    for (int row = 0; row < total; ++row)
    {
        ierr = MatSetValue(matrix, row, row, 1.0, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
#if (NDIM == 2)
    const CouplingAwareASMSeedTraversalOrder default_order = CouplingAwareASMSeedTraversalOrder::I_J;
    const CouplingAwareASMSeedTraversalOrder alternate_order = CouplingAwareASMSeedTraversalOrder::J_I;
#else
    const CouplingAwareASMSeedTraversalOrder default_order = CouplingAwareASMSeedTraversalOrder::I_J_K;
    const CouplingAwareASMSeedTraversalOrder alternate_order = CouplingAwareASMSeedTraversalOrder::J_K_I;
#endif
    const auto assemble = [&]()
    {
        int error = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(error);
        error = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(error);
    };
    const auto verify_partition =
        [&](const std::vector<std::set<int>>& overlap, const std::vector<std::set<int>>& partition)
    {
        std::set<int> owned;
        if (overlap.size() != partition.size())
        {
            ++failures;
            return;
        }
        for (std::size_t k = 0; k < partition.size(); ++k)
        {
            for (const int dof : partition[k])
            {
                if (!overlap[k].count(dof) || !owned.insert(dof).second)
                {
                    ++failures;
                }
            }
        }
        if (owned.size() != static_cast<std::size_t>(total))
        {
            ++failures;
        }
    };
    std::vector<std::set<int>> overlap, partition;
    if (sparse)
    {
        const CACell origin{};
        CACell strong{}, weak{}, zero{};
        strong.fill(2);
        weak.fill(1);
        zero.fill(3);
        const int seed = fields.at(origin)[0];
        // Pressure entries must not influence either scale or stored-entry count.
        for (const auto& cell : fields)
        {
            ierr = MatSetValue(matrix, seed, cell.second[NDIM], 1.0e12, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
        for (const double tolerance : { 0.0, 1.0e-3 })
        {
            const double strong_value = tolerance == 0.0 ? 8 * std::numeric_limits<double>::epsilon() : 2.0e-3;
            // Four stored velocity entries, including zero: 3.5 epsilon is rejected.
            const double weak_value = tolerance == 0.0 ? 3.5 * std::numeric_limits<double>::epsilon() : 5.0e-4;
            ierr = MatSetValue(matrix, seed, fields.at(strong)[0], strong_value, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValue(matrix, seed, fields.at(weak)[0], weak_value, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetValue(matrix, seed, fields.at(zero)[0], 0.0, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
            assemble();
            StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains(
                overlap,
                partition,
                counts,
                udi,
                pdi,
                level,
                matrix,
                0,
                1,
                default_order,
                CouplingAwareASMClosurePolicy::RELAXED,
                tolerance);
            std::vector<std::set<int>> expected;
            for (const auto& cell : fields)
            {
                std::set<int> stencil = ca_face_stencil(fields, cell.first, 0, n);
                if (cell.first == origin)
                {
                    const std::set<int> extra = ca_face_stencil(fields, strong, 0, n);
                    stencil.insert(extra.begin(), extra.end());
                }
                expected.push_back(std::move(stencil));
            }
            if (overlap != expected)
            {
                ++failures;
            }
            verify_partition(overlap, partition);
        }
    }
    else if (scenario == "patch_level_maps")
    {
        assemble();
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (const CouplingAwareASMSeedTraversalOrder order : { default_order, alternate_order })
            {
                StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains(
                    overlap, partition, counts, udi, pdi, level, matrix, axis, 1, order);
                std::vector<CACell> ordered;
                for (const auto& cell : fields)
                {
                    ordered.push_back(cell.first);
                }
                if (order == alternate_order)
                {
                    std::sort(ordered.begin(),
                              ordered.end(),
                              [](const CACell& a, const CACell& b)
                              {
                                  for (int d = 0; d < NDIM; ++d)
                                  {
                                      const int axis = (d + 1) % NDIM;
                                      if (a[axis] != b[axis])
                                      {
                                          return a[axis] < b[axis];
                                      }
                                  }
                                  return false;
                              });
                }
                std::vector<std::set<int>> expected;
                for (const CACell& cell : ordered)
                {
                    expected.push_back(ca_face_stencil(fields, cell, axis, n));
                }
                if (overlap != expected)
                {
                    ++failures;
                }
                verify_partition(overlap, partition);
            }
        }
        StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains(
            overlap, partition, counts, udi, pdi, level, matrix, NDIM - 1, 2, default_order);
        std::vector<std::set<int>> expected;
        std::size_t count = 0;
        for (const auto& cell : fields)
        {
            if (count++ % 2 == 0)
            {
                expected.push_back(ca_face_stencil(fields, cell.first, NDIM - 1, n));
            }
        }
        if (overlap != expected)
        {
            ++failures;
        }
        verify_partition(overlap, partition);
    }
    else
    {
        // Same-component nearest-neighbor rows give one STRICT Vanka cell per seed.
        for (const auto& cell : fields)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                ierr = MatSetValue(matrix, cell.second[axis], cell.second[axis], 2 * NDIM + 1.0, INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
                for (int direction = 0; direction < NDIM; ++direction)
                {
                    for (const int sign : { -1, 1 })
                    {
                        CACell neighbor = cell.first;
                        neighbor[direction] += sign;
                        ierr = MatSetValue(
                            matrix, cell.second[axis], fields.at(ca_wrap(neighbor, n))[axis], -1.0, INSERT_VALUES);
                        IBTK_CHKERRQ(ierr);
                    }
                }
            }
        }
        assemble();
        for (const CouplingAwareASMClosurePolicy policy :
             { CouplingAwareASMClosurePolicy::RELAXED, CouplingAwareASMClosurePolicy::STRICT })
        {
            StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains(
                overlap, partition, counts, udi, pdi, level, matrix, 0, 1, default_order, policy);
            verify_partition(overlap, partition);
            if (policy == CouplingAwareASMClosurePolicy::STRICT)
            {
                std::vector<std::set<int>> expected;
                for (const auto& cell : fields)
                {
                    expected.push_back(ca_cell_stencil(fields, cell.first, n));
                }
                if (overlap != expected)
                {
                    ++failures;
                }
                std::vector<std::set<int>> geometric, geometric_partition;
                StaggeredStokesPETScMatUtilities::constructPatchLevelASMSubdomains(geometric,
                                                                                   geometric_partition,
                                                                                   IntVector<NDIM>(1),
                                                                                   IntVector<NDIM>(0),
                                                                                   counts,
                                                                                   udi,
                                                                                   pdi,
                                                                                   level,
                                                                                   nullptr);
                std::sort(geometric.begin(), geometric.end());
                std::sort(expected.begin(), expected.end());
                std::vector<std::set<int>> canonical = overlap;
                std::sort(canonical.begin(), canonical.end());
                if (canonical != geometric || canonical != expected)
                {
                    ++failures;
                }
            }
        }
        Pointer<SideVariable<NDIM, double>> u = new SideVariable<NDIM, double>("ca_u");
        Pointer<CellVariable<NDIM, double>> p = new CellVariable<NDIM, double>("ca_p");
        const int ui = variables->registerVariableAndContext(u, context, IntVector<NDIM>(1));
        const int pi = variables->registerVariableAndContext(p, context, IntVector<NDIM>(1));
        level->allocatePatchData(ui);
        level->allocatePatchData(pi);
        SAMRAIVectorReal<NDIM, double> x("ca_x", hierarchy, 0, 0), b("ca_b", hierarchy, 0, 0);
        x.addComponent(u, ui);
        x.addComponent(p, pi);
        b.addComponent(u, ui);
        b.addComponent(p, pi);
        x.setToScalar(0.0);
        // Dense velocity coupling closes the entire level in either policy.
        Mat changed_matrix = nullptr;
        ierr = MatDuplicate(matrix, MAT_COPY_VALUES, &changed_matrix);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetOption(changed_matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
        IBTK_CHKERRQ(ierr);
        for (const auto& row_cell : fields)
        {
            for (int row_axis = 0; row_axis < NDIM; ++row_axis)
            {
                const int row = row_cell.second[row_axis];
                for (const auto& column_cell : fields)
                {
                    for (int column_axis = 0; column_axis < NDIM; ++column_axis)
                    {
                        const int column = column_cell.second[column_axis];
                        ierr = MatSetValue(
                            changed_matrix, row, column, row == column ? 2 * NDIM + 1.0 : 0.01, INSERT_VALUES);
                        IBTK_CHKERRQ(ierr);
                    }
                }
            }
        }
        ierr = MatAssemblyBegin(changed_matrix, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        ierr = MatAssemblyEnd(changed_matrix, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(ierr);
        std::set<int> all_dofs;
        for (int dof = 0; dof < total; ++dof)
        {
            all_dofs.insert(dof);
        }
        const std::vector<std::set<int>> changed_overlap(fields.size(), all_dofs);
        std::vector<std::set<int>> changed_partition(fields.size());
        changed_partition.front() = all_dofs;
        for (const std::string pc_type : { "asm", "shell" })
        {
            for (const std::string policy : { "RELAXED", "STRICT" })
            {
                Pointer<MemoryDatabase> db = new MemoryDatabase("ca_solver");
                db->putString("pc_type", pc_type);
                db->putString("shell_pc_type", "multiplicative");
                db->putString("asm_subdomain_construction_mode", "COUPLING_AWARE");
                db->putString("coupling_aware_asm_closure_policy", policy);
                StaggeredStokesPETScLevelSolver solver("ca_solver_" + pc_type + "_" + policy, db, "ca_");
                std::vector<std::set<int>> expected, expected_partition;
                StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains(
                    expected,
                    expected_partition,
                    counts,
                    udi,
                    pdi,
                    level,
                    matrix,
                    0,
                    1,
                    default_order,
                    IBAMR::string_to_enum<CouplingAwareASMClosurePolicy>(policy));
                TBOX_ASSERT(expected != changed_overlap && expected_partition != changed_partition);
                // Replace the operator after teardown, then restore it in a third lifetime.
                for (int cycle = 0; cycle < 3; ++cycle)
                {
                    solver.setOperatorMat(cycle == 1 ? changed_matrix : matrix);
                    solver.initializeSolverState(x, b);
                    std::vector<IS>* actual_overlap = nullptr;
                    std::vector<IS>* actual_partition = nullptr;
                    solver.getASMSubdomains(&actual_partition, &actual_overlap);
                    const bool overlap_matches =
                        ca_read_sets(*actual_overlap) == (cycle == 1 ? changed_overlap : expected);
                    const bool partition_matches =
                        ca_read_sets(*actual_partition) == (cycle == 1 ? changed_partition : expected_partition);
                    if (!overlap_matches || !partition_matches)
                    {
                        pout << "CA lifetime mismatch: " << pc_type << " " << policy << " cycle=" << cycle
                             << " overlap_matches=" << overlap_matches << " partition_matches=" << partition_matches
                             << '\n';
                        ++failures;
                    }
                    solver.deallocateSolverState();
                }
            }
        }
        ierr = MatDestroy(&changed_matrix);
        IBTK_CHKERRQ(ierr);
        level->deallocatePatchData(ui);
        level->deallocatePatchData(pi);
    }
    ierr = MatDestroy(&matrix);
    IBTK_CHKERRQ(ierr);
    level->deallocatePatchData(udi);
    level->deallocatePatchData(pdi);
    plog << "failures = " << failures << '\n';
    return failures ? 1 : 0;
}

/*! \brief Check pressure-seeded patches against small periodic cell stencils. */
int
check_cav_construction(Pointer<AppInitializer> app, Pointer<Database> test)
{
    const auto hierarchy_data = setup_hierarchy<NDIM>(app);
    Pointer<PatchHierarchy<NDIM>> hierarchy = std::get<0>(hierarchy_data);
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = variables->getContext("cav_test");
    Pointer<SideVariable<NDIM, int>> u_dof = new SideVariable<NDIM, int>("cav_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof = new CellVariable<NDIM, int>("cav_p_dof");
    const int udi = variables->registerVariableAndContext(u_dof, context, IntVector<NDIM>(1));
    const int pdi = variables->registerVariableAndContext(p_dof, context, IntVector<NDIM>(1));
    level->allocatePatchData(udi);
    level->allocatePatchData(pdi);
    std::vector<int> counts;
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(counts, udi, pdi, level);
    std::vector<std::set<int>> patches;
    std::vector<int> seeds;
    const std::string scenario = test->getString("cav_construction");
    if (scenario == "parallel")
    {
        // Valid distributed DOF setup reaches the serial guard before matrix access.
        StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches(
            patches, seeds, counts, udi, pdi, level, nullptr);
        level->deallocatePatchData(udi);
        level->deallocatePatchData(pdi);
        return 0;
    }
    const int total = counts.front();
    const int n = app->getInputDatabase()->getInteger("N");
    CAFields fields;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        Pointer<SideData<NDIM, int>> u = patch->getPatchData(udi);
        Pointer<CellData<NDIM, int>> pressure = patch->getPatchData(pdi);
        // Interleave the full coupled IDs, including all shared/periodic ghosts.
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> sides = SideGeometry<NDIM>::toSideBox(u->getGhostBox(), axis);
            for (Box<NDIM>::Iterator b(sides); b; b++)
            {
                int& dof = (*u)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower));
                if (dof >= 0)
                {
                    dof = (5 * dof + 1) % total;
                }
            }
        }
        for (Box<NDIM>::Iterator b(pressure->getGhostBox()); b; b++)
        {
            int& dof = (*pressure)(b());
            if (dof >= 0)
            {
                dof = (5 * dof + 1) % total;
            }
        }
        for (Box<NDIM>::Iterator b(patch->getBox()); b; b++)
        {
            CACell cell{};
            for (int d = 0; d < NDIM; ++d)
            {
                cell[d] = b()(d);
            }
            for (int axis = 0; axis < NDIM; ++axis)
            {
                fields[cell][axis] = (*u)(SideIndex<NDIM>(b(), axis, SideIndex<NDIM>::Lower));
            }
            fields[cell][NDIM] = (*pressure)(b());
        }
    }
    int failures = 0;
    Mat elasticity = nullptr;
    int ierr = MatCreateSeqAIJ(PETSC_COMM_WORLD, total, total, 2 * NDIM + 3, nullptr, &elasticity);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(elasticity, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    IBTK_CHKERRQ(ierr);
    ierr = MatSetOption(elasticity, MAT_IGNORE_ZERO_ENTRIES, PETSC_FALSE);
    IBTK_CHKERRQ(ierr);
    const auto assemble = [&]()
    {
        int error = MatAssemblyBegin(elasticity, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(error);
        error = MatAssemblyEnd(elasticity, MAT_FINAL_ASSEMBLY);
        IBTK_CHKERRQ(error);
    };
    const CACell origin{};
    const int origin_velocity = fields.at(origin)[0];
    if (scenario == "pressure_row" || scenario == "pressure_column")
    {
        const int pressure = fields.at(origin)[NDIM];
        ierr = MatSetValue(elasticity,
                           scenario == "pressure_row" ? pressure : origin_velocity,
                           scenario == "pressure_row" ? origin_velocity : pressure,
                           1.0,
                           INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
        assemble();
        StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches(
            patches, seeds, counts, udi, pdi, level, elasticity);
        ierr = MatDestroy(&elasticity);
        IBTK_CHKERRQ(ierr);
        level->deallocatePatchData(udi);
        level->deallocatePatchData(pdi);
        return 0;
    }
    assemble();
#if (NDIM == 2)
    const CouplingAwareASMSeedTraversalOrder default_order = CouplingAwareASMSeedTraversalOrder::I_J;
    const CouplingAwareASMSeedTraversalOrder alternate_order = CouplingAwareASMSeedTraversalOrder::J_I;
#else
    const CouplingAwareASMSeedTraversalOrder default_order = CouplingAwareASMSeedTraversalOrder::I_J_K;
    const CouplingAwareASMSeedTraversalOrder alternate_order = CouplingAwareASMSeedTraversalOrder::J_K_I;
#endif
    std::vector<std::set<int>> standard;
    std::vector<int> expected_seeds;
    for (const auto& cell : fields)
    {
        standard.push_back(ca_cell_stencil(fields, cell.first, n));
        expected_seeds.push_back(cell.second[NDIM]);
    }
    // Exact outer-vector equality detects missing, duplicated, or reordered patches.
    for (const CouplingAwareASMClosurePolicy policy :
         { CouplingAwareASMClosurePolicy::RELAXED, CouplingAwareASMClosurePolicy::STRICT })
    {
        for (const CouplingAwareASMSeedTraversalOrder order : { default_order, alternate_order })
        {
            std::vector<CACell> ordered;
            for (const auto& cell : fields)
            {
                ordered.push_back(cell.first);
            }
            if (order == alternate_order)
            {
                std::sort(ordered.begin(),
                          ordered.end(),
                          [](const CACell& a, const CACell& b)
                          {
                              for (int d = 0; d < NDIM; ++d)
                              {
                                  const int axis = (d + 1) % NDIM;
                                  if (a[axis] != b[axis])
                                  {
                                      return a[axis] < b[axis];
                                  }
                              }
                              return false;
                          });
            }
            for (const int stride : { 1, 2 })
            {
                std::vector<int> ordered_seeds;
                std::vector<std::set<int>> expected;
                for (std::size_t k = 0; k < ordered.size(); k += stride)
                {
                    ordered_seeds.push_back(fields.at(ordered[k])[NDIM]);
                    expected.push_back(ca_cell_stencil(fields, ordered[k], n));
                }
                StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches(
                    patches, seeds, counts, udi, pdi, level, elasticity, stride, order, policy);
                if (patches != expected || seeds != ordered_seeds)
                {
                    ++failures;
                }
            }
        }
    }
    std::vector<std::set<int>> velocity_patches, partition;
    StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains(
        velocity_patches, partition, counts, udi, pdi, level, elasticity);
    if (velocity_patches == standard)
    {
        ++failures;
    }

    CACell remote{}, weak{};
    remote.fill(2);
    weak.fill(1);
    const int remote_velocity = fields.at(remote)[0];
    const int weak_velocity = fields.at(weak)[0];
    std::vector<std::set<int>> expected_relaxed;
    // A single distant face edge affects precisely the cells touching either end.
    for (const auto& cell : fields)
    {
        std::set<int> expected = ca_cell_stencil(fields, cell.first, n);
        const bool at_origin = expected.count(origin_velocity);
        const bool at_remote = expected.count(remote_velocity);
        if (at_origin || at_remote)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                for (const int sign : { -1, 1 })
                {
                    CACell neighbor = cell.first;
                    neighbor[axis] += sign;
                    const std::set<int> stencil = ca_cell_stencil(fields, neighbor, n);
                    expected.insert(stencil.begin(), stencil.end());
                }
            }
            const std::set<int> extra = ca_face_stencil(fields, at_origin ? remote : origin, 0, n);
            expected.insert(extra.begin(), extra.end());
        }
        expected_relaxed.push_back(std::move(expected));
    }
    for (const bool column_edge : { false, true })
    {
        ierr = MatZeroEntries(elasticity);
        IBTK_CHKERRQ(ierr);
        const int row = column_edge ? remote_velocity : origin_velocity;
        const int column = column_edge ? origin_velocity : remote_velocity;
        ierr = MatSetValue(elasticity, row, column, 1.0, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
        // Omit this row-relative weak edge, even though it is structurally present.
        ierr = MatSetValue(elasticity, row, weak_velocity, 5.0e-4, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
        assemble();
        for (const CouplingAwareASMClosurePolicy policy :
             { CouplingAwareASMClosurePolicy::RELAXED, CouplingAwareASMClosurePolicy::STRICT })
        {
            StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches(
                patches, seeds, counts, udi, pdi, level, elasticity, 1, default_order, policy, 1.0e-3);
            if (seeds != expected_seeds ||
                patches != (policy == CouplingAwareASMClosurePolicy::RELAXED ? expected_relaxed : standard))
            {
                ++failures;
            }
        }
    }
    if (expected_relaxed == standard || !expected_relaxed.front().count(remote_velocity) ||
        standard.front().count(remote_velocity))
    {
        ++failures;
    }
    // Changing the same borrowed matrix must remove the previous expansion.
    ierr = MatZeroEntries(elasticity);
    IBTK_CHKERRQ(ierr);
    assemble();
    StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches(
        patches, seeds, counts, udi, pdi, level, elasticity);
    if (patches != standard || seeds != expected_seeds)
    {
        ++failures;
    }
    const std::set<int> remote_patch = ca_cell_stencil(fields, remote, n);
    for (const int dof : remote_patch)
    {
        if (dof != fields.at(remote)[NDIM])
        {
            ierr = MatSetValue(elasticity, origin_velocity, dof, 1.0, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }
    // A second-hop edge must not expand the origin seed recursively.
    ierr = MatSetValue(elasticity, remote_velocity, weak_velocity, 1.0, INSERT_VALUES);
    IBTK_CHKERRQ(ierr);
    assemble();
    StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches(
        patches, seeds, counts, udi, pdi, level, elasticity, 1, default_order, CouplingAwareASMClosurePolicy::STRICT);
    std::set<int> complete = standard.front();
    complete.insert(remote_patch.begin(), remote_patch.end());
    if (patches.size() != standard.size() || seeds != expected_seeds || patches.front() != complete)
    {
        ++failures;
    }
    // In RELAXED, the one-hop stencil closes neighboring cells around both cells.
    std::set<int> complete_relaxed = complete;
    for (const CACell center : { origin, remote })
    {
        for (int axis = 0; axis < NDIM; ++axis)
        {
            for (const int sign : { -1, 1 })
            {
                CACell neighbor = center;
                neighbor[axis] += sign;
                const std::set<int> stencil = ca_cell_stencil(fields, neighbor, n);
                complete_relaxed.insert(stencil.begin(), stencil.end());
            }
        }
    }
    StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches(
        patches, seeds, counts, udi, pdi, level, elasticity);
    if (patches.size() != standard.size() || seeds != expected_seeds || patches.front() != complete_relaxed)
    {
        ++failures;
    }

    PoissonSpecifications coefficients("cav_mac");
    coefficients.setCConstant(1.0);
    coefficients.setDConstant(-1.0);
    std::vector<RobinBcCoefStrategy<NDIM>*> boundary(NDIM, nullptr);
    Mat mac = nullptr;
    StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(
        mac, coefficients, boundary, 0.0, counts, udi, pdi, level);
    double pressure_row_error = 0.0;
    for (const auto& cell : fields)
    {
        std::map<int, double> expected;
        for (int axis = 0; axis < NDIM; ++axis)
        {
            CACell upper = cell.first;
            ++upper[axis];
            expected[cell.second[axis]] = n;
            expected[fields.at(ca_wrap(upper, n))[axis]] = -n;
        }
        PetscInt count = 0;
        const PetscInt* columns = nullptr;
        const PetscScalar* values = nullptr;
        ierr = MatGetRow(mac, cell.second[NDIM], &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < count; ++k)
        {
            pressure_row_error =
                std::max(pressure_row_error, static_cast<double>(PetscAbsScalar(values[k] - expected[columns[k]])));
            expected.erase(columns[k]);
        }
        for (const auto& missing : expected)
        {
            pressure_row_error = std::max(pressure_row_error, std::abs(missing.second));
        }
        ierr = MatRestoreRow(mac, cell.second[NDIM], &count, &columns, &values);
        IBTK_CHKERRQ(ierr);
    }
    plog << "pressure_row_error = " << pressure_row_error << '\n';
#if (NDIM == 3)
    Vec pressure_shift = nullptr, action = nullptr;
    ierr = MatCreateVecs(mac, &pressure_shift, &action);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(pressure_shift, 0.0);
    IBTK_CHKERRQ(ierr);
    for (const auto& cell : fields)
    {
        ierr = VecSetValue(pressure_shift, cell.second[NDIM], 1.0, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(pressure_shift);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(pressure_shift);
    IBTK_CHKERRQ(ierr);
    ierr = VecNormalize(pressure_shift, nullptr);
    IBTK_CHKERRQ(ierr);
    ierr = MatMult(mac, pressure_shift, action);
    IBTK_CHKERRQ(ierr);
    PetscReal shift_error = 0.0;
    ierr = VecNorm(action, NORM_INFINITY, &shift_error);
    IBTK_CHKERRQ(ierr);
    MatNullSpace nullspace = nullptr;
    ierr = MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_FALSE, 1, &pressure_shift, &nullspace);
    IBTK_CHKERRQ(ierr);
    PetscBool valid = PETSC_FALSE;
    ierr = MatNullSpaceTest(nullspace, mac, &valid);
    IBTK_CHKERRQ(ierr);
    if (!valid)
    {
        ++failures;
    }
    plog << "pressure_shift_error = " << shift_error << '\n';
    ierr = MatNullSpaceDestroy(&nullspace);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&pressure_shift);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&action);
    IBTK_CHKERRQ(ierr);
#endif
    ierr = MatDestroy(&mac);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&elasticity);
    IBTK_CHKERRQ(ierr);
    level->deallocatePatchData(udi);
    level->deallocatePatchData(pdi);
    plog << "pressure_patches = " << standard.size() << '\n';
    plog << "failures = " << failures << '\n';
    return failures ? 1 : 0;
}

double
norm_inf(Vec x)
{
    PetscReal value = 0.0;
    const int ierr = VecNorm(x, NORM_INFINITY, &value);
    IBTK_CHKERRQ(ierr);
    return value;
}

struct RowMatrix
{
    std::vector<std::vector<PetscInt>> columns{ { 0, 1 }, { 0, 1, 2 }, { 1, 2, 3 }, { 2, 3 } };
    std::vector<std::vector<PetscScalar>> values{ { 9, 9 }, { 9, 9, 9 }, { 9, 9, 9 }, { 9, 9 } };
    int row_reads = 0, multiplies = 0;
};

PetscErrorCode
get_test_row(Mat mat, PetscInt row, PetscInt* n, const PetscInt** columns, const PetscScalar** values)
{
    RowMatrix* context = nullptr;
    int ierr = MatShellGetContext(mat, &context);
    IBTK_CHKERRQ(ierr);
    ++context->row_reads;
    *n = static_cast<PetscInt>(context->columns[row].size());
    *columns = context->columns[row].data();
    if (values)
    {
        *values = context->values[row].data();
    }
    return 0;
}

PetscErrorCode
restore_test_row(Mat, PetscInt, PetscInt*, const PetscInt**, const PetscScalar**)
{
    return 0;
}

PetscErrorCode
multiply_test_matrix(Mat mat, Vec x, Vec y)
{
    RowMatrix* context = nullptr;
    int ierr = MatShellGetContext(mat, &context);
    IBTK_CHKERRQ(ierr);
    ++context->multiplies;
    const PetscScalar* input = nullptr;
    PetscScalar* output = nullptr;
    ierr = VecGetArrayRead(x, &input);
    IBTK_CHKERRQ(ierr);
    ierr = VecGetArray(y, &output);
    IBTK_CHKERRQ(ierr);
    for (std::size_t row = 0; row < context->columns.size(); ++row)
    {
        output[row] = 0.0;
        for (std::size_t j = 0; j < context->columns[row].size(); ++j)
        {
            output[row] += context->values[row][j] * input[context->columns[row][j]];
        }
    }
    ierr = VecRestoreArray(y, &output);
    IBTK_CHKERRQ(ierr);
    ierr = VecRestoreArrayRead(x, &input);
    IBTK_CHKERRQ(ierr);
    return 0;
}

// Prescribed local corrections isolate the real shared composer's residual action.
class StageBackend : public PETScLevelSolverShellBackend
{
public:
    /*! \brief Record residual samples in test-owned storage. */
    StageBackend(std::vector<PetscScalar>& samples, std::vector<std::size_t>& visits, std::size_t stages)
        : d_samples(samples), d_visits(visits), d_stages(stages)
    {
    }
    ~StageBackend() override
    {
        deallocateSolverState();
    }
    void
    initializeSolverState(Mat mat,
                          Vec x,
                          Vec b,
                          const std::vector<IS>&,
                          const std::vector<IS>&,
                          const std::string&,
                          bool multiplicative = false,
                          PETScLevelSolverShellTraversal traversal = PETScLevelSolverShellTraversal::FORWARD) override
    {
        deallocateSolverState();
        initializeComposition(mat, x, b, multiplicative, traversal);
        finalizeComposition();
    }
    void deallocateSolverState() override
    {
        deallocateComposition();
    }

protected:
    std::size_t getNumberOfSubdomains() const override
    {
        return d_stages;
    }
    void beginSubdomainRhs(std::size_t i, Vec source) override
    {
        PetscScalar value = 0.0;
        int ierr = VecGetValues(source, 1, d_dofs[i].data(), &value);
        IBTK_CHKERRQ(ierr);
        d_samples.push_back(value);
        d_visits.push_back(i);
    }
    void endSubdomainRhs(std::size_t, Vec) override
    {
    }
    void solveSubdomain(std::size_t) override
    {
    }
    void accumulateSubdomainCorrection(std::size_t i, Vec y) override
    {
        int ierr = VecSetValue(y, d_dofs[i][0], 1.0, ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyBegin(y);
        IBTK_CHKERRQ(ierr);
        ierr = VecAssemblyEnd(y);
        IBTK_CHKERRQ(ierr);
    }
    const std::vector<PetscInt>& getSubdomainCorrectionDofs(std::size_t i) const override
    {
        return d_dofs[i];
    }
    void copySubdomainCorrection(std::size_t, PetscScalar* values) override
    {
        values[0] = 1.0;
    }

private:
    std::vector<PetscScalar>& d_samples;
    std::vector<std::size_t>& d_visits;
    const std::size_t d_stages;
    const std::vector<std::vector<PetscInt>> d_dofs{ { 0 }, { 1 }, { 2 } };
};

int
check_stages(const bool fallback,
             const bool invalidate,
             const PETScLevelSolverShellTraversal traversal = PETScLevelSolverShellTraversal::FORWARD,
             const std::size_t stages = 3)
{
    RowMatrix context;
    Mat mat = nullptr;
    Vec x = nullptr, y = nullptr;
    int ierr = MatCreateShell(PETSC_COMM_SELF, 4, 4, 4, 4, &context, &mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatShellSetOperation(mat, MATOP_MULT, reinterpret_cast<void (*)(void)>(multiply_test_matrix));
    IBTK_CHKERRQ(ierr);
    if (!fallback)
    {
        ierr = MatShellSetOperation(mat, MATOP_GET_ROW, reinterpret_cast<void (*)(void)>(get_test_row));
        IBTK_CHKERRQ(ierr);
        ierr = MatShellSetOperation(mat, MATOP_RESTORE_ROW, reinterpret_cast<void (*)(void)>(restore_test_row));
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecCreateSeq(PETSC_COMM_SELF, 4, &x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(x, &y);
    IBTK_CHKERRQ(ierr);
    const PetscInt indices[] = { 0, 1, 2, 3 };
    const PetscScalar rhs[] = { 10, 20, 30, 40 };
    ierr = VecSetValues(x, 4, indices, rhs, INSERT_VALUES);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(x);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    IBTK_CHKERRQ(ierr);
    std::vector<PetscScalar> samples;
    std::vector<std::size_t> visits;
    StageBackend backend(samples, visits, stages);
    std::vector<std::size_t> expected_visits;
    std::vector<PetscScalar> expected_samples;
    std::vector<PetscScalar> expected_correction(4, 0.0);
    int expected_reads = 0;
    if (stages == 1)
    {
        expected_visits = { 0 };
        expected_samples = { 10 };
        expected_correction[0] = 1;
    }
    else if (stages == 3)
    {
        if (traversal == PETScLevelSolverShellTraversal::FORWARD)
        {
            expected_visits = { 0, 1, 2 };
            expected_samples = { 10, 21, 31 };
            expected_correction = { 1, 1, 1, 0 };
            expected_reads = 5;
        }
        else if (traversal == PETScLevelSolverShellTraversal::REVERSE)
        {
            expected_visits = { 2, 1, 0 };
            expected_samples = { 30, 21, 11 };
            expected_correction = { 1, 1, 1, 0 };
            expected_reads = 6;
        }
        else
        {
            expected_visits = { 0, 1, 2, 1, 0 };
            expected_samples = { 10, 21, 31, 20, 10 };
            expected_correction = { 2, 2, 1, 0 };
            expected_reads = 11;
        }
    }
    const int expected_multiplies = expected_visits.empty() ? 0 : static_cast<int>(expected_visits.size()) - 1;
    int failures = 0;
    for (int cycle = 0; cycle < 2; ++cycle)
    {
        backend.initializeSolverState(mat, x, x, {}, {}, "", true, traversal);
        // Application must read values at cached offsets, not an update-matrix copy.
        context.values = { { 2, -1 }, { -1, 2, -1 }, { -1, 2, -1 }, { -1, 2 } };
        if (invalidate)
        {
            std::swap(context.columns[0][0], context.columns[0][1]);
        }
        for (int application = 0; application < 2; ++application)
        {
            samples.clear();
            visits.clear();
            context.row_reads = context.multiplies = 0;
            backend.apply(x, y);
            if (invalidate)
            {
                backend.deallocateSolverState();
                ierr = MatDestroy(&mat);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&x);
                IBTK_CHKERRQ(ierr);
                ierr = VecDestroy(&y);
                IBTK_CHKERRQ(ierr);
                return 0;
            }
            if (visits != expected_visits || samples != expected_samples ||
                context.row_reads != (fallback ? 0 : expected_reads) ||
                context.multiplies != (fallback ? expected_multiplies : 0))
            {
                ++failures;
            }
            const PetscScalar* values = nullptr;
            ierr = VecGetArrayRead(y, &values);
            IBTK_CHKERRQ(ierr);
            if (!std::equal(expected_correction.begin(), expected_correction.end(), values))
            {
                ++failures;
            }
            if (stages == 3 && cycle == 0 && application == 0)
            {
                plog << "stage_rhs =";
                for (PetscScalar sample : samples)
                {
                    plog << ' ' << sample;
                }
                plog << "\ncorrection = " << values[0] << ' ' << values[1] << ' ' << values[2] << ' ' << values[3]
                     << "\nrow_reads = " << context.row_reads << "\nmatmult_calls = " << context.multiplies << '\n';
            }
            ierr = VecRestoreArrayRead(y, &values);
            IBTK_CHKERRQ(ierr);
        }
        backend.deallocateSolverState();
        backend.deallocateSolverState();
    }
    ierr = MatDestroy(&mat);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&y);
    IBTK_CHKERRQ(ierr);
    return failures;
}

int
check_hand_solve(const bool blas,
                 const PETScLevelSolverShellTraversal traversal = PETScLevelSolverShellTraversal::FORWARD,
                 const bool traversal_case = false,
                 const std::string& eigen_backend = "",
                 const bool additive = false)
{
    const int rank = IBTK_MPI::getRank();
    Mat mat = nullptr;
    Vec x = nullptr, y = nullptr, expected = nullptr;
    int ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 3, 3, 3, nullptr, 3, nullptr, &mat);
    IBTK_CHKERRQ(ierr);
    PetscInt lo = 0, hi = 0;
    ierr = MatGetOwnershipRange(mat, &lo, &hi);
    IBTK_CHKERRQ(ierr);
    for (PetscInt row = lo; row < hi; ++row)
    {
        for (PetscInt col = std::max<PetscInt>(0, row - 1); col <= std::min<PetscInt>(2, row + 1); ++col)
        {
            ierr = MatSetValue(mat, row, col, row == col ? 2.0 : -1.0, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatCreateVecs(mat, &x, &y);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(y, &expected);
    IBTK_CHKERRQ(ierr);
    for (PetscInt row = lo; row < hi; ++row)
    {
        ierr = VecSetValue(x, row, row + 1.0, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
        const bool parallel = traversal_case && IBTK_MPI::getNodes() == 2;
        const std::vector<PetscScalar> target =
            additive ? std::vector<PetscScalar>{ 4.0 / 3.0, 5.0 / 3.0, 8.0 / 3.0 } :
            traversal == PETScLevelSolverShellTraversal::REVERSE ?
                       (parallel ? std::vector<PetscScalar>{ 40.0 / 9.0, 41.0 / 9.0, 8.0 / 3.0 } :
                                   std::vector<PetscScalar>{ 20.0 / 9.0, 31.0 / 9.0, 8.0 / 3.0 }) :
            traversal == PETScLevelSolverShellTraversal::SYMMETRIC ?
                       (parallel ? std::vector<PetscScalar>{ 64.0 / 27.0, 107.0 / 27.0, 32.0 / 9.0 } :
                                   std::vector<PetscScalar>{ 64.0 / 27.0, 101.0 / 27.0, 28.0 / 9.0 }) :
                       (parallel ? std::vector<PetscScalar>{ 8.0 / 3.0, 37.0 / 9.0, 32.0 / 9.0 } :
                                   std::vector<PetscScalar>{ 4.0 / 3.0, 29.0 / 9.0, 28.0 / 9.0 });
        ierr = VecSetValue(expected, row, target[row], INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(x);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(x);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyBegin(expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(expected);
    IBTK_CHKERRQ(ierr);
    // Rank 1 contributes only to stage zero; traversal cases duplicate rank zero
    // there to distinguish the summed parallel correction from the serial result.
    const std::vector<std::vector<PetscInt>> dofs =
        rank == 0 ? std::vector<std::vector<PetscInt>>{ { 0, 1 }, { 1, 2 } } :
                    std::vector<std::vector<PetscInt>>{ traversal_case ? std::vector<PetscInt>{ 0, 1 } :
                                                                         std::vector<PetscInt>{ 1, 2 } };
    std::vector<IS> overlap(dofs.size()), partition(dofs.size());
    for (std::size_t i = 0; i < dofs.size(); ++i)
    {
        ierr = ISCreateGeneral(PETSC_COMM_SELF, 2, dofs[i].data(), PETSC_COPY_VALUES, &overlap[i]);
        IBTK_CHKERRQ(ierr);
        ierr = ISCreateGeneral(PETSC_COMM_SELF,
                               additive ? (i == 0 ? 2 : 1) : 0,
                               additive ? dofs[i].data() + (i == 0 ? 0 : 1) : nullptr,
                               PETSC_COPY_VALUES,
                               &partition[i]);
        IBTK_CHKERRQ(ierr);
    }
    Pointer<MemoryDatabase> db = new MemoryDatabase("hand");
    db->putString("blas_lapack_subdomain_solver_type", "lu");
    std::unique_ptr<PETScLevelSolverShellBackend> backend =
        PETScLevelSolverShellBackendManager::get_manager().allocateBackend(
            eigen_backend.empty() ? (blas ? "blas-lapack" : "petsc") : eigen_backend, db);
    int failures = 0;
    for (int cycle = 0; cycle < 2; ++cycle)
    {
        backend->initializeSolverState(mat, x, y, overlap, partition, "r09_hand", !additive, traversal);
        for (int application = 0; application < 2; ++application)
        {
            backend->apply(x, y);
            ierr = VecAXPY(y, -1.0, expected);
            IBTK_CHKERRQ(ierr);
            const double error = norm_inf(y);
            if (!std::isfinite(error) || error > 1.0e-12)
            {
                ++failures;
            }
            if (cycle == 0 && application == 0)
            {
                plog << "hand_error = " << error << '\n';
            }
        }
        backend->deallocateSolverState();
        backend->deallocateSolverState();
    }
    for (IS& is : overlap)
    {
        ierr = ISDestroy(&is);
        IBTK_CHKERRQ(ierr);
    }
    for (IS& is : partition)
    {
        ierr = ISDestroy(&is);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatDestroy(&mat);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&x);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&y);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&expected);
    IBTK_CHKERRQ(ierr);
    return failures;
}

// Gather the RHS independently of the backend's restriction/prolongation
// scatters. Additive writes use the partition; multiplicative visits add full corrections.
void
reference_action(Mat mat,
                 Vec rhs,
                 Vec result,
                 const std::vector<IS>& overlap,
                 const std::vector<IS>& partition,
                 const bool multiplicative,
                 const PETScLevelSolverShellTraversal traversal)
{
    Vec residual = nullptr, gathered = nullptr;
    VecScatter gather = nullptr;
    int ierr = VecDuplicate(rhs, &residual);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterCreateToAll(rhs, &gather, &gathered);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(result, 0.0);
    IBTK_CHKERRQ(ierr);
    Mat* submat = nullptr;
    ierr = MatCreateSubMatrices(
        mat, static_cast<PetscInt>(overlap.size()), overlap.data(), overlap.data(), MAT_INITIAL_MATRIX, &submat);
    IBTK_CHKERRQ(ierr);
    if (!multiplicative)
    {
        ierr = VecScatterBegin(gather, rhs, gathered, INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
        ierr = VecScatterEnd(gather, rhs, gathered, INSERT_VALUES, SCATTER_FORWARD);
        IBTK_CHKERRQ(ierr);
    }
    const int n_stages = IBTK_MPI::maxReduction(static_cast<int>(overlap.size()));
    std::vector<int> visits(n_stages);
    std::iota(visits.begin(), visits.end(), 0);
    if (multiplicative && traversal == PETScLevelSolverShellTraversal::REVERSE)
    {
        std::reverse(visits.begin(), visits.end());
    }
    else if (multiplicative && traversal == PETScLevelSolverShellTraversal::SYMMETRIC)
    {
        for (int i = n_stages - 2; i >= 0; --i)
        {
            visits.push_back(i);
        }
    }
    for (const int i : visits)
    {
        if (multiplicative)
        {
            // Recompute the global original residual independently before each stage.
            ierr = MatMult(mat, result, residual);
            IBTK_CHKERRQ(ierr);
            ierr = VecAYPX(residual, -1.0, rhs);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterBegin(gather, residual, gathered, INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
            ierr = VecScatterEnd(gather, residual, gathered, INSERT_VALUES, SCATTER_FORWARD);
            IBTK_CHKERRQ(ierr);
        }
        if (i < static_cast<int>(overlap.size()))
        {
            PetscInt n = 0, m = 0;
            const PetscInt* indices = nullptr;
            const PetscInt* owned = nullptr;
            ierr = ISGetLocalSize(overlap[i], &n);
            IBTK_CHKERRQ(ierr);
            if (!multiplicative)
            {
                ierr = ISGetLocalSize(partition[i], &m);
                IBTK_CHKERRQ(ierr);
            }
            ierr = ISGetIndices(overlap[i], &indices);
            IBTK_CHKERRQ(ierr);
            if (!multiplicative)
            {
                ierr = ISGetIndices(partition[i], &owned);
                IBTK_CHKERRQ(ierr);
            }
            Vec local_rhs = nullptr, local_solution = nullptr;
            ierr = MatCreateVecs(submat[i], &local_solution, &local_rhs);
            IBTK_CHKERRQ(ierr);
            const PetscScalar* global_values = nullptr;
            PetscScalar* local_values = nullptr;
            ierr = VecGetArrayRead(gathered, &global_values);
            IBTK_CHKERRQ(ierr);
            ierr = VecGetArray(local_rhs, &local_values);
            IBTK_CHKERRQ(ierr);
            for (PetscInt j = 0; j < n; ++j)
            {
                local_values[j] = global_values[indices[j]];
            }
            ierr = VecRestoreArray(local_rhs, &local_values);
            IBTK_CHKERRQ(ierr);
            ierr = VecRestoreArrayRead(gathered, &global_values);
            IBTK_CHKERRQ(ierr);
            KSP ksp = nullptr;
            PC pc = nullptr;
            ierr = KSPCreate(PETSC_COMM_SELF, &ksp);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetOperators(ksp, submat[i], submat[i]);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSetType(ksp, KSPPREONLY);
            IBTK_CHKERRQ(ierr);
            ierr = KSPGetPC(ksp, &pc);
            IBTK_CHKERRQ(ierr);
            ierr = PCSetType(pc, PCSVD);
            IBTK_CHKERRQ(ierr);
            ierr = KSPSolve(ksp, local_rhs, local_solution);
            IBTK_CHKERRQ(ierr);
            const PetscScalar* solution_values = nullptr;
            ierr = VecGetArrayRead(local_solution, &solution_values);
            IBTK_CHKERRQ(ierr);
            for (PetscInt j = 0; j < (multiplicative ? n : m); ++j)
            {
                const PetscInt position =
                    multiplicative ? j :
                                     static_cast<PetscInt>(std::lower_bound(indices, indices + n, owned[j]) - indices);
                TBOX_ASSERT(position < n && (multiplicative || indices[position] == owned[j]));
                ierr = VecSetValue(result,
                                   multiplicative ? indices[j] : owned[j],
                                   solution_values[position],
                                   multiplicative ? ADD_VALUES : INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            ierr = VecRestoreArrayRead(local_solution, &solution_values);
            IBTK_CHKERRQ(ierr);
            if (!multiplicative)
            {
                ierr = ISRestoreIndices(partition[i], &owned);
                IBTK_CHKERRQ(ierr);
            }
            ierr = ISRestoreIndices(overlap[i], &indices);
            IBTK_CHKERRQ(ierr);
            ierr = KSPDestroy(&ksp);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&local_rhs);
            IBTK_CHKERRQ(ierr);
            ierr = VecDestroy(&local_solution);
            IBTK_CHKERRQ(ierr);
        }
        if (multiplicative)
        {
            ierr = VecAssemblyBegin(result);
            IBTK_CHKERRQ(ierr);
            ierr = VecAssemblyEnd(result);
            IBTK_CHKERRQ(ierr);
        }
    }
    ierr = VecAssemblyBegin(result);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(result);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroyMatrices(static_cast<PetscInt>(overlap.size()), &submat);
    IBTK_CHKERRQ(ierr);
    ierr = VecScatterDestroy(&gather);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&gathered);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&residual);
    IBTK_CHKERRQ(ierr);
}

// The application fixture has one distant face edge initially and a complete
// remote velocity stencil after reinitialization. Enumerate logical cell supports
// independently of the production geometry owner and its global-index closure.
std::vector<std::set<int>>
cav_application_patches(const CAFields& fields, const int n, const bool strict, const int cycle, Mat elasticity)
{
    const CACell origin{};
    CACell remote{};
    remote.fill(2);
    const int source = fields.at(origin)[0];
    std::set<int> targets{ fields.at(remote)[0] };
    if (cycle == 1)
    {
        targets = ca_cell_stencil(fields, remote, n);
        targets.erase(fields.at(remote)[NDIM]);
    }
    int ierr = MatZeroEntries(elasticity);
    IBTK_CHKERRQ(ierr);
    for (const int target : targets)
    {
        // A symmetric positive semidefinite spring couples the two velocities.
        ierr = MatSetValue(elasticity, source, source, 0.25, ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetValue(elasticity, target, target, 0.25, ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetValue(elasticity, source, target, -0.25, ADD_VALUES);
        IBTK_CHKERRQ(ierr);
        ierr = MatSetValue(elasticity, target, source, -0.25, ADD_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(elasticity, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(elasticity, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    std::vector<std::set<int>> patches;
    for (const auto& seed : fields)
    {
        std::set<int> patch = ca_cell_stencil(fields, seed.first, n);
        std::set<int> velocities = patch;
        velocities.erase(seed.second[NDIM]);
        const std::set<int> original = velocities;
        if (original.count(source))
        {
            velocities.insert(targets.begin(), targets.end());
        }
        if (std::any_of(targets.begin(), targets.end(), [&](int dof) { return original.count(dof); }))
        {
            velocities.insert(source);
        }
        if (velocities != original)
        {
            for (const auto& candidate : fields)
            {
                const std::set<int> stencil = ca_cell_stencil(fields, candidate.first, n);
                int incident = 0;
                for (const int dof : stencil)
                {
                    incident += velocities.count(dof);
                }
                if (strict ? incident == 2 * NDIM : incident > 0)
                {
                    patch.insert(stencil.begin(), stencil.end());
                }
            }
        }
        patches.push_back(std::move(patch));
    }
    return patches;
}

class FieldCheckingSolver : public StaggeredStokesPETScLevelSolver
{
public:
    using StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver;
    void checkFields(const std::set<int>& velocity, const std::set<int>& pressure)
    {
        std::vector<std::string> names;
        std::vector<std::set<int>> fields;
        generateFieldSplitSubdomains(names, fields);
        TBOX_ASSERT(names.size() == 2 && fields.size() == 2);
        for (std::size_t i = 0; i < names.size(); ++i)
        {
            TBOX_ASSERT(names[i] == "velocity" || names[i] == "pressure");
            TBOX_ASSERT(fields[i] == (names[i] == "velocity" ? velocity : pressure));
        }
        plog << "velocity_dofs = " << velocity.size() << " pressure_dofs = " << pressure.size() << '\n';
    }
};

// Algebraic cases use a known solution and deliberately interleaved Stokes fields.
int
check_eigen_local(Pointer<Database> test)
{
    const std::string mode = test->getString("eigen_local");
    const bool rank_test = mode == "rank";
    const bool nonsymmetric = mode == "nonsymmetric";
    const bool basic = mode == "basic";
    const std::vector<std::string> types =
        basic ?
            std::vector<std::string>{
                "FULL_PIV_LU", "COL_PIV_HOUSEHOLDER_QR", "FULL_PIV_HOUSEHOLDER_QR", "cod", "JACOBI_SVD", "BDC_SVD"
            } :
        mode == "all" ?
            std::vector<std::string>{ "LLT",
                                      "ldlt",
                                      "PartialPivLU",
                                      "full-piv-lu",
                                      "HOUSEHOLDER_QR",
                                      "COL_PIV_HOUSEHOLDER_QR",
                                      "cod",
                                      "FULL_PIV_HOUSEHOLDER_QR",
                                      "JACOBI_SVD",
                                      "BDC_SVD" } :
        rank_test ?
            std::vector<std::string>{ "cod", "JACOBI_SVD", "BDC_SVD", "FULL_PIV_LU", "FULL_PIV_HOUSEHOLDER_QR" } :
            std::vector<std::string>{ "FULL_PIV_HOUSEHOLDER_QR" };
    const std::vector<std::string> backends =
        nonsymmetric ? std::vector<std::string>{ "eigen", "eigen-pseudoinverse", "eigen-schur-complement" } :
        mode == "all" || rank_test || basic ?
                       std::vector<std::string>{ "eigen", "eigen-pseudoinverse" } :
                       std::vector<std::string>{ test->getStringWithDefault("backend", "eigen-schur-complement") };
    Mat mat = nullptr;
    Vec rhs = nullptr, result = nullptr;
    int ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 3, 3, 3, nullptr, 3, nullptr, &mat);
    IBTK_CHKERRQ(ierr);
    PetscInt lo = 0, hi = 0;
    ierr = MatGetOwnershipRange(mat, &lo, &hi);
    IBTK_CHKERRQ(ierr);
    const double entries[3][3] = { { basic ? 1.0 : 4.0, basic || nonsymmetric ? 2.0 : 1.0, basic ? 0.0 : 1.0 },
                                   { basic ? 0.0 : 1.0, basic ? 0.0 : 3.0, basic ? 0.0 : -1.0 },
                                   { basic || nonsymmetric ? 0.0 : 1.0, basic ? 0.0 : -1.0, basic ? 2.0 : 5.0 } };
    const double target[3] = { 1, -2, 3 };
    for (PetscInt i = lo; i < hi; ++i)
    {
        for (PetscInt j = 0; j < 3; ++j)
        {
            const double value = rank_test ? (i == j ? (i == 0 ? 2.0 : i == 1 ? 0.01 : 0.0) : 0.0) : entries[i][j];
            ierr = MatSetValue(mat, i, j, value, INSERT_VALUES);
            IBTK_CHKERRQ(ierr);
        }
    }
    ierr = MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    IBTK_CHKERRQ(ierr);
    ierr = MatCreateVecs(mat, &rhs, &result);
    IBTK_CHKERRQ(ierr);
    for (PetscInt i = lo; i < hi; ++i)
    {
        double value = basic ? (i == 0 ? 2.0 : i == 1 ? 0.0 : 4.0) : rank_test ? i + 2.0 : 0.0;
        if (!rank_test && !basic)
        {
            for (PetscInt j = 0; j < 3; ++j)
            {
                value += entries[i][j] * target[j];
            }
        }
        ierr = VecSetValue(rhs, i, value, INSERT_VALUES);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecAssemblyBegin(rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(rhs);
    IBTK_CHKERRQ(ierr);
    int failures = 0;
    for (const std::string& backend_name : backends)
    {
        for (const std::string& type : types)
        {
            Pointer<MemoryDatabase> db = new MemoryDatabase("eigen");
            db->putString("eigen_subdomain_solver_type", type);
            db->putString("eigen_subdomain_pseudoinverse_type", type);
            db->putDouble("eigen_subdomain_solver_threshold", rank_test ? 0.1 : -1.0);
            db->putDouble("eigen_subdomain_pseudoinverse_threshold", rank_test ? 0.1 : -1.0);
            std::unique_ptr<PETScLevelSolverShellBackend> backend =
                PETScLevelSolverShellBackendManager::get_manager().allocateBackend(backend_name, db);
            for (int cycle = 0; cycle < 2; ++cycle)
            {
                std::vector<IS> overlap(2), partition(2);
                ierr = ISCreateStride(PETSC_COMM_SELF, 3, 0, 1, &overlap[0]);
                IBTK_CHKERRQ(ierr);
                ierr = ISCreateStride(PETSC_COMM_SELF, 0, 0, 1, &overlap[1]);
                IBTK_CHKERRQ(ierr);
                ierr = ISDuplicate(overlap[0], &partition[0]);
                IBTK_CHKERRQ(ierr);
                ierr = ISDuplicate(overlap[1], &partition[1]);
                IBTK_CHKERRQ(ierr);
                StaggeredStokesEigenSchurComplementShellBackend* schur =
                    dynamic_cast<StaggeredStokesEigenSchurComplementShellBackend*>(backend.get());
                if (schur && mode != "missing_fields")
                {
                    std::set<int> velocity = cycle == 0 ? std::set<int>{ 0, 2 } : std::set<int>{ 1 };
                    std::set<int> pressure = cycle == 0 ? std::set<int>{ 1 } : std::set<int>{ 0, 2 };
                    if (mode == "invalid_fields")
                    {
                        pressure.insert(0);
                    }
                    if (mode == "velocity_only")
                    {
                        velocity = { 0, 1, 2 };
                        pressure.clear();
                    }
                    if (mode == "pressure_only")
                    {
                        pressure = { 0, 1, 2 };
                        velocity.clear();
                    }
                    schur->initializeSolverState(mat,
                                                 rhs,
                                                 result,
                                                 overlap,
                                                 partition,
                                                 velocity,
                                                 pressure,
                                                 "test",
                                                 false,
                                                 PETScLevelSolverShellTraversal::FORWARD);
                }
                else
                {
                    backend->initializeSolverState(mat, rhs, result, overlap, partition, "test");
                }
                // The initializer borrows these only during setup. Destroy them before applying.
                for (IS& is : overlap)
                {
                    ierr = ISDestroy(&is);
                    IBTK_CHKERRQ(ierr);
                }
                for (IS& is : partition)
                {
                    ierr = ISDestroy(&is);
                    IBTK_CHKERRQ(ierr);
                }
                if (mode == "missing_fields" || mode == "invalid_fields" || mode == "parallel")
                {
                    return 0;
                }
                for (int application = 0; application < 3; ++application)
                {
                    backend->apply(rhs, result);
                    const PetscScalar* values = nullptr;
                    ierr = VecGetArrayRead(result, &values);
                    IBTK_CHKERRQ(ierr);
                    double error = 0.0;
                    for (PetscInt j = 0; j < 3; ++j)
                    {
                        const bool minimum_norm = type == "cod" || type == "JACOBI_SVD" || type == "BDC_SVD";
                        const double expected = basic     ? (j == 2       ? 2.0 :
                                                             minimum_norm ? (j == 0 ? 0.4 : 0.8) :
                                                             j == 0       ? 0.0 :
                                                                            1.0) :
                                                rank_test ? (j == 0 ? 1.0 : 0.0) :
                                                            target[j];
                        if (!std::isfinite(PetscRealPart(values[j])))
                        {
                            ++failures;
                        }
                        error = std::max(error, std::abs(PetscRealPart(values[j]) - expected));
                    }
                    if (!std::isfinite(error) || error > 1.e-10)
                    {
                        ++failures;
                    }
                    if (cycle == 0 && application == 0)
                    {
                        plog << backend_name << " " << type << " solution = " << values[0] << " " << values[1] << " "
                             << values[2] << '\n';
                    }
                    ierr = VecRestoreArrayRead(result, &values);
                    IBTK_CHKERRQ(ierr);
                }
                backend->deallocateSolverState();
                backend->deallocateSolverState();
                // Reinitialization must rebuild factors for changed matrix values.
                ierr = MatScale(mat, 2.0);
                IBTK_CHKERRQ(ierr);
                ierr = VecScale(rhs, 2.0);
                IBTK_CHKERRQ(ierr);
            }
            ierr = MatScale(mat, 0.25);
            IBTK_CHKERRQ(ierr);
            ierr = VecScale(rhs, 0.25);
            IBTK_CHKERRQ(ierr);
        }
    }
    ierr = MatDestroy(&mat);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&result);
    IBTK_CHKERRQ(ierr);
    plog << "failures = " << failures << '\n';
    return failures ? 1 : 0;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTKInit init(argc, argv, MPI_COMM_WORLD);
    Pointer<Logger::Appender> appender = new TestAppender();
    Logger::getInstance()->setAbortAppender(appender);
    Pointer<AppInitializer> app = new AppInitializer(argc, argv, "output");
    Pointer<Database> input = app->getInputDatabase();
    Pointer<Database> test = input->getDatabase("test");
    if (test->keyExists("cav_construction"))
    {
        return check_cav_construction(app, test);
    }
    if (test->keyExists("ca_construction"))
    {
        return check_ca_construction(app, test);
    }
    if (test->keyExists("eigen_local"))
    {
        plog << std::setprecision(12);
        return check_eigen_local(test);
    }
    const bool cav = test->keyExists("cav_application");
    const std::string cav_scenario = cav ? test->getString("cav_application") : "";
    const bool cav_override = cav_scenario == "pc_override";
    const bool cav_strict = test->getStringWithDefault("coupling_aware_asm_closure_policy", "RELAXED") == "STRICT";
    const bool boundary = test->getBoolWithDefault("boundary", false);
    const bool lifetime = (cav && !cav_override) || test->getBoolWithDefault("lifetime", false);
    const bool invalid = test->getBoolWithDefault("invalid", false);
    const std::string shell_type = test->getStringWithDefault("shell_pc_type", "multiplicative");
    const bool multiplicative = shell_type.rfind("multiplicative", 0) == 0;
    std::string traversal_name = test->getStringWithDefault("shell_pc_subdomain_traversal", "FORWARD");
    std::transform(traversal_name.begin(),
                   traversal_name.end(),
                   traversal_name.begin(),
                   [](const unsigned char c) { return static_cast<char>(std::toupper(c)); });
    const PETScLevelSolverShellTraversal traversal =
        traversal_name == "REVERSE"   ? PETScLevelSolverShellTraversal::REVERSE :
        traversal_name == "SYMMETRIC" ? PETScLevelSolverShellTraversal::SYMMETRIC :
                                        PETScLevelSolverShellTraversal::FORWARD;
    if (test->getBoolWithDefault("traversal_stages", false))
    {
        int failures = 0;
        for (PETScLevelSolverShellTraversal order : { PETScLevelSolverShellTraversal::FORWARD,
                                                      PETScLevelSolverShellTraversal::REVERSE,
                                                      PETScLevelSolverShellTraversal::SYMMETRIC })
        {
            for (const bool fallback : { false, true })
            {
                failures += check_stages(fallback, false, order);
                failures += check_stages(fallback, false, order, 0);
                failures += check_stages(fallback, false, order, 1);
            }
        }
        plog << "failures = " << failures << '\n';
        return failures ? 1 : 0;
    }
    if (test->getBoolWithDefault("stages", false))
    {
        return check_stages(test->getBoolWithDefault("fallback", false), test->getBoolWithDefault("invalidate", false));
    }
    const int hand_failures = test->getBoolWithDefault("hand_solve", false) ?
                                  check_hand_solve(shell_type == "multiplicative-blas-lapack",
                                                   traversal,
                                                   test->keyExists("shell_pc_subdomain_traversal"),
                                                   test->getStringWithDefault("hand_backend", ""),
                                                   !multiplicative) :
                                  0;
    const bool diagonal_operator = test->getBoolWithDefault("diagonal_operator", false);
    const double small_diagonal = test->getDoubleWithDefault("small_diagonal", 1.0e-12);
    const bool all_blas_modes = test->getBoolWithDefault("all_blas_modes", false);
    const std::vector<std::string> solver_types =
        all_blas_modes ?
            std::vector<std::string>{ "", "svd", "lu", "symmetric-indefinite", "qr" } :
            std::vector<std::string>{ test->getStringWithDefault("blas_lapack_subdomain_solver_type", "") };
    const auto hierarchy_data = setup_hierarchy<NDIM>(app);
    Pointer<PatchHierarchy<NDIM>> hierarchy = std::get<0>(hierarchy_data);
    Pointer<PatchLevel<NDIM>> level = hierarchy->getPatchLevel(0);
    VariableDatabase<NDIM>* variables = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> context = variables->getContext("shell_test");
    Pointer<SideVariable<NDIM, double>> u = new SideVariable<NDIM, double>("u");
    Pointer<CellVariable<NDIM, double>> p = new CellVariable<NDIM, double>("p");
    Pointer<SideVariable<NDIM, double>> f = new SideVariable<NDIM, double>("f");
    Pointer<CellVariable<NDIM, double>> h = new CellVariable<NDIM, double>("h");
    Pointer<SideVariable<NDIM, int>> u_dof = new SideVariable<NDIM, int>("u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof = new CellVariable<NDIM, int>("p_dof");
    const int ui = variables->registerVariableAndContext(u, context, IntVector<NDIM>(1));
    const int pi = variables->registerVariableAndContext(p, context, IntVector<NDIM>(1));
    const int fi = variables->registerVariableAndContext(f, context, IntVector<NDIM>(1));
    const int hi = variables->registerVariableAndContext(h, context, IntVector<NDIM>(1));
    const int udi = variables->registerVariableAndContext(u_dof, context, IntVector<NDIM>(1));
    const int pdi = variables->registerVariableAndContext(p_dof, context, IntVector<NDIM>(1));
    for (int index : { ui, pi, fi, hi, udi, pdi })
    {
        level->allocatePatchData(index);
    }
    SAMRAIVectorReal<NDIM, double> x("x", hierarchy, 0, 0), b("b", hierarchy, 0, 0);
    x.addComponent(u, ui);
    x.addComponent(p, pi);
    b.addComponent(f, fi);
    b.addComponent(h, hi);
    x.setToScalar(0.0);
    b.setToScalar(0.0);
    const double wavenumber = 2.0 * std::acos(-1.0) / input->getInteger("N");
    for (PatchLevel<NDIM>::Iterator patch_number(level); patch_number; patch_number++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(patch_number());
        Pointer<SideData<NDIM, double>> force = patch->getPatchData(fi);
        for (int axis = 0; axis < NDIM; ++axis)
        {
            const Box<NDIM> box = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
            for (Box<NDIM>::Iterator it(box); it; it++)
            {
                (*force)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)) =
                    boundary ? axis + 1.0 : std::sin(wavenumber * it()(axis)) + 0.25 * (axis + 1.0);
            }
        }
    }
    CAFields cav_fields;
    std::vector<int> dofs;
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(dofs, udi, pdi, level);
    if (cav)
    {
        for (PatchLevel<NDIM>::Iterator patch_number(level); patch_number; patch_number++)
        {
            Pointer<Patch<NDIM>> patch = level->getPatch(patch_number());
            Pointer<SideData<NDIM, int>> velocity = patch->getPatchData(udi);
            Pointer<CellData<NDIM, int>> pressure = patch->getPatchData(pdi);
            Pointer<CellData<NDIM, double>> pressure_rhs = patch->getPatchData(hi);
            for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
            {
                CACell cell{};
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    cell[axis] = it()(axis);
                }
                for (int axis = 0; axis < NDIM; ++axis)
                {
                    cav_fields[cell][axis] = (*velocity)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower));
                }
                cav_fields[cell][NDIM] = (*pressure)(it());
                (*pressure_rhs)(it()) = 0.125 * std::sin(wavenumber * it()(0));
            }
        }
    }
    Vec rhs = nullptr, expected = nullptr, actual = nullptr;
    int ierr = VecCreateMPI(PETSC_COMM_WORLD, dofs[IBTK_MPI::getRank()], PETSC_DETERMINE, &rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(rhs, &expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecDuplicate(rhs, &actual);
    IBTK_CHKERRQ(ierr);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(rhs, fi, udi, hi, pdi, level);
    Pointer<MemoryDatabase> db = new MemoryDatabase("solver");
    db->putString("ksp_type", "preonly");
    db->putString("pc_type", test->getStringWithDefault("pc_type", "shell"));
    if (test->keyExists("shell_pc_type"))
    {
        db->putString("shell_pc_type", shell_type);
    }
    for (const std::string key : { "eigen_subdomain_solver_type",
                                   "eigen_subdomain_pseudoinverse_type",
                                   "a00_solver_type",
                                   "schur_solver_type" })
    {
        if (test->keyExists(key))
        {
            db->putString(key, test->getString(key));
        }
    }
    db->putBool("initial_guess_nonzero", false);
    db->putInteger("max_iterations", 1);
    if (cav)
    {
        db->putString("asm_subdomain_construction_mode", "COUPLING_AWARE");
        db->putString("coupling_aware_asm_patch_seed_type", "PRESSURE_CELL");
        db->putString("coupling_aware_asm_closure_policy", cav_strict ? "STRICT" : "RELAXED");
    }
    int box_size[NDIM];
    std::fill_n(box_size, NDIM, test->getIntegerWithDefault("box_size", boundary ? input->getInteger("N") : 4));
    db->putIntegerArray("subdomain_box_size", box_size, NDIM);
    int failures = hand_failures;
    plog << std::setprecision(12);
    for (const std::string& solver_type : solver_types)
    {
        if (!solver_type.empty())
        {
            db->putString("blas_lapack_subdomain_solver_type", solver_type);
        }
        if (test->keyExists("blas_lapack_subdomain_solver_rcond"))
        {
            db->putDouble("blas_lapack_subdomain_solver_rcond", test->getDouble("blas_lapack_subdomain_solver_rcond"));
        }
        FieldCheckingSolver solver("shell_solver", db, "shell_");
        PoissonSpecifications coefficients("coefficients");
        coefficients.setCConstant(1.0);
        coefficients.setDConstant(-1.0);
        solver.setVelocityPoissonSpecifications(coefficients);
        solver.setComponentsHaveNullSpace(false, !diagonal_operator);
        std::vector<std::unique_ptr<LocationIndexRobinBcCoefs<NDIM>>> bc_storage;
        std::vector<RobinBcCoefStrategy<NDIM>*> bcs(NDIM, nullptr);
        if (boundary)
        {
            for (int axis = 0; axis < NDIM; ++axis)
            {
                bc_storage.push_back(std::make_unique<LocationIndexRobinBcCoefs<NDIM>>("bc", nullptr));
                for (int face = 0; face < 2 * NDIM; ++face)
                {
                    bc_storage.back()->setBoundaryValue(face, axis + 1.0);
                }
                bcs[axis] = bc_storage.back().get();
            }
            Pointer<StaggeredStokesPhysicalBoundaryHelper> helper = new StaggeredStokesPhysicalBoundaryHelper();
            helper->cacheBcCoefData(bcs, 0.0, hierarchy);
            solver.setPhysicalBcCoefs(bcs, nullptr);
            solver.setPhysicalBoundaryHelper(helper);
            solver.setHomogeneousBc(false);
        }
        Mat elasticity = nullptr;
        std::vector<std::set<int>> cav_expected, previous_patches;
        if (cav && !cav_override)
        {
            ierr = MatCreateSeqDense(PETSC_COMM_WORLD, dofs.front(), dofs.front(), nullptr, &elasticity);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyBegin(elasticity, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(elasticity, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            cav_expected = cav_application_patches(cav_fields, input->getInteger("N"), cav_strict, 0, elasticity);
            if (cav_scenario != "missing_matrix")
            {
                PetscInt before = 0, after = 0;
                ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(elasticity), &before);
                IBTK_CHKERRQ(ierr);
                solver.setCouplingAwareASMConstructionMat(elasticity);
                ierr = PetscObjectGetReference(reinterpret_cast<PetscObject>(elasticity), &after);
                IBTK_CHKERRQ(ierr);
                if (before != after)
                {
                    ++failures;
                }
            }
            solver.setAugmentedOperatorMat(elasticity);
        }
        Mat supplied = nullptr;
        if (diagonal_operator)
        {
            // A diagonal operator has an analytic truncated pseudoinverse: entries
            // below the requested cutoff contribute zero, and the others divide by two.
            PetscInt n = 0;
            ierr = VecGetSize(rhs, &n);
            IBTK_CHKERRQ(ierr);
            ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 1, nullptr, &supplied);
            IBTK_CHKERRQ(ierr);
            for (PetscInt j = 0; j < n; ++j)
            {
                ierr = MatSetValue(supplied, j, j, j % 2 == 0 ? 2.0 : small_diagonal, INSERT_VALUES);
                IBTK_CHKERRQ(ierr);
            }
            ierr = MatAssemblyBegin(supplied, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(supplied, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            solver.setOperatorMat(supplied);
        }
        Vec default_forward = nullptr;
        if (test->getBoolWithDefault("compare_forward", false))
        {
            ierr = VecDuplicate(rhs, &default_forward);
            IBTK_CHKERRQ(ierr);
            solver.initializeSolverState(x, b);
            if (!solver.solveSystem(x, b))
            {
                ++failures;
            }
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(default_forward, ui, udi, pi, pdi, level);
            solver.deallocateSolverState();
        }
        if (test->keyExists("shell_pc_subdomain_traversal"))
        {
            db->putString("shell_pc_subdomain_traversal", test->getString("shell_pc_subdomain_traversal"));
        }
        solver.initializeSolverState(x, b);
        if (cav_scenario == "initialized_setter")
        {
            solver.setCouplingAwareASMConstructionMat(elasticity);
        }
        if (cav && !cav_override && cav_scenario != "apply")
        {
            solver.deallocateSolverState();
            ierr = MatDestroy(&elasticity);
            IBTK_CHKERRQ(ierr);
            return 0;
        }
        if (invalid)
        {
            solver.deallocateSolverState();
            return 0;
        }
        if (lifetime && !diagonal_operator && !cav)
        {
            Mat assembled = nullptr;
            ierr = KSPGetOperators(solver.getPETScKSP(), &assembled, nullptr);
            IBTK_CHKERRQ(ierr);
            ierr = MatDuplicate(assembled, MAT_COPY_VALUES, &supplied);
            IBTK_CHKERRQ(ierr);
            solver.deallocateSolverState();
            solver.setOperatorMat(supplied);
            solver.setOperatorMat(supplied);
            solver.initializeSolverState(x, b);
        }
        if (all_blas_modes)
        {
            plog << "solver_type = " << (solver_type.empty() ? "default" : solver_type) << '\n';
        }
        for (int cycle = 0; cycle < (lifetime ? 2 : 1); ++cycle)
        {
            if (shell_type.find("eigen-schur-complement") != std::string::npos)
            {
                std::set<int> velocity, pressure;
                for (PatchLevel<NDIM>::Iterator patch_number(level); patch_number; patch_number++)
                {
                    Pointer<Patch<NDIM>> patch = level->getPatch(patch_number());
                    Pointer<SideData<NDIM, int>> u_indices = patch->getPatchData(udi);
                    Pointer<CellData<NDIM, int>> p_indices = patch->getPatchData(pdi);
                    for (Box<NDIM>::Iterator it(patch->getBox()); it; it++)
                    {
                        pressure.insert((*p_indices)(it()));
                    }
                    for (int axis = 0; axis < NDIM; ++axis)
                    {
                        const Box<NDIM> box = SideGeometry<NDIM>::toSideBox(patch->getBox(), axis);
                        for (Box<NDIM>::Iterator it(box); it; it++)
                        {
                            velocity.insert((*u_indices)(SideIndex<NDIM>(it(), axis, SideIndex<NDIM>::Lower)));
                        }
                    }
                }
                solver.checkFields(velocity, pressure);
            }
            Mat mat = nullptr;
            PC pc = nullptr;
            ierr = KSPGetOperators(solver.getPETScKSP(), &mat, nullptr);
            IBTK_CHKERRQ(ierr);
            ierr = KSPGetPC(solver.getPETScKSP(), &pc);
            IBTK_CHKERRQ(ierr);
            PCType pc_type = nullptr;
            ierr = PCGetType(pc, &pc_type);
            IBTK_CHKERRQ(ierr);
            if (std::string(pc_type) != (cav_override ? "none" : "shell") || (lifetime && !cav && mat != supplied))
            {
                ++failures;
            }
            if (diagonal_operator)
            {
                const PetscScalar* rhs_values = nullptr;
                PetscScalar* expected_values = nullptr;
                PetscInt n = 0;
                ierr = VecGetSize(rhs, &n);
                IBTK_CHKERRQ(ierr);
                ierr = VecGetArrayRead(rhs, &rhs_values);
                IBTK_CHKERRQ(ierr);
                ierr = VecGetArray(expected, &expected_values);
                IBTK_CHKERRQ(ierr);
                for (PetscInt j = 0; j < n; ++j)
                {
                    expected_values[j] = j % 2 == 0 ? rhs_values[j] / (cycle == 0 ? 2.0 : 4.0) : 0.0;
                }
                ierr = VecRestoreArray(expected, &expected_values);
                IBTK_CHKERRQ(ierr);
                ierr = VecRestoreArrayRead(rhs, &rhs_values);
                IBTK_CHKERRQ(ierr);
            }
            else if (!boundary)
            {
                std::vector<IS>* overlap = nullptr;
                std::vector<IS>* partition = nullptr;
                solver.getASMSubdomains(&partition, &overlap);
                if (cav_override)
                {
                    if (!partition->empty() || !overlap->empty())
                    {
                        ++failures;
                    }
                    plog << "effective_pc = " << pc_type << '\n';
                }
                else if (cav)
                {
                    const std::vector<std::set<int>> patches = ca_read_sets(*overlap);
                    if (!partition->empty() || patches != cav_expected || (cycle == 1 && patches == previous_patches))
                    {
                        ++failures;
                    }
                    previous_patches = patches;
                    plog << "pressure_patches = " << patches.size() << "\npartition_size = " << partition->size()
                         << '\n';
                }
                if (cav_override)
                {
                    ierr = VecCopy(rhs, expected);
                    IBTK_CHKERRQ(ierr);
                }
                else
                {
                    reference_action(mat, rhs, expected, *overlap, *partition, multiplicative, traversal);
                }
                // Left-preconditioned PETSc KSP removes the operator nullspace after PCApply.
                MatNullSpace nullspace = nullptr;
                ierr = MatGetNullSpace(mat, &nullspace);
                IBTK_CHKERRQ(ierr);
                if (nullspace)
                {
                    ierr = MatNullSpaceRemove(nullspace, expected);
                    IBTK_CHKERRQ(ierr);
                }
            }
            x.setToScalar(0.0);
            if (!solver.solveSystem(x, b))
            {
                ++failures;
            }
            StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(actual, ui, udi, pi, pdi, level);
            const double action_norm = norm_inf(actual);
            if (cav)
            {
                PetscScalar pressure_sum = 0.0;
                const PetscScalar* values = nullptr;
                ierr = VecGetArrayRead(actual, &values);
                IBTK_CHKERRQ(ierr);
                for (const auto& cell : cav_fields)
                {
                    pressure_sum += values[cell.second[NDIM]];
                }
                ierr = VecRestoreArrayRead(actual, &values);
                IBTK_CHKERRQ(ierr);
                const double pressure_mean = static_cast<double>(PetscRealPart(pressure_sum)) / cav_fields.size();
                if (!std::isfinite(pressure_mean) || std::abs(pressure_mean) > 1.0e-9)
                {
                    ++failures;
                }
                plog << "pressure_mean = " << pressure_mean << '\n';
            }
            if (default_forward)
            {
                ierr = VecAXPY(default_forward, -1.0, actual);
                IBTK_CHKERRQ(ierr);
                const double default_error = norm_inf(default_forward);
                if (!std::isfinite(default_error) || default_error > 1.0e-12)
                {
                    ++failures;
                }
                plog << "default_forward_error = " << default_error << '\n';
                ierr = VecCopy(actual, default_forward);
                IBTK_CHKERRQ(ierr);
            }
            if (boundary)
            {
                // Constant velocity (1,2), pressure zero solves -Laplace(u)+u+grad(p)=u.
                // Its Dirichlet values are prescribed on every face, independently of RHS packing.
                StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(expected, fi, udi, hi, pdi, level);
            }
            ierr = VecAXPY(actual, -1.0, expected);
            IBTK_CHKERRQ(ierr);
            const double error = norm_inf(actual);
            if (!std::isfinite(error) || error > 1.0e-9 || action_norm <= 0.0)
            {
                ++failures;
            }
            plog << "action_norm = " << action_norm << "\nerror = " << error << '\n';
            if (lifetime)
            {
                if (!solver.solveSystem(x, b))
                {
                    ++failures;
                }
                StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(actual, ui, udi, pi, pdi, level);
                ierr = VecAXPY(actual, -1.0, expected);
                IBTK_CHKERRQ(ierr);
                const double repeated_error = norm_inf(actual);
                if (!std::isfinite(repeated_error) || repeated_error > 1.0e-9)
                {
                    ++failures;
                }
                plog << "repeated_error = " << repeated_error << '\n';
            }
            std::vector<IS>* cav_overlap = nullptr;
            std::vector<IS>* cav_partition = nullptr;
            if (cav)
            {
                // Query while initialized; the returned containers outlive solver state.
                solver.getASMSubdomains(&cav_partition, &cav_overlap);
            }
            solver.deallocateSolverState();
            if (cav)
            {
                if (!cav_overlap->empty() || !cav_partition->empty())
                {
                    ++failures;
                }
                if (!cav_override)
                {
                    PetscReal matrix_norm = 0.0;
                    ierr = MatNorm(elasticity, NORM_INFINITY, &matrix_norm);
                    IBTK_CHKERRQ(ierr);
                    if (matrix_norm != 0.5 * (cycle == 0 ? 1 : 2 * NDIM))
                    {
                        ++failures;
                    }
                }
                if (cycle == 0 && !cav_override)
                {
                    solver.setAugmentedOperatorMat(nullptr);
                    cav_expected =
                        cav_application_patches(cav_fields, input->getInteger("N"), cav_strict, 1, elasticity);
                    solver.setCouplingAwareASMConstructionMat(elasticity);
                    solver.setAugmentedOperatorMat(elasticity);
                    solver.initializeSolverState(x, b);
                }
            }
            else if (lifetime)
            {
                PetscReal matrix_norm = 0.0;
                ierr = MatNorm(supplied, NORM_INFINITY, &matrix_norm);
                IBTK_CHKERRQ(ierr);
                if (!(matrix_norm > 0.0))
                {
                    ++failures;
                }
                if (cycle == 0)
                {
                    if (shell_type == "additive-blas-lapack")
                    {
                        ierr = MatScale(supplied, 2.0);
                        IBTK_CHKERRQ(ierr);
                    }
                    solver.initializeSolverState(x, b);
                }
            }
        }
        ierr = MatDestroy(&elasticity);
        IBTK_CHKERRQ(ierr);
        ierr = VecDestroy(&default_forward);
        IBTK_CHKERRQ(ierr);
        ierr = MatDestroy(&supplied);
        IBTK_CHKERRQ(ierr);
    }
    ierr = VecDestroy(&rhs);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&expected);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&actual);
    IBTK_CHKERRQ(ierr);
    plog << "failures = " << failures << '\n';
    return failures ? 1 : 0;
}
