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

#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/MemoryDatabase.h>

#include <CellData.h>
#include <CellVariable.h>
#include <PoissonSpecifications.h>
#include <SideData.h>
#include <SideGeometry.h>
#include <SideVariable.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <set>

#include "../tests.h"
#include "coupling_aware_asm_test_utilities.h"

#include <ibamr/app_namespaces.h>

namespace
{
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
    const std::string scenario = test->getString("ca_construction");
    if (scenario == "parallel")
    {
        // The construction aborts on more than one rank before it reads the matrix.
        std::vector<std::set<int>> overlap, partition;
        StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains(
            overlap, partition, counts, udi, pdi, level, nullptr);
        return 0;
    }
    TBOX_ASSERT(counts.size() == 1);
    const int total = counts.front();
    plog << "dofs = " << total << '\n';
    const int n = app->getInputDatabase()->getInteger("N");
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
            TBOX_ERROR("Failed check: overlap.size() != partition.size().\n");
        }
        for (std::size_t k = 0; k < partition.size(); ++k)
        {
            for (const int dof : partition[k])
            {
                if (!overlap[k].count(dof) || !owned.insert(dof).second)
                {
                    TBOX_ERROR("Failed check: !overlap[k].count(dof) || !owned.insert(dof).second.\n");
                }
            }
        }
        if (owned.size() != static_cast<std::size_t>(total))
        {
            TBOX_ERROR("Failed check: owned.size() != static_cast<std::size_t>(total).\n");
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
                TBOX_ERROR("Failed check: overlap != expected.\n");
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
                    TBOX_ERROR("Failed check: overlap != expected.\n");
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
            TBOX_ERROR("Failed check: overlap != expected.\n");
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
                    TBOX_ERROR("Failed check: overlap != expected.\n");
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
                    TBOX_ERROR("Failed check: canonical != geometric || canonical != expected.\n");
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
                // The type must be checked both at construction and after a command-line override.
                db->putString("pc_type", scenario == "pc_type_ilu" ? "ilu" : pc_type);
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
                    if (scenario == "pc_type_option")
                    {
                        ierr = PetscOptionsSetValue(nullptr, "-ca_pc_type", "ilu");
                        IBTK_CHKERRQ(ierr);
                    }
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
                        TBOX_ERROR("CA lifetime mismatch: " << pc_type << " " << policy << " cycle=" << cycle
                                                            << " overlap_matches=" << overlap_matches
                                                            << " partition_matches=" << partition_matches << '\n');
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
    return 0;
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
                    TBOX_ERROR("Failed check: patches != expected || seeds != ordered_seeds.\n");
                }
            }
        }
    }
    std::vector<std::set<int>> velocity_patches, partition;
    StaggeredStokesPETScMatUtilities::construct_patch_level_coupling_aware_asm_subdomains(
        velocity_patches, partition, counts, udi, pdi, level, elasticity);
    // Every velocity-seeded patch closes at least one cell, and coupling through the elasticity enlarges some of them
    // beyond the standard patch of a cell.
    bool enlarged = false;
    for (const std::set<int>& patch : velocity_patches)
    {
        const bool closes_a_cell =
            std::any_of(standard.begin(),
                        standard.end(),
                        [&](const std::set<int>& cell_patch)
                        { return std::includes(patch.begin(), patch.end(), cell_patch.begin(), cell_patch.end()); });
        if (!closes_a_cell)
        {
            TBOX_ERROR("A velocity-seeded patch of " << patch.size() << " DOFs contains no standard patch.\n");
        }
        enlarged = enlarged || patch.size() > 2 * NDIM + 1;
    }
    if (!enlarged)
    {
        TBOX_ERROR("No velocity-seeded patch exceeds the " << 2 * NDIM + 1 << " DOFs of a standard patch.\n");
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
                TBOX_ERROR(
                    "Failed check: seeds != expected_seeds || patches != (policy == "
                    "CouplingAwareASMClosurePolicy::RELAXED ? expected_relaxed : standard).\n");
            }
        }
    }
    if (expected_relaxed == standard || !expected_relaxed.front().count(remote_velocity) ||
        standard.front().count(remote_velocity))
    {
        TBOX_ERROR(
            "Failed check: expected_relaxed == standard || !expected_relaxed.front().count(remote_velocity) || "
            "standard.front().count(remote_velocity).\n");
    }
    // Changing the same borrowed matrix must remove the previous expansion.
    ierr = MatZeroEntries(elasticity);
    IBTK_CHKERRQ(ierr);
    assemble();
    StaggeredStokesPETScMatUtilities::construct_patch_level_pressure_cell_seeded_cav_patches(
        patches, seeds, counts, udi, pdi, level, elasticity);
    if (patches != standard || seeds != expected_seeds)
    {
        TBOX_ERROR("Failed check: patches != standard || seeds != expected_seeds.\n");
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
        TBOX_ERROR(
            "Failed check: patches.size() != standard.size() || seeds != expected_seeds || patches.front() != "
            "complete.\n");
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
        TBOX_ERROR(
            "Failed check: patches.size() != standard.size() || seeds != expected_seeds || patches.front() != "
            "complete_relaxed.\n");
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
            const double entry_error = static_cast<double>(PetscAbsScalar(values[k] - expected[columns[k]]));
            // std::max would discard a NaN entry error.
            pressure_row_error = std::isfinite(entry_error) ? std::max(pressure_row_error, entry_error) :
                                                              std::numeric_limits<double>::infinity();
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
        TBOX_ERROR("Failed check: !valid.\n");
    }
    plog << "pressure_shift_error = " << shift_error << '\n';
    ierr = MatNullSpaceDestroy(&nullspace);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&pressure_shift);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&action);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&mac);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&elasticity);
    IBTK_CHKERRQ(ierr);
    level->deallocatePatchData(udi);
    level->deallocatePatchData(pdi);
    plog << "pressure_patches = " << standard.size() << '\n';
    return 0;
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
    if (test->keyExists("ca_construction"))
    {
        return check_ca_construction(app, test);
    }
    return check_cav_construction(app, test);
}
