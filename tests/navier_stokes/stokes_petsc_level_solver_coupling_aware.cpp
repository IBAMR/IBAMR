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
        // Exercise ordinary initialization, lazy generation, repeated construction,
        // and deallocation. Subdomain inspection uses the existing PETSc IS query.
        for (const std::string policy : { "RELAXED", "STRICT" })
        {
            Pointer<MemoryDatabase> db = new MemoryDatabase("ca_solver");
            // The type must be checked both at construction and after a command-line override.
            db->putString("pc_type", scenario == "pc_type_ilu" ? "ilu" : "asm");
            db->putString("asm_subdomain_construction_mode", "COUPLING_AWARE");
            db->putString("coupling_aware_asm_closure_policy", policy);
            StaggeredStokesPETScLevelSolver solver("ca_solver_" + policy, db, "ca_");
            solver.setOperatorMat(matrix);
            if (scenario == "pc_type_option")
            {
                ierr = PetscOptionsSetValue(nullptr, "-ca_pc_type", "ilu");
                IBTK_CHKERRQ(ierr);
            }
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
            for (int cycle = 0; cycle < 2; ++cycle)
            {
                solver.initializeSolverState(x, b);
                std::vector<IS>* actual_overlap = nullptr;
                std::vector<IS>* actual_partition = nullptr;
                solver.getASMSubdomains(&actual_partition, &actual_overlap);
                if (ca_read_sets(*actual_overlap) != expected || ca_read_sets(*actual_partition) != expected_partition)
                {
                    TBOX_ERROR(
                        "Failed check: ca_read_sets(*actual_overlap) != expected || ca_read_sets(*actual_partition) != "
                        "expected_partition.\n");
                }
                solver.deallocateSolverState();
            }
        }
        level->deallocatePatchData(ui);
        level->deallocatePatchData(pi);
    }
    ierr = MatDestroy(&matrix);
    IBTK_CHKERRQ(ierr);
    level->deallocatePatchData(udi);
    level->deallocatePatchData(pdi);
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
    return 1;
}
