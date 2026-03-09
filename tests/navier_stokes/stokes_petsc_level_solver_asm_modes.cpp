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
#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/MemoryDatabase.h>

#include <petscis.h>

#include <CellVariable.h>
#include <PoissonSpecifications.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "../tests.h"

#include <ibtk/app_namespaces.h>

namespace
{
std::set<int>
is_to_set(IS is)
{
    int n = 0;
    int ierr = ISGetLocalSize(is, &n);
    IBTK_CHKERRQ(ierr);
    const int* idxs = nullptr;
    ierr = ISGetIndices(is, &idxs);
    IBTK_CHKERRQ(ierr);
    std::set<int> dofs;
    dofs.insert(idxs, idxs + n);
    ierr = ISRestoreIndices(is, &idxs);
    IBTK_CHKERRQ(ierr);
    return dofs;
}

void
check_partition_invariants(int& test_failures,
                           const std::vector<std::set<int>>& overlap_sets,
                           const std::vector<std::set<int>>& nonoverlap_sets,
                           const int n_local_dofs)
{
    if (overlap_sets.empty() || nonoverlap_sets.empty() || overlap_sets.size() != nonoverlap_sets.size())
    {
        ++test_failures;
        return;
    }

    std::set<int> overlap_union, nonoverlap_union;
    std::map<int, int> nonoverlap_counts;
    for (std::size_t k = 0; k < overlap_sets.size(); ++k)
    {
        overlap_union.insert(overlap_sets[k].begin(), overlap_sets[k].end());
        nonoverlap_union.insert(nonoverlap_sets[k].begin(), nonoverlap_sets[k].end());
        for (const int dof : nonoverlap_sets[k]) ++nonoverlap_counts[dof];
        if (!std::includes(
                overlap_sets[k].begin(), overlap_sets[k].end(), nonoverlap_sets[k].begin(), nonoverlap_sets[k].end()))
        {
            ++test_failures;
        }
    }

    if (static_cast<int>(overlap_union.size()) != n_local_dofs) ++test_failures;
    if (static_cast<int>(nonoverlap_union.size()) != n_local_dofs) ++test_failures;
    for (const auto& kv : nonoverlap_counts)
    {
        if (kv.second != 1)
        {
            ++test_failures;
            break;
        }
    }
}

struct SubdomainStats
{
    std::size_t n_subdomains = 0;
    std::size_t overlap_min = 0;
    std::size_t overlap_max = 0;
    double overlap_avg = 0.0;
    long long overlap_hash = 0;
    std::size_t nonoverlap_min = 0;
    std::size_t nonoverlap_max = 0;
    double nonoverlap_avg = 0.0;
    long long nonoverlap_hash = 0;
};

SubdomainStats
compute_stats(const std::vector<std::set<int>>& overlap_sets, const std::vector<std::set<int>>& nonoverlap_sets)
{
    SubdomainStats stats;
    stats.n_subdomains = overlap_sets.size();
    if (stats.n_subdomains == 0) return stats;

    std::size_t overlap_sum = 0;
    std::size_t nonoverlap_sum = 0;
    std::size_t overlap_min = std::numeric_limits<std::size_t>::max();
    std::size_t nonoverlap_min = std::numeric_limits<std::size_t>::max();
    std::size_t overlap_max = 0;
    std::size_t nonoverlap_max = 0;
    long long overlap_hash = 0;
    long long nonoverlap_hash = 0;
    for (std::size_t k = 0; k < stats.n_subdomains; ++k)
    {
        const std::size_t n_overlap = overlap_sets[k].size();
        const std::size_t n_nonoverlap = nonoverlap_sets[k].size();
        overlap_sum += n_overlap;
        nonoverlap_sum += n_nonoverlap;
        overlap_min = std::min(overlap_min, n_overlap);
        nonoverlap_min = std::min(nonoverlap_min, n_nonoverlap);
        overlap_max = std::max(overlap_max, n_overlap);
        nonoverlap_max = std::max(nonoverlap_max, n_nonoverlap);
        overlap_hash = overlap_hash * 1315423911LL + static_cast<long long>((k + 1) * (n_overlap + 3));
        nonoverlap_hash = nonoverlap_hash * 1315423911LL + static_cast<long long>((k + 1) * (n_nonoverlap + 7));
    }
    stats.overlap_min = overlap_min;
    stats.overlap_max = overlap_max;
    stats.overlap_avg = static_cast<double>(overlap_sum) / static_cast<double>(stats.n_subdomains);
    stats.overlap_hash = overlap_hash;
    stats.nonoverlap_min = nonoverlap_min;
    stats.nonoverlap_max = nonoverlap_max;
    stats.nonoverlap_avg = static_cast<double>(nonoverlap_sum) / static_cast<double>(stats.n_subdomains);
    stats.nonoverlap_hash = nonoverlap_hash;
    return stats;
}

std::vector<std::set<int>>
collect_is_sets(const std::vector<IS>& is_vec)
{
    std::vector<std::set<int>> out(is_vec.size());
    for (std::size_t k = 0; k < is_vec.size(); ++k) out[k] = is_to_set(is_vec[k]);
    return out;
}

IntVector<NDIM>
read_int_vector_with_default(Pointer<Database> db, const std::string& key, const IntVector<NDIM>& default_value)
{
    IntVector<NDIM> value = default_value;
    if (db->keyExists(key))
    {
        int vals[NDIM];
        db->getIntegerArray(key, vals, NDIM);
        for (int d = 0; d < NDIM; ++d) value(d) = vals[d];
    }
    return value;
}

std::vector<std::vector<int>>
canonicalize_subdomain_sets(const std::vector<std::set<int>>& subdomains)
{
    std::vector<std::vector<int>> canonical;
    canonical.reserve(subdomains.size());
    for (const auto& subdomain : subdomains)
    {
        canonical.emplace_back(subdomain.begin(), subdomain.end());
    }
    std::sort(canonical.begin(), canonical.end());
    return canonical;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();
    Pointer<Database> test_db = input_db->keyExists("test") ? input_db->getDatabase("test") : input_db;
    const std::string asm_mode = test_db->getStringWithDefault("asm_subdomain_construction_mode", "GEOMETRICAL");
    const std::string closure_policy = test_db->getStringWithDefault("coupling_aware_asm_closure_policy", "RELAXED");
    const IntVector<NDIM> box_size = read_int_vector_with_default(test_db, "subdomain_box_size", IntVector<NDIM>(2));
    const IntVector<NDIM> overlap_size =
        read_int_vector_with_default(test_db, "subdomain_overlap_size", IntVector<NDIM>(1));
    const int seed_axis = test_db->getIntegerWithDefault("coupling_aware_asm_seed_axis", 0);
    const int seed_stride = test_db->getIntegerWithDefault("coupling_aware_asm_seed_stride", 1);
    const bool run_comparison = test_db->getBoolWithDefault("compare_against_mode", false);
    const std::string compare_asm_mode =
        test_db->getStringWithDefault("compare_asm_subdomain_construction_mode", "GEOMETRICAL");
    const std::string compare_closure_policy =
        test_db->getStringWithDefault("compare_coupling_aware_asm_closure_policy", "RELAXED");
    const IntVector<NDIM> compare_box_size =
        read_int_vector_with_default(test_db, "compare_subdomain_box_size", box_size);
    const IntVector<NDIM> compare_overlap_size =
        read_int_vector_with_default(test_db, "compare_subdomain_overlap_size", overlap_size);
    const int compare_seed_axis = test_db->getIntegerWithDefault("compare_coupling_aware_asm_seed_axis", seed_axis);
    const int compare_seed_stride =
        test_db->getIntegerWithDefault("compare_coupling_aware_asm_seed_stride", seed_stride);

    int test_failures = 0;
    const auto hierarchy_tuple = setup_hierarchy<NDIM>(app_initializer);
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = std::get<0>(hierarchy_tuple);
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(0);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("stokes_petsc_level_solver_asm_modes_ctx");
    Pointer<SideVariable<NDIM, double>> u_var = new SideVariable<NDIM, double>("solver_u");
    Pointer<CellVariable<NDIM, double>> p_var = new CellVariable<NDIM, double>("solver_p");
    Pointer<SideVariable<NDIM, double>> f_u_var = new SideVariable<NDIM, double>("solver_f_u");
    Pointer<CellVariable<NDIM, double>> f_p_var = new CellVariable<NDIM, double>("solver_f_p");
    const int u_idx = var_db->registerVariableAndContext(u_var, ctx, IntVector<NDIM>(1));
    const int p_idx = var_db->registerVariableAndContext(p_var, ctx, IntVector<NDIM>(1));
    const int f_u_idx = var_db->registerVariableAndContext(f_u_var, ctx, IntVector<NDIM>(1));
    const int f_p_idx = var_db->registerVariableAndContext(f_p_var, ctx, IntVector<NDIM>(1));
    level->allocatePatchData(u_idx);
    level->allocatePatchData(p_idx);
    level->allocatePatchData(f_u_idx);
    level->allocatePatchData(f_p_idx);

    SAMRAIVectorReal<NDIM, double> x_vec("x", patch_hierarchy, 0, 0);
    SAMRAIVectorReal<NDIM, double> b_vec("b", patch_hierarchy, 0, 0);
    x_vec.addComponent(u_var, u_idx);
    x_vec.addComponent(p_var, p_idx);
    b_vec.addComponent(f_u_var, f_u_idx);
    b_vec.addComponent(f_p_var, f_p_idx);

    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("solver_mode_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("solver_mode_p_dof");
    const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, ctx, IntVector<NDIM>(1));
    const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, ctx, IntVector<NDIM>(1));
    level->allocatePatchData(u_dof_index_idx);
    level->allocatePatchData(p_dof_index_idx);
    std::vector<int> num_dofs_per_proc;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);
    const int rank = IBTK_MPI::getRank();
    const int n_local_dofs = num_dofs_per_proc[rank];
    level->deallocatePatchData(u_dof_index_idx);
    level->deallocatePatchData(p_dof_index_idx);

    auto run_mode = [&](const std::string& mode,
                        const std::string& policy,
                        const IntVector<NDIM>& mode_box_size,
                        const IntVector<NDIM>& mode_overlap_size,
                        const int mode_seed_axis,
                        const int mode_seed_stride) -> std::pair<std::vector<std::set<int>>, std::vector<std::set<int>>>
    {
        Pointer<MemoryDatabase> solver_db = new MemoryDatabase("solver_db");
        solver_db->putString("pc_type", "asm");
        solver_db->putString("ksp_type", "richardson");
        solver_db->putInteger("max_iterations", 1);
        int box_arr[NDIM];
        int overlap_arr[NDIM];
        for (int d = 0; d < NDIM; ++d)
        {
            box_arr[d] = mode_box_size(d);
            overlap_arr[d] = mode_overlap_size(d);
        }
        solver_db->putIntegerArray("subdomain_box_size", box_arr, NDIM);
        solver_db->putIntegerArray("subdomain_overlap_size", overlap_arr, NDIM);
        solver_db->putString("asm_subdomain_construction_mode", mode);
        solver_db->putInteger("coupling_aware_asm_seed_axis", mode_seed_axis);
        solver_db->putInteger("coupling_aware_asm_seed_stride", mode_seed_stride);
        solver_db->putString("coupling_aware_asm_closure_policy", policy);

        Pointer<IBAMR::StaggeredStokesPETScLevelSolver> solver =
            new IBAMR::StaggeredStokesPETScLevelSolver("solver_" + mode, solver_db, "stokes_mode_");
        PoissonSpecifications problem_coefs("stokes_mode_poisson");
        problem_coefs.setCConstant(1.0);
        problem_coefs.setDConstant(-1.0);
        solver->setVelocityPoissonSpecifications(problem_coefs);
        solver->initializeSolverState(x_vec, b_vec);

        std::vector<IS>* nonoverlap_is_ptr = nullptr;
        std::vector<IS>* overlap_is_ptr = nullptr;
        solver->getASMSubdomains(&nonoverlap_is_ptr, &overlap_is_ptr);
        std::vector<std::set<int>> overlap_sets = collect_is_sets(*overlap_is_ptr);
        std::vector<std::set<int>> nonoverlap_sets = collect_is_sets(*nonoverlap_is_ptr);

        solver->deallocateSolverState();
        return std::make_pair(overlap_sets, nonoverlap_sets);
    };

    const auto run = run_mode(asm_mode, closure_policy, box_size, overlap_size, seed_axis, seed_stride);
    check_partition_invariants(test_failures, run.first, run.second, n_local_dofs);
    const SubdomainStats stats = compute_stats(run.first, run.second);
    bool compare_overlap_equal = false;
    bool compare_nonoverlap_equal = false;
    SubdomainStats compare_stats;
    if (run_comparison)
    {
        const auto compare_run = run_mode(compare_asm_mode,
                                          compare_closure_policy,
                                          compare_box_size,
                                          compare_overlap_size,
                                          compare_seed_axis,
                                          compare_seed_stride);
        check_partition_invariants(test_failures, compare_run.first, compare_run.second, n_local_dofs);
        compare_stats = compute_stats(compare_run.first, compare_run.second);
        compare_overlap_equal =
            (canonicalize_subdomain_sets(run.first) == canonicalize_subdomain_sets(compare_run.first));
        compare_nonoverlap_equal =
            (canonicalize_subdomain_sets(run.second) == canonicalize_subdomain_sets(compare_run.second));
        if (!compare_overlap_equal || !compare_nonoverlap_equal) ++test_failures;
    }

    level->deallocatePatchData(u_idx);
    level->deallocatePatchData(p_idx);
    level->deallocatePatchData(f_u_idx);
    level->deallocatePatchData(f_p_idx);

    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "asm_subdomain_construction_mode = " << asm_mode << "\n";
    pout << "coupling_aware_asm_closure_policy = " << closure_policy << "\n";
    pout << "subdomain_box_size = " << box_size << "\n";
    pout << "subdomain_overlap_size = " << overlap_size << "\n";
    pout << "coupling_aware_asm_seed_axis = " << seed_axis << "\n";
    pout << "coupling_aware_asm_seed_stride = " << seed_stride << "\n";
    pout << "n_local_dofs = " << n_local_dofs << "\n";
    pout << "n_subdomains = " << stats.n_subdomains << "\n";
    pout << "overlap_sizes = " << stats.overlap_min << " " << stats.overlap_max << " " << stats.overlap_avg << "\n";
    pout << "nonoverlap_sizes = " << stats.nonoverlap_min << " " << stats.nonoverlap_max << " " << stats.nonoverlap_avg
         << "\n";
    pout << "overlap_hash = " << stats.overlap_hash << "\n";
    pout << "nonoverlap_hash = " << stats.nonoverlap_hash << "\n";
    pout << "compare_against_mode = " << (run_comparison ? 1 : 0) << "\n";
    if (run_comparison)
    {
        pout << "compare_asm_subdomain_construction_mode = " << compare_asm_mode << "\n";
        pout << "compare_coupling_aware_asm_closure_policy = " << compare_closure_policy << "\n";
        pout << "compare_subdomain_box_size = " << compare_box_size << "\n";
        pout << "compare_subdomain_overlap_size = " << compare_overlap_size << "\n";
        pout << "compare_coupling_aware_asm_seed_axis = " << compare_seed_axis << "\n";
        pout << "compare_coupling_aware_asm_seed_stride = " << compare_seed_stride << "\n";
        pout << "compare_n_subdomains = " << compare_stats.n_subdomains << "\n";
        pout << "compare_overlap_sizes = " << compare_stats.overlap_min << " " << compare_stats.overlap_max << " "
             << compare_stats.overlap_avg << "\n";
        pout << "compare_nonoverlap_sizes = " << compare_stats.nonoverlap_min << " " << compare_stats.nonoverlap_max
             << " " << compare_stats.nonoverlap_avg << "\n";
        pout << "compare_overlap_hash = " << compare_stats.overlap_hash << "\n";
        pout << "compare_nonoverlap_hash = " << compare_stats.nonoverlap_hash << "\n";
        pout << "compare_overlap_equal = " << (compare_overlap_equal ? 1 : 0) << "\n";
        pout << "compare_nonoverlap_equal = " << (compare_nonoverlap_equal ? 1 : 0) << "\n";
    }
    pout << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
