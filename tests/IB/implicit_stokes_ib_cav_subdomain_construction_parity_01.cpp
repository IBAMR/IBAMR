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

#include <ibamr/IBImplicitStaggeredHierarchyIntegrator.h>
#include <ibamr/IBMethod.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/StaggeredStokesPETScMatUtilities.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LData.h>
#include <ibtk/LDataManager.h>

#include <petscmat.h>

#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <CellVariable.h>
#include <IntVector.h>
#include <LoadBalancer.h>
#include <PatchHierarchy.h>
#include <PatchLevel.h>
#include <PoissonSpecifications.h>
#include <SAMRAI_config.h>
#include <SideVariable.h>
#include <StandardTagAndInitialize.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <vector>

#include <ibamr/app_namespaces.h>

namespace
{
constexpr double L_DOMAIN = 1.0;
constexpr int N_GRID = 32;
constexpr double R_CYL = 0.25;
constexpr double ALPHA = 0.23;
constexpr double BETA = (R_CYL * R_CYL) / ALPHA;
constexpr double RHO = 1.0;
constexpr double MU = 1.0e-2;
constexpr double K_SPRING = 1.0e2;
constexpr double DX = L_DOMAIN / static_cast<double>(N_GRID);
constexpr double DT = 0.5 * DX;

int s_finest_ln = 0;

void
ib4_interp_fcn(const double r, double* const w)
{
    const double q = std::sqrt(-7.0 + 12.0 * r - 4.0 * r * r);
    w[0] = 0.125 * (5.0 - 2.0 * r - q);
    w[1] = 0.125 * (5.0 - 2.0 * r + q);
    w[2] = 0.125 * (-1.0 + 2.0 * r + q);
    w[3] = 0.125 * (-1.0 + 2.0 * r - q);
}

int
compute_num_lag_nodes()
{
    const double ds0 = DX / R_CYL;
    return static_cast<int>(std::llround((2.0 * M_PI) / ds0));
}

void
generate_structure(const unsigned int& strct_num,
                   const int& ln,
                   int& num_vertices,
                   std::vector<IBTK::Point>& vertex_posn,
                   void*)
{
    if (ln != s_finest_ln || strct_num != 0)
    {
        num_vertices = 0;
        vertex_posn.clear();
        return;
    }

    const int n_nodes = compute_num_lag_nodes();
    num_vertices = n_nodes;
    vertex_posn.resize(num_vertices);
    for (int k = 0; k < num_vertices; ++k)
    {
        const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_vertices);
        vertex_posn[k](0) = 0.5 + ALPHA * std::cos(theta);
        vertex_posn[k](1) = 0.5 + BETA * std::sin(theta);
    }
}

void
generate_springs(
    const unsigned int& strct_num,
    const int& ln,
    std::multimap<int, IBRedundantInitializer::Edge>& spring_map,
    std::map<IBRedundantInitializer::Edge, IBRedundantInitializer::SpringSpec, IBRedundantInitializer::EdgeComp>&
        spring_spec,
    void*)
{
    if (ln != s_finest_ln || strct_num != 0) return;

    const int n_nodes = compute_num_lag_nodes();
    const double ds = (2.0 * M_PI) / static_cast<double>(n_nodes);
    const double spring_k = K_SPRING / (ds * ds);
    for (int k = 0; k < n_nodes; ++k)
    {
        IBRedundantInitializer::Edge edge = { k, (k + 1) % n_nodes };
        if (edge.first > edge.second) std::swap(edge.first, edge.second);
        spring_map.insert(std::make_pair(edge.first, edge));

        IBRedundantInitializer::SpringSpec spec_data;
        spec_data.force_fcn_idx = 0;
        spec_data.parameters.resize(2);
        spec_data.parameters[0] = spring_k;
        spec_data.parameters[1] = 0.0;
        spring_spec.insert(std::make_pair(edge, spec_data));
    }
}

std::set<int>
matlab_extract_coupled_dofs_relaxed(const int seed_velocity_dof,
                                    Mat A00_mat,
                                    const std::map<int, std::set<int>>& velocity_dof_to_adjacent_cell_dofs,
                                    const std::map<int, std::set<int>>& cell_dof_to_closure_dofs,
                                    const std::map<int, int>& velocity_dof_to_component_axis)
{
    std::set<int> initial_velocity_dofs = { seed_velocity_dof };
    PetscInt ncols = 0;
    const PetscInt* cols = nullptr;
    const PetscInt row = static_cast<PetscInt>(seed_velocity_dof);
    int ierr = MatGetRow(A00_mat, row, &ncols, &cols, nullptr);
    IBTK_CHKERRQ(ierr);
    for (PetscInt k = 0; k < ncols; ++k) initial_velocity_dofs.insert(static_cast<int>(cols[k]));
    ierr = MatRestoreRow(A00_mat, row, &ncols, &cols, nullptr);
    IBTK_CHKERRQ(ierr);

    std::set<int> involved_cells;
    for (const int vel_dof : initial_velocity_dofs)
    {
        const auto it = velocity_dof_to_adjacent_cell_dofs.find(vel_dof);
        if (it == velocity_dof_to_adjacent_cell_dofs.end()) continue;
        involved_cells.insert(it->second.begin(), it->second.end());
    }

    std::set<int> closure;
    for (const int cell_dof : involved_cells)
    {
        const auto it = cell_dof_to_closure_dofs.find(cell_dof);
        if (it == cell_dof_to_closure_dofs.end()) continue;
        closure.insert(it->second.begin(), it->second.end());
    }

    std::set<int> out = closure;
    for (const int dof : initial_velocity_dofs)
    {
        if (velocity_dof_to_component_axis.find(dof) != velocity_dof_to_component_axis.end()) out.insert(dof);
    }
    return out;
}

std::set<int>
matlab_extract_coupled_dofs_strict(const int seed_velocity_dof,
                                   Mat A00_mat,
                                   const std::map<int, std::set<int>>& velocity_dof_to_adjacent_cell_dofs,
                                   const std::map<int, std::set<int>>& cell_dof_to_closure_dofs,
                                   const std::map<int, int>& velocity_dof_to_component_axis,
                                   const std::map<int, std::set<int>>& velocity_dof_to_paired_seed_velocity_dofs)
{
    std::set<int> initial_seed_components = { seed_velocity_dof };
    const auto pair_it = velocity_dof_to_paired_seed_velocity_dofs.find(seed_velocity_dof);
    if (pair_it != velocity_dof_to_paired_seed_velocity_dofs.end())
    {
        initial_seed_components.insert(pair_it->second.begin(), pair_it->second.end());
    }

    std::set<int> initial_velocity_dofs = initial_seed_components;
    for (const int seed_component_dof : initial_seed_components)
    {
        PetscInt ncols = 0;
        const PetscInt* cols = nullptr;
        const PetscInt row = static_cast<PetscInt>(seed_component_dof);
        int ierr = MatGetRow(A00_mat, row, &ncols, &cols, nullptr);
        IBTK_CHKERRQ(ierr);
        for (PetscInt k = 0; k < ncols; ++k) initial_velocity_dofs.insert(static_cast<int>(cols[k]));
        ierr = MatRestoreRow(A00_mat, row, &ncols, &cols, nullptr);
        IBTK_CHKERRQ(ierr);
    }

    std::set<int> candidate_cells;
    for (const int vel_dof : initial_velocity_dofs)
    {
        const auto it = velocity_dof_to_adjacent_cell_dofs.find(vel_dof);
        if (it == velocity_dof_to_adjacent_cell_dofs.end()) continue;
        candidate_cells.insert(it->second.begin(), it->second.end());
    }

    std::set<int> accepted_cells;
    for (const int cell_dof : candidate_cells)
    {
        const auto closure_it = cell_dof_to_closure_dofs.find(cell_dof);
        if (closure_it == cell_dof_to_closure_dofs.end()) continue;
        bool valid = true;
        for (const int dof : closure_it->second)
        {
            if (velocity_dof_to_component_axis.find(dof) == velocity_dof_to_component_axis.end()) continue;
            if (initial_velocity_dofs.find(dof) == initial_velocity_dofs.end())
            {
                valid = false;
                break;
            }
        }
        if (valid) accepted_cells.insert(cell_dof);
    }

    std::set<int> out;
    for (const int cell_dof : accepted_cells)
    {
        const auto it = cell_dof_to_closure_dofs.find(cell_dof);
        if (it == cell_dof_to_closure_dofs.end()) continue;
        out.insert(it->second.begin(), it->second.end());
    }
    return out;
}
} // namespace

int
main(int argc, char* argv[])
{
    IBTK::IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    if (IBTK_MPI::getNodes() != 1)
    {
        TBOX_ERROR("implicit_stokes_ib_cav_subdomain_construction_parity_01 requires serial execution (np=1).\n");
    }

#ifndef IBTK_HAVE_SILO
    SAMRAI::tbox::Logger::getInstance()->setWarning(false);
#endif

    Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "output");
    Pointer<Database> input_db = app_initializer->getInputDatabase();

    Pointer<INSStaggeredHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
        "INSStaggeredHierarchyIntegrator", app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    Pointer<IBMethod> ib_method_ops = new IBMethod("IBMethod", app_initializer->getComponentDatabase("IBMethod"));
    ib_method_ops->setUseFixedLEOperators(true);

    Pointer<IBHierarchyIntegrator> time_integrator =
        new IBImplicitStaggeredHierarchyIntegrator("IBHierarchyIntegrator",
                                                   app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                                   ib_method_ops,
                                                   navier_stokes_integrator);

    Pointer<CartesianGridGeometry<NDIM>> grid_geometry = new CartesianGridGeometry<NDIM>(
        "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM>> patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM>> error_detector = new StandardTagAndInitialize<NDIM>(
        "StandardTagAndInitialize", time_integrator, app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM>> box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM>> load_balancer =
        new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM>> gridding_algorithm =
        new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                    app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                    error_detector,
                                    box_generator,
                                    load_balancer);

    Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
        "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
    s_finest_ln = input_db->getInteger("MAX_LEVELS") - 1;
    ib_initializer->setStructureNamesOnLevel(s_finest_ln, { "parity_curve2d" });
    ib_initializer->registerInitStructureFunction(generate_structure);
    ib_initializer->registerInitSpringDataFunction(generate_springs);
    ib_method_ops->registerLInitStrategy(ib_initializer);
    Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
    ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

    time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
    const double current_time = time_integrator->getIntegratorTime();
    const double new_time = current_time + time_integrator->getMaximumTimeStepSize();
    const int num_cycles = time_integrator->getNumberOfCycles();
    time_integrator->preprocessIntegrateHierarchy(current_time, new_time, num_cycles);

    const int finest_ln = patch_hierarchy->getFinestLevelNumber();
    Pointer<PatchLevel<NDIM>> level = patch_hierarchy->getPatchLevel(finest_ln);
    const double data_time = new_time;

    Mat A = nullptr;
    ib_method_ops->updateFixedLEOperators();
    ib_method_ops->constructLagrangianForceJacobian(A, MATAIJ, data_time);

    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    Pointer<VariableContext> ctx = var_db->getContext("cav_subdomain_parity_ctx");
    Pointer<SideVariable<NDIM, int>> u_dof_index_var = new SideVariable<NDIM, int>("cav_parity_u_dof");
    Pointer<CellVariable<NDIM, int>> p_dof_index_var = new CellVariable<NDIM, int>("cav_parity_p_dof");
    const int u_dof_index_idx = var_db->registerVariableAndContext(u_dof_index_var, ctx, IntVector<NDIM>(1));
    const int p_dof_index_idx = var_db->registerVariableAndContext(p_dof_index_var, ctx, IntVector<NDIM>(1));
    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> ln_level = patch_hierarchy->getPatchLevel(ln);
        ln_level->allocatePatchData(u_dof_index_idx, data_time);
        ln_level->allocatePatchData(p_dof_index_idx, data_time);
    }

    std::vector<int> num_dofs_per_proc;
    IBAMR::StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);

    Mat J = nullptr;
    ib_method_ops->constructInterpOp(J, ib4_interp_fcn, 4, num_dofs_per_proc, u_dof_index_idx, data_time);

    Mat SAJ = nullptr;
    int ierr = MatPtAP(A, J, MAT_INITIAL_MATRIX, 1.0, &SAJ);
    IBTK_CHKERRQ(ierr);

    SAMRAI::solv::PoissonSpecifications U_problem_coefs("U_problem_coefs");
    U_problem_coefs.setCConstant(RHO / DT);
    U_problem_coefs.setDConstant(-MU);
    std::vector<SAMRAI::solv::RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM, nullptr);

    Mat stokes_mat = nullptr;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(
        stokes_mat, U_problem_coefs, u_bc_coefs, data_time, num_dofs_per_proc, u_dof_index_idx, p_dof_index_idx, level);

    Mat A00_full = nullptr;
    ierr = MatDuplicate(stokes_mat, MAT_COPY_VALUES, &A00_full);
    IBTK_CHKERRQ(ierr);
    ierr = MatAXPY(A00_full, 1.0, SAJ, DIFFERENT_NONZERO_PATTERN);
    IBTK_CHKERRQ(ierr);

    IBAMR::StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData map_data;
    IBAMR::StaggeredStokesPETScMatUtilities::buildPatchLevelCellClosureMaps(
        map_data, u_dof_index_idx, p_dof_index_idx, level);

    std::vector<std::set<int>> overlap_relaxed, nonoverlap_relaxed;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
        overlap_relaxed,
        nonoverlap_relaxed,
        num_dofs_per_proc,
        u_dof_index_idx,
        level,
        Pointer<CoarseFineBoundary<NDIM>>(nullptr),
        A00_full,
        map_data,
        0,
        1,
        IBAMR::CouplingAwareASMClosurePolicy::RELAXED);

    std::vector<std::set<int>> overlap_strict, nonoverlap_strict;
    IBAMR::StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
        overlap_strict,
        nonoverlap_strict,
        num_dofs_per_proc,
        u_dof_index_idx,
        level,
        Pointer<CoarseFineBoundary<NDIM>>(nullptr),
        A00_full,
        map_data,
        0,
        1,
        IBAMR::CouplingAwareASMClosurePolicy::STRICT);

    int test_failures = 0;
    if (overlap_relaxed.size() != overlap_strict.size())
    {
        ++test_failures;
    }

    std::set<int> axis_velocity_dofs;
    for (PatchLevel<NDIM>::Iterator p(level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = level->getPatch(p());
        const Box<NDIM>& patch_box = patch->getBox();
        const Box<NDIM> side_patch_box = SideGeometry<NDIM>::toSideBox(patch_box, 0);
        Pointer<SideData<NDIM, int>> u_dof_data = patch->getPatchData(u_dof_index_idx);
        for (Box<NDIM>::Iterator b(side_patch_box); b; b++)
        {
            const SideIndex<NDIM> i_s(b(), 0, SideIndex<NDIM>::Lower);
            const int dof = (*u_dof_data)(i_s);
            if (dof < 0) continue;
            const auto axis_it = map_data.velocity_dof_to_component_axis.find(dof);
            if (axis_it == map_data.velocity_dof_to_component_axis.end()) continue;
            if (axis_it->second != 0) continue;
            axis_velocity_dofs.insert(dof);
        }
    }

    std::vector<int> ordered_seeds;
    ordered_seeds.reserve(axis_velocity_dofs.size());
    for (const int dof : axis_velocity_dofs) ordered_seeds.push_back(dof);

    if (overlap_relaxed.size() != ordered_seeds.size() || overlap_strict.size() != ordered_seeds.size())
    {
        ++test_failures;
    }

    for (std::size_t k = 0; k < overlap_relaxed.size() && k < ordered_seeds.size(); ++k)
    {
        const int seed = ordered_seeds[k];
        const std::set<int> matlab_relaxed =
            matlab_extract_coupled_dofs_relaxed(seed,
                                                A00_full,
                                                map_data.velocity_dof_to_adjacent_cell_dofs,
                                                map_data.cell_dof_to_closure_dofs,
                                                map_data.velocity_dof_to_component_axis);
        if (overlap_relaxed[k] != matlab_relaxed)
        {
            ++test_failures;
            break;
        }
    }

    for (std::size_t k = 0; k < overlap_strict.size() && k < ordered_seeds.size(); ++k)
    {
        const int seed = ordered_seeds[k];
        const std::set<int> matlab_strict =
            matlab_extract_coupled_dofs_strict(seed,
                                               A00_full,
                                               map_data.velocity_dof_to_adjacent_cell_dofs,
                                               map_data.cell_dof_to_closure_dofs,
                                               map_data.velocity_dof_to_component_axis,
                                               map_data.velocity_dof_to_paired_seed_velocity_dofs);
        if (overlap_strict[k] != matlab_strict)
        {
            ++test_failures;
            break;
        }
    }

    ierr = MatDestroy(&A00_full);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&stokes_mat);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&SAJ);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&J);
    IBTK_CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    IBTK_CHKERRQ(ierr);

    for (int ln = 0; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM>> ln_level = patch_hierarchy->getPatchLevel(ln);
        ln_level->deallocatePatchData(u_dof_index_idx);
        ln_level->deallocatePatchData(p_dof_index_idx);
    }

    plog << "Input database:\n";
    input_db->printClassData(plog);
    pout << "test_failures = " << test_failures << "\n";
    return test_failures > 0 ? 1 : 0;
}
