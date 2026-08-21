// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2026 by the IBAMR developers
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
#include <ibamr/private/StaggeredStokesEigenSchurComplementShellBackend.h>

#include <tbox/Database.h>

#include <petscvec.h>

namespace IBAMR
{
namespace
{
std::unique_ptr<IBTK::PETScLevelSolverShellBackend>
allocate_eigen_schur_complement_shell_backend(IBTK::PETScLevelSolver& solver)
{
    auto* stokes_solver = dynamic_cast<StaggeredStokesPETScLevelSolver*>(&solver);
    if (!stokes_solver) return nullptr;
    return std::make_unique<StaggeredStokesEigenSchurComplementShellBackend>(*stokes_solver);
}

Eigen::MatrixXd
extract_dense_block(const Eigen::MatrixXd& matrix,
                    const std::vector<int>& row_positions,
                    const std::vector<int>& col_positions)
{
    Eigen::MatrixXd block(static_cast<Eigen::Index>(row_positions.size()),
                          static_cast<Eigen::Index>(col_positions.size()));
    for (std::size_t row_idx = 0; row_idx < row_positions.size(); ++row_idx)
    {
        for (std::size_t col_idx = 0; col_idx < col_positions.size(); ++col_idx)
        {
            block(static_cast<Eigen::Index>(row_idx), static_cast<Eigen::Index>(col_idx)) =
                matrix(row_positions[row_idx], col_positions[col_idx]);
        }
    }
    return block;
}

class StaggeredStokesEigenSchurComplementShellBackendRegistrar
{
public:
    StaggeredStokesEigenSchurComplementShellBackendRegistrar()
    {
        IBTK::PETScLevelSolverShellBackendManager::getManager()->registerShellBackendFactoryFunction(
            "eigen-schur-complement", allocate_eigen_schur_complement_shell_backend);
    }
};

static StaggeredStokesEigenSchurComplementShellBackendRegistrar eigen_schur_complement_shell_backend_registrar;
} // namespace

const std::string StaggeredStokesEigenSchurComplementShellBackend::s_backend_name = "EigenSchurComplement";

StaggeredStokesEigenSchurComplementShellBackend::StaggeredStokesEigenSchurComplementShellBackend(
    StaggeredStokesPETScLevelSolver& solver)
    : d_stokes_solver(solver)
{
}

void
StaggeredStokesEigenSchurComplementShellBackend::configureFromInputDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    d_a00_solver_type = EigenSubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR;
    d_a00_solver_threshold = -1.0;
    d_schur_solver_type = EigenSubdomainSolverType::FULL_PIV_HOUSEHOLDER_QR;
    d_schur_solver_threshold = -1.0;
    if (input_db && input_db->keyExists("a00_solver_type"))
    {
        d_a00_solver_type = parseBuiltinSolverType(input_db->getString("a00_solver_type"));
    }
    if (input_db && input_db->keyExists("a00_solver_threshold"))
    {
        d_a00_solver_threshold = input_db->getDouble("a00_solver_threshold");
    }
    if (input_db && input_db->keyExists("schur_solver_type"))
    {
        d_schur_solver_type = parseBuiltinSolverType(input_db->getString("schur_solver_type"));
    }
    if (input_db && input_db->keyExists("schur_solver_threshold"))
    {
        d_schur_solver_threshold = input_db->getDouble("schur_solver_threshold");
    }
}

const std::string&
StaggeredStokesEigenSchurComplementShellBackend::getName() const
{
    return s_backend_name;
}

StaggeredStokesEigenSchurComplementShellBackend::EigenSubdomainSolverType
StaggeredStokesEigenSchurComplementShellBackend::parseBuiltinSolverType(const std::string& type) const
{
    const auto solver_type = parseEigenSubdomainSolverType(
        type, "StaggeredStokesEigenSchurComplementShellBackend::parseBuiltinSolverType()");
    return solver_type;
}

Eigen::MatrixXd
StaggeredStokesEigenSchurComplementShellBackend::buildSchurSolveMatrix(const Eigen::MatrixXd& schur) const
{
    return buildDenseEigenSolveMatrix(d_schur_solver_type,
                                      schur,
                                      d_schur_solver_threshold,
                                      "StaggeredStokesEigenSchurComplementShellBackend::buildSchurSolveMatrix()");
}

void
StaggeredStokesEigenSchurComplementShellBackend::initializeSolverState(
    const IBTK::PETScLevelSolverShellBackendState& solver_state)
{
    checkSerialEigenShellBackend(solver_state,
                                 "StaggeredStokesEigenSchurComplementShellBackend::initializeSolverState()");
    initializeCorrectionCompositionState(solver_state);
    configureFromInputDatabase(d_solver_state.input_db);
    const std::size_t n_subdomains = d_solver_state.subdomain_dofs->size();
    if (d_subdomain_caches.size() != n_subdomains)
    {
        d_subdomain_caches.resize(n_subdomains);
    }

    auto initialize_impl = [this, n_subdomains](auto a00_solver_tag)
    {
        using SolverType = typename decltype(a00_solver_tag)::type;
        initializeCustomEigenA00SolveStorage<SolverType>(n_subdomains);
        initializeCommonDataWithLocalOperatorHook(
            [this](const Eigen::MatrixXd& local_operator, const std::size_t subdomain_num)
            {
                const auto& overlap_dofs = getCommonSubdomains()[subdomain_num].overlap_dofs;
                auto& cache = d_subdomain_caches[subdomain_num];
                cache = CustomEigenSchurSubdomainCache();

                for (std::size_t local_pos = 0; local_pos < overlap_dofs.size(); ++local_pos)
                {
                    const int dof = overlap_dofs[local_pos];
                    if (d_stokes_solver.isVelocityDOF(dof))
                    {
                        cache.velocity_positions.push_back(static_cast<int>(local_pos));
                    }
                    else if (d_stokes_solver.isPressureDOF(dof))
                    {
                        cache.pressure_positions.push_back(static_cast<int>(local_pos));
                    }
                    else
                    {
                        TBOX_ERROR("StaggeredStokesEigenSchurComplementShellBackend::initializeSolverState():\n"
                                   << "  unable to classify local overlap DOF " << dof
                                   << " as velocity or pressure.\n");
                    }
                }

                cache.A00 = extract_dense_block(local_operator, cache.velocity_positions, cache.velocity_positions);
                cache.A01 = extract_dense_block(local_operator, cache.velocity_positions, cache.pressure_positions);
                cache.A10 = extract_dense_block(local_operator, cache.pressure_positions, cache.velocity_positions);
                cache.A11 = extract_dense_block(local_operator, cache.pressure_positions, cache.pressure_positions);

                if (!cache.A00.size() && !cache.A11.size())
                {
                    TBOX_ERROR("StaggeredStokesEigenSchurComplementShellBackend::initializeSolverState():\n"
                               << "  local custom Schur subdomain has no velocity or pressure DOFs.\n");
                }

                if (cache.A00.rows() > 0)
                {
                    auto& solver = getCustomEigenA00SolveStorage<SolverType>().solvers[subdomain_num];
                    initializeCustomEigenA00Solver(solver, cache.A00, subdomain_num);
                }

                if (cache.A11.rows() > 0)
                {
                    if (cache.A00.rows() > 0)
                    {
                        const auto& solver = getCustomEigenA00SolveStorage<SolverType>().solvers[subdomain_num];
                        cache.A00_inv_A01 =
                            StaggeredStokesEigenSchurComplementShellBackendDetail::solve_custom_a00(solver, cache.A01);
                        cache.schur = cache.A11 - cache.A10 * cache.A00_inv_A01;
                    }
                    else
                    {
                        cache.schur = cache.A11;
                    }
                    cache.schur_solve_matrix = buildSchurSolveMatrix(cache.schur);
                }

                cache.velocity_rhs_workspace.resize(static_cast<Eigen::Index>(cache.velocity_positions.size()));
                cache.pressure_rhs_workspace.resize(static_cast<Eigen::Index>(cache.pressure_positions.size()));
                cache.velocity_solution_workspace.resize(static_cast<Eigen::Index>(cache.velocity_positions.size()));
                cache.pressure_solution_workspace.resize(static_cast<Eigen::Index>(cache.pressure_positions.size()));
            });
    };

    dispatchEigenSolverType(d_a00_solver_type, [&initialize_impl](auto solver_tag) { initialize_impl(solver_tag); });
}

void
StaggeredStokesEigenSchurComplementShellBackend::deallocateSolverState()
{
    deallocateCorrectionCompositionState();
    clearCommonData();
    d_subdomain_caches.clear();
    d_a00_solver_storage.reset();
}

void
StaggeredStokesEigenSchurComplementShellBackend::solveSubdomain(const std::size_t subdomain_num)
{
    auto& common_cache = getCommonSubdomains()[subdomain_num];
    auto& custom_cache = d_subdomain_caches[subdomain_num];
    d_a00_solver_storage->solveSubdomain(
        *this, custom_cache, common_cache.rhs_workspace, common_cache.delta_workspace, subdomain_num);
}
} // namespace IBAMR
