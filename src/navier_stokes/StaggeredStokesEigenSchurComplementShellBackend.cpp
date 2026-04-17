// ---------------------------------------------------------------------
//
// Copyright (c) 2014 - 2025 by the IBAMR developers
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

#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/Database.h>

#include <petscvec.h>

namespace IBAMR
{
namespace
{
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

} // namespace

StaggeredStokesEigenSchurComplementShellBackend::StaggeredStokesEigenSchurComplementShellBackend(
    StaggeredStokesPETScLevelSolver& solver,
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
    : PETScLevelSolverEigenShellBackendBase(solver), d_stokes_solver(solver)
{
    configure(input_db);
}

void
StaggeredStokesEigenSchurComplementShellBackend::configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
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
StaggeredStokesEigenSchurComplementShellBackend::getTypeKey() const
{
    return d_type_key;
}

const char*
StaggeredStokesEigenSchurComplementShellBackend::getPCNameSuffixAdditive() const
{
    return "PC_EigenSchurComplementAdditive";
}

const char*
StaggeredStokesEigenSchurComplementShellBackend::getPCNameSuffixMultiplicative() const
{
    return "PC_EigenSchurComplement";
}

void
StaggeredStokesEigenSchurComplementShellBackend::reset()
{
    clearCommonData();
    d_subdomain_caches.clear();
    d_a00_solver_storage.reset();
}

void
StaggeredStokesEigenSchurComplementShellBackend::deallocate()
{
    reset();
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
    Eigen::MatrixXd solve_matrix;
    dispatchEigenSolverType(
        d_schur_solver_type,
        [this, &schur, &solve_matrix](auto solver_tag)
        {
            using SolverType = typename decltype(solver_tag)::type;
            if constexpr (std::is_same_v<SolverType, Eigen::LLT<Eigen::MatrixXd>>)
            {
                solve_matrix = buildLLTSolveMatrix(schur);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::LDLT<Eigen::MatrixXd>>)
            {
                solve_matrix = buildLDLTSolveMatrix(schur);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::PartialPivLU<Eigen::MatrixXd>>)
            {
                solve_matrix = buildPartialPivLUSolveMatrix(schur);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::FullPivLU<Eigen::MatrixXd>>)
            {
                solve_matrix = buildFullPivLUSolveMatrix(schur, d_schur_solver_threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::HouseholderQR<Eigen::MatrixXd>>)
            {
                solve_matrix = buildHouseholderQRSolveMatrix(schur);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>)
            {
                solve_matrix = buildQRSolveMatrix<SolverType>(schur, d_schur_solver_threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>>)
            {
                solve_matrix = buildCompleteOrthogonalDecompositionPseudoinverse(schur, d_schur_solver_threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>)
            {
                solve_matrix = buildQRSolveMatrix<SolverType>(schur, d_schur_solver_threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::JacobiSVD<Eigen::MatrixXd>> ||
                               std::is_same_v<SolverType, Eigen::BDCSVD<Eigen::MatrixXd>>)
            {
                solve_matrix = buildSVDPseudoinverse<SolverType>(schur, d_schur_solver_threshold);
            }
        });
    return solve_matrix;
}

void
StaggeredStokesEigenSchurComplementShellBackend::initialize()
{
    const std::size_t n_subdomains = d_context.getSubdomainDOFsForBackend().size();
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
                cache.overlap_size = static_cast<int>(overlap_dofs.size());
                cache.overlap_dofs = &overlap_dofs;

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
                        TBOX_ERROR("StaggeredStokesEigenSchurComplementShellBackend::initialize():\n"
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
                    TBOX_ERROR("StaggeredStokesEigenSchurComplementShellBackend::initialize():\n"
                               << "  local custom Schur subdomain has no velocity or pressure DOFs.\n");
                }

                if (cache.A00.rows() > 0)
                {
                    auto& solver = getCustomEigenA00SolveStorage<SolverType>().solvers[subdomain_num];
                    StaggeredStokesEigenSchurComplementShellBackendDetail::initialize_custom_a00_solver(
                        solver, cache.A00, d_a00_solver_threshold, subdomain_num);
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

                cache.rhs_workspace.resize(static_cast<Eigen::Index>(cache.overlap_size));
                cache.delta_workspace.resize(static_cast<Eigen::Index>(cache.overlap_size));
                cache.velocity_rhs_workspace.resize(static_cast<Eigen::Index>(cache.velocity_positions.size()));
                cache.pressure_rhs_workspace.resize(static_cast<Eigen::Index>(cache.pressure_positions.size()));
                cache.velocity_solution_workspace.resize(static_cast<Eigen::Index>(cache.velocity_positions.size()));
                cache.pressure_solution_workspace.resize(static_cast<Eigen::Index>(cache.pressure_positions.size()));
            });

        const auto& common_subdomains = getCommonSubdomains();
        for (std::size_t subdomain_num = 0; subdomain_num < n_subdomains; ++subdomain_num)
        {
            auto& cache = d_subdomain_caches[subdomain_num];
            const auto& common_cache = common_subdomains[subdomain_num];
            cache.update_dofs = &common_cache.update_dofs;
            cache.update_local_positions = &common_cache.update_local_positions;
            cache.active_residual_update_rows = &common_cache.active_residual_update_rows;
            cache.active_residual_update_mat = &common_cache.active_residual_update_mat;
            cache.residual_input_workspace.resize(static_cast<Eigen::Index>(cache.update_local_positions->size()));
            cache.residual_delta_workspace.resize(static_cast<Eigen::Index>(cache.active_residual_update_rows->size()));
        }
    };

    dispatchEigenSolverType(d_a00_solver_type, [&initialize_impl](auto solver_tag) { initialize_impl(solver_tag); });
}

void
StaggeredStokesEigenSchurComplementShellBackend::applyAdditive(Vec x, Vec y)
{
    dispatchEigenSolverType(d_a00_solver_type,
                            [this, x, y](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                applyAdditiveImpl<SolverType>(x, y);
                            });
}

void
StaggeredStokesEigenSchurComplementShellBackend::applyMultiplicative(Vec x, Vec y)
{
    dispatchEigenSolverType(d_a00_solver_type,
                            [this, x, y](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                applyMultiplicativeImpl<SolverType>(x, y);
                            });
}
} // namespace IBAMR
