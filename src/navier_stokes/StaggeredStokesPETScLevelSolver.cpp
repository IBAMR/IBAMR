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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ibamr/StaggeredStokesPETScLevelSolver.h>
#include <ibamr/StaggeredStokesPETScVecUtilities.h>
#include <ibamr/StaggeredStokesPhysicalBoundaryHelper.h>

#include <ibtk/GeneralSolver.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>
#include <ibtk/LinearSolver.h>
#include <ibtk/PoissonUtilities.h>

#include <tbox/Array.h>
#include <tbox/Database.h>
#include <tbox/Pointer.h>

#include <petscvec.h>

#include <BoundaryBox.h>
#include <CellData.h>
#include <CellVariable.h>
#include <CoarseFineBoundary.h>
#include <IntVector.h>
#include <Patch.h>
#include <PatchGeometry.h>
#include <PatchLevel.h>
#include <SAMRAIVectorReal.h>
#include <SideData.h>
#include <SideVariable.h>
#include <VariableContext.h>
#include <VariableDatabase.h>

#include <algorithm>
#include <string>
#include <type_traits>
#include <vector>

#include <ibamr/namespaces.h> // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
// Number of ghosts cells used for each variable quantity.
static const int CELLG = 1;
static const int SIDEG = 1;
static const int NOGHOST = 0;

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
} // extract_dense_block

Eigen::VectorXd
extract_subvector(const Eigen::VectorXd& vector, const std::vector<int>& positions)
{
    Eigen::VectorXd subvector(static_cast<Eigen::Index>(positions.size()));
    for (std::size_t k = 0; k < positions.size(); ++k)
    {
        subvector(static_cast<Eigen::Index>(k)) = vector(positions[k]);
    }
    return subvector;
} // extract_subvector

void
extract_subvector_into(Eigen::VectorXd& subvector, const Eigen::VectorXd& vector, const std::vector<int>& positions)
{
    TBOX_ASSERT(subvector.size() == static_cast<Eigen::Index>(positions.size()));
    for (std::size_t k = 0; k < positions.size(); ++k)
    {
        subvector(static_cast<Eigen::Index>(k)) = vector(positions[k]);
    }
    return;
} // extract_subvector_into

void
insert_subvector(Eigen::VectorXd& vector, const std::vector<int>& positions, const Eigen::VectorXd& subvector)
{
    for (std::size_t k = 0; k < positions.size(); ++k)
    {
        vector(positions[k]) = subvector(static_cast<Eigen::Index>(k));
    }
    return;
} // insert_subvector

Eigen::MatrixXd
build_llt_solve_matrix(const Eigen::MatrixXd& matrix)
{
    Eigen::LLT<Eigen::MatrixXd> solver(matrix);
    return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
} // build_llt_solve_matrix

Eigen::MatrixXd
build_ldlt_solve_matrix(const Eigen::MatrixXd& matrix)
{
    Eigen::LDLT<Eigen::MatrixXd> solver(matrix);
    return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
} // build_ldlt_solve_matrix

template <class SolverType>
Eigen::MatrixXd
build_qr_solve_matrix(const Eigen::MatrixXd& matrix, const double threshold)
{
    SolverType solver(matrix);
    if (threshold >= 0.0) solver.setThreshold(threshold);
    return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
} // build_qr_solve_matrix

Eigen::MatrixXd
build_householder_qr_solve_matrix(const Eigen::MatrixXd& matrix)
{
    Eigen::HouseholderQR<Eigen::MatrixXd> solver(matrix);
    return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
} // build_householder_qr_solve_matrix

Eigen::MatrixXd
build_partial_piv_lu_solve_matrix(const Eigen::MatrixXd& matrix)
{
    Eigen::PartialPivLU<Eigen::MatrixXd> solver(matrix);
    return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
} // build_partial_piv_lu_solve_matrix

Eigen::MatrixXd
build_full_piv_lu_solve_matrix(const Eigen::MatrixXd& matrix, const double threshold)
{
    Eigen::FullPivLU<Eigen::MatrixXd> solver(matrix);
    if (threshold >= 0.0) solver.setThreshold(threshold);
    return solver.solve(Eigen::MatrixXd::Identity(matrix.rows(), matrix.cols()));
} // build_full_piv_lu_solve_matrix

template <class SVDType>
Eigen::MatrixXd
build_svd_pseudoinverse(const Eigen::MatrixXd& matrix, const double threshold)
{
    SVDType svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    if (threshold >= 0.0) svd.setThreshold(threshold);

    const auto& singular_values = svd.singularValues();
    Eigen::VectorXd inv_singular_values = Eigen::VectorXd::Zero(singular_values.size());
    if (singular_values.size() > 0)
    {
        const double cutoff = singular_values[0] * svd.threshold();
        for (Eigen::Index k = 0; k < singular_values.size(); ++k)
        {
            if (singular_values[k] > cutoff) inv_singular_values[k] = 1.0 / singular_values[k];
        }
    }

    return svd.matrixV() * inv_singular_values.asDiagonal() * svd.matrixU().adjoint();
} // build_svd_pseudoinverse

Eigen::MatrixXd
build_complete_orthogonal_decomposition_pseudoinverse(const Eigen::MatrixXd& matrix, const double threshold)
{
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> solver(matrix);
    if (threshold >= 0.0) solver.setThreshold(threshold);
    Eigen::MatrixXd pseudoinverse = solver.pseudoInverse();
    return pseudoinverse;
} // build_complete_orthogonal_decomposition_pseudoinverse

template <class SolverType>
void
initialize_custom_a00_solver(SolverType& solver,
                             const Eigen::MatrixXd& matrix,
                             const double threshold,
                             const std::size_t subdomain_num)
{
    if constexpr (std::is_same_v<SolverType, Eigen::LLT<Eigen::MatrixXd>>)
    {
        solver.compute(matrix);
        if (solver.info() != Eigen::Success)
        {
            TBOX_ERROR("initialize_custom_a00_solver():\n"
                       << "  LLT factorization failed for the local A00 block on subdomain " << subdomain_num << ".\n");
        }
    }
    else if constexpr (std::is_same_v<SolverType, Eigen::LDLT<Eigen::MatrixXd>>)
    {
        solver.compute(matrix);
        if (solver.info() != Eigen::Success)
        {
            TBOX_ERROR("initialize_custom_a00_solver():\n"
                       << "  LDLT factorization failed for the local A00 block on subdomain " << subdomain_num
                       << ".\n");
        }
    }
    else if constexpr (std::is_same_v<SolverType, Eigen::JacobiSVD<Eigen::MatrixXd>> ||
                       std::is_same_v<SolverType, Eigen::BDCSVD<Eigen::MatrixXd>>)
    {
        solver.compute(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
        if (threshold >= 0.0) solver.setThreshold(threshold);
    }
    else
    {
        solver.compute(matrix);
        if constexpr (std::is_same_v<SolverType, Eigen::FullPivLU<Eigen::MatrixXd>> ||
                      std::is_same_v<SolverType, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>> ||
                      std::is_same_v<SolverType, Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>> ||
                      std::is_same_v<SolverType, Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>)
        {
            if (threshold >= 0.0) solver.setThreshold(threshold);
        }
    }
}

template <class SolverType, class RhsType>
auto
solve_custom_a00(const SolverType& solver, const RhsType& rhs) -> decltype(solver.solve(rhs))
{
    return solver.solve(rhs);
}
} // namespace

/////////////////////////////// PUBLIC ///////////////////////////////////////

template <class SolverType>
StaggeredStokesPETScLevelSolver::CustomEigenA00TypedSolveStorage<SolverType>&
StaggeredStokesPETScLevelSolver::getCustomEigenA00SolveStorage()
{
    auto* storage = static_cast<CustomEigenA00TypedSolveStorage<SolverType>*>(d_custom_eigen_a00_solver_storage.get());
    TBOX_ASSERT(storage);
    return *storage;
}

template <class SolverType>
const StaggeredStokesPETScLevelSolver::CustomEigenA00TypedSolveStorage<SolverType>&
StaggeredStokesPETScLevelSolver::getCustomEigenA00SolveStorage() const
{
    const auto* storage =
        static_cast<const CustomEigenA00TypedSolveStorage<SolverType>*>(d_custom_eigen_a00_solver_storage.get());
    TBOX_ASSERT(storage);
    return *storage;
}

template <class SolverType>
void
StaggeredStokesPETScLevelSolver::initializeCustomEigenA00SolveStorage(const std::size_t n_subdomains)
{
    auto storage = std::make_unique<CustomEigenA00TypedSolveStorage<SolverType>>();
    storage->solvers.resize(n_subdomains);
    d_custom_eigen_a00_solver_storage = std::move(storage);
}

template <class SolverType>
void
StaggeredStokesPETScLevelSolver::solveCustomEigenSubdomain(CustomEigenSchurSubdomainCache& custom_cache,
                                                           const SolverType& a00_solver) const
{
    extract_subvector_into(
        custom_cache.velocity_rhs_workspace, custom_cache.rhs_workspace, custom_cache.velocity_positions);
    extract_subvector_into(
        custom_cache.pressure_rhs_workspace, custom_cache.rhs_workspace, custom_cache.pressure_positions);
    custom_cache.delta_workspace.setZero();

    if (custom_cache.pressure_positions.empty())
    {
        custom_cache.velocity_solution_workspace =
            solve_custom_a00(a00_solver, custom_cache.velocity_rhs_workspace).eval();
        insert_subvector(
            custom_cache.delta_workspace, custom_cache.velocity_positions, custom_cache.velocity_solution_workspace);
    }
    else if (custom_cache.velocity_positions.empty())
    {
        custom_cache.pressure_solution_workspace.noalias() =
            custom_cache.schur_solve_matrix * custom_cache.pressure_rhs_workspace;
        insert_subvector(
            custom_cache.delta_workspace, custom_cache.pressure_positions, custom_cache.pressure_solution_workspace);
    }
    else
    {
        custom_cache.velocity_solution_workspace =
            solve_custom_a00(a00_solver, custom_cache.velocity_rhs_workspace).eval();
        custom_cache.pressure_rhs_workspace -= custom_cache.A10 * custom_cache.velocity_solution_workspace;
        custom_cache.pressure_solution_workspace.noalias() =
            custom_cache.schur_solve_matrix * custom_cache.pressure_rhs_workspace;
        custom_cache.velocity_solution_workspace -= custom_cache.A00_inv_A01 * custom_cache.pressure_solution_workspace;
        insert_subvector(
            custom_cache.delta_workspace, custom_cache.velocity_positions, custom_cache.velocity_solution_workspace);
        insert_subvector(
            custom_cache.delta_workspace, custom_cache.pressure_positions, custom_cache.pressure_solution_workspace);
    }
    return;
} // solveCustomEigenSubdomain

template <class SolverType>
void
StaggeredStokesPETScLevelSolver::applyAdditiveCustomEigenImpl(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    auto& typed_storage = getCustomEigenA00SolveStorage<SolverType>();
    PetscInt n_local = 0;
    int ierr = VecGetLocalSize(x, &n_local);
    IBTK_CHKERRQ(ierr);
    const Eigen::Index n = static_cast<Eigen::Index>(n_local);
    TBOX_ASSERT(n > 0);

    {
        ConstPetscVecArrayMap x_array(x, n);
        PetscVecArrayMap y_array(y, n);
        const auto x_map = x_array.getMap();
        auto y_map = y_array.getMap();
        y_map.setZero();
        for (std::size_t subdomain_num = 0; subdomain_num < d_custom_eigen_schur_subdomain_caches.size();
             ++subdomain_num)
        {
            auto& custom_cache = d_custom_eigen_schur_subdomain_caches[subdomain_num];
            const auto& a00_solver = typed_storage.solvers[subdomain_num];
            const auto& overlap_dofs = *custom_cache.overlap_dofs;
            const auto& update_dofs = *custom_cache.update_dofs;
            const auto& update_local_positions = *custom_cache.update_local_positions;
            std::size_t rhs_idx = 0;
            for (const int dof : overlap_dofs)
            {
                custom_cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = x_map[dof];
            }

            solveCustomEigenSubdomain(custom_cache, a00_solver);

            std::size_t update_idx = 0;
            for (const int dof : update_dofs)
            {
                y_map[dof] +=
                    custom_cache.delta_workspace[static_cast<Eigen::Index>(update_local_positions[update_idx++])];
            }
        }
    }
    postprocessShellResult(y);
    return;
} // applyAdditiveCustomEigenImpl

template <class SolverType>
void
StaggeredStokesPETScLevelSolver::applyMultiplicativeCustomEigenImpl(Vec x, Vec y)
{
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    auto& typed_storage = getCustomEigenA00SolveStorage<SolverType>();
    PetscInt n_local = 0;
    int ierr = VecGetLocalSize(x, &n_local);
    IBTK_CHKERRQ(ierr);
    const Eigen::Index n = static_cast<Eigen::Index>(n_local);
    TBOX_ASSERT(n > 0);

    {
        ConstPetscVecArrayMap x_array(x, n);
        PetscVecArrayMap y_array(y, n);
        const auto x_map = x_array.getMap();
        auto y_map = y_array.getMap();
        Eigen::VectorXd residual(n);
        y_map.setZero();
        residual = x_map;
        const std::size_t n_subdomains = d_custom_eigen_schur_subdomain_caches.size();
        for (std::size_t subdomain_num = 0; subdomain_num + 1 < n_subdomains; ++subdomain_num)
        {
            auto& custom_cache = d_custom_eigen_schur_subdomain_caches[subdomain_num];
            const auto& a00_solver = typed_storage.solvers[subdomain_num];
            const auto& overlap_dofs = *custom_cache.overlap_dofs;
            const auto& update_dofs = *custom_cache.update_dofs;
            const auto& update_local_positions = *custom_cache.update_local_positions;
            const auto& active_residual_update_rows = *custom_cache.active_residual_update_rows;
            const auto& active_residual_update_mat = *custom_cache.active_residual_update_mat;
            std::size_t rhs_idx = 0;
            for (const int dof : overlap_dofs)
            {
                custom_cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }

            solveCustomEigenSubdomain(custom_cache, a00_solver);

            std::size_t update_idx = 0;
            for (const int dof : update_dofs)
            {
                y_map[dof] +=
                    custom_cache.delta_workspace[static_cast<Eigen::Index>(update_local_positions[update_idx++])];
            }
            if (active_residual_update_mat.rows() > 0)
            {
                std::size_t residual_input_idx = 0;
                for (const int local_pos : update_local_positions)
                {
                    custom_cache.residual_input_workspace[static_cast<Eigen::Index>(residual_input_idx++)] =
                        custom_cache.delta_workspace[static_cast<Eigen::Index>(local_pos)];
                }
                custom_cache.residual_delta_workspace.noalias() =
                    active_residual_update_mat * custom_cache.residual_input_workspace;
                std::size_t row_idx = 0;
                for (const int row : active_residual_update_rows)
                {
                    residual[row] -= custom_cache.residual_delta_workspace[static_cast<Eigen::Index>(row_idx++)];
                }
            }
        }
        if (n_subdomains > 0)
        {
            auto& custom_cache = d_custom_eigen_schur_subdomain_caches.back();
            const auto& a00_solver = typed_storage.solvers.back();
            const auto& overlap_dofs = *custom_cache.overlap_dofs;
            const auto& update_dofs = *custom_cache.update_dofs;
            const auto& update_local_positions = *custom_cache.update_local_positions;
            std::size_t rhs_idx = 0;
            for (const int dof : overlap_dofs)
            {
                custom_cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }

            solveCustomEigenSubdomain(custom_cache, a00_solver);

            std::size_t update_idx = 0;
            for (const int dof : update_dofs)
            {
                y_map[dof] +=
                    custom_cache.delta_workspace[static_cast<Eigen::Index>(update_local_positions[update_idx++])];
            }
        }
    }
    postprocessShellResult(y);
    return;
} // applyMultiplicativeCustomEigenImpl

PETScLevelSolver::EigenSubdomainSolverType
StaggeredStokesPETScLevelSolver::parseCustomBuiltinEigenSolverType(const std::string& type,
                                                                   const char* option_name) const
{
    const auto solver_type = parseEigenSubdomainSolverType(type);
    if (solver_type == EigenSubdomainSolverType::CUSTOM)
    {
        TBOX_ERROR("StaggeredStokesPETScLevelSolver::parseCustomBuiltinEigenSolverType():\n"
                   << "  invalid " << option_name << " = " << type << "\n"
                   << "  CUSTOM is reserved for the outer Eigen shell solver selection.\n");
    }
    return solver_type;
} // parseCustomBuiltinEigenSolverType

Eigen::MatrixXd
StaggeredStokesPETScLevelSolver::buildCustomSchurSolveMatrix(const Eigen::MatrixXd& schur) const
{
    Eigen::MatrixXd solve_matrix;
    dispatchEigenSolverType(
        d_custom_eigen_schur_solver_type,
        [&](auto solver_tag)
        {
            using SolverType = typename decltype(solver_tag)::type;
            if constexpr (std::is_same_v<SolverType, Eigen::LLT<Eigen::MatrixXd>>)
            {
                solve_matrix = build_llt_solve_matrix(schur);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::LDLT<Eigen::MatrixXd>>)
            {
                solve_matrix = build_ldlt_solve_matrix(schur);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::PartialPivLU<Eigen::MatrixXd>>)
            {
                solve_matrix = build_partial_piv_lu_solve_matrix(schur);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::FullPivLU<Eigen::MatrixXd>>)
            {
                solve_matrix = build_full_piv_lu_solve_matrix(schur, d_custom_eigen_schur_solver_threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::HouseholderQR<Eigen::MatrixXd>>)
            {
                solve_matrix = build_householder_qr_solve_matrix(schur);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::ColPivHouseholderQR<Eigen::MatrixXd>>)
            {
                solve_matrix = build_qr_solve_matrix<SolverType>(schur, d_custom_eigen_schur_solver_threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>>)
            {
                solve_matrix =
                    build_complete_orthogonal_decomposition_pseudoinverse(schur, d_custom_eigen_schur_solver_threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::FullPivHouseholderQR<Eigen::MatrixXd>>)
            {
                solve_matrix = build_qr_solve_matrix<SolverType>(schur, d_custom_eigen_schur_solver_threshold);
            }
            else if constexpr (std::is_same_v<SolverType, Eigen::JacobiSVD<Eigen::MatrixXd>> ||
                               std::is_same_v<SolverType, Eigen::BDCSVD<Eigen::MatrixXd>>)
            {
                solve_matrix = build_svd_pseudoinverse<SolverType>(schur, d_custom_eigen_schur_solver_threshold);
            }
        });
    return solve_matrix;
} // buildCustomSchurSolveMatrix

StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver(const std::string& object_name,
                                                                 Pointer<Database> input_db,
                                                                 const std::string& default_options_prefix)
{
    GeneralSolver::init(object_name, /*homogeneous_bc*/ false);
    PETScLevelSolver::init(input_db, default_options_prefix);
    if (input_db && input_db->keyExists("subdomain_box_size"))
    {
        input_db->getIntegerArray("subdomain_box_size", d_box_size, NDIM);
    }
    if (input_db && input_db->keyExists("subdomain_overlap_size"))
    {
        input_db->getIntegerArray("subdomain_overlap_size", d_overlap_size, NDIM);
    }
    if (input_db && input_db->keyExists("asm_subdomain_construction_mode"))
    {
        const std::string mode = input_db->getString("asm_subdomain_construction_mode");
        d_asm_subdomain_construction_mode = string_to_enum<ASMSubdomainConstructionMode>(mode);
        if (d_asm_subdomain_construction_mode == ASMSubdomainConstructionMode::UNKNOWN)
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver():\n"
                       << "  invalid asm_subdomain_construction_mode = " << mode << "\n"
                       << "  expected values are " << enum_to_string(ASMSubdomainConstructionMode::GEOMETRICAL)
                       << " and " << enum_to_string(ASMSubdomainConstructionMode::COUPLING_AWARE) << ".\n");
        }
    }
    if (input_db && input_db->keyExists("coupling_aware_asm_seed_axis"))
    {
        d_coupling_aware_asm_seed_axis = input_db->getInteger("coupling_aware_asm_seed_axis");
        if (d_coupling_aware_asm_seed_axis < 0 || d_coupling_aware_asm_seed_axis >= NDIM)
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver():\n"
                       << "  invalid coupling_aware_asm_seed_axis = " << d_coupling_aware_asm_seed_axis << "\n"
                       << "  expected value in [0, " << NDIM - 1 << "].\n");
        }
    }
    if (input_db && input_db->keyExists("coupling_aware_asm_seed_stride"))
    {
        d_coupling_aware_asm_seed_stride = input_db->getInteger("coupling_aware_asm_seed_stride");
        if (d_coupling_aware_asm_seed_stride < 1)
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver():\n"
                       << "  invalid coupling_aware_asm_seed_stride = " << d_coupling_aware_asm_seed_stride << "\n"
                       << "  expected value >= 1.\n");
        }
    }
    if (input_db && input_db->keyExists("coupling_aware_asm_closure_policy"))
    {
        const std::string closure_policy = input_db->getString("coupling_aware_asm_closure_policy");
        d_coupling_aware_asm_closure_policy = string_to_enum<CouplingAwareASMClosurePolicy>(closure_policy);
        if (d_coupling_aware_asm_closure_policy == CouplingAwareASMClosurePolicy::UNKNOWN)
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::StaggeredStokesPETScLevelSolver():\n"
                       << "  invalid coupling_aware_asm_closure_policy = " << closure_policy << "\n"
                       << "  expected values are " << enum_to_string(CouplingAwareASMClosurePolicy::RELAXED) << " and "
                       << enum_to_string(CouplingAwareASMClosurePolicy::STRICT) << ".\n");
        }
    }
    if (input_db && input_db->keyExists("coupling_aware_asm_relative_zero_tol"))
    {
        d_coupling_aware_asm_relative_zero_tol = input_db->getDouble("coupling_aware_asm_relative_zero_tol");
    }
    if (input_db && input_db->keyExists("log_ASM_subdomains"))
    {
        d_log_asm_subdomains = input_db->getBool("log_ASM_subdomains");
    }
    if (input_db && input_db->keyExists("custom_eigen_a00_solver_type"))
    {
        d_custom_eigen_a00_solver_type = parseCustomBuiltinEigenSolverType(
            input_db->getString("custom_eigen_a00_solver_type"), "custom_eigen_a00_solver_type");
    }
    if (input_db && input_db->keyExists("custom_eigen_a00_solver_threshold"))
    {
        d_custom_eigen_a00_solver_threshold = input_db->getDouble("custom_eigen_a00_solver_threshold");
    }
    if (input_db && input_db->keyExists("custom_eigen_schur_solver_type"))
    {
        d_custom_eigen_schur_solver_type = parseCustomBuiltinEigenSolverType(
            input_db->getString("custom_eigen_schur_solver_type"), "custom_eigen_schur_solver_type");
    }
    if (input_db && input_db->keyExists("custom_eigen_schur_solver_threshold"))
    {
        d_custom_eigen_schur_solver_threshold = input_db->getDouble("custom_eigen_schur_solver_threshold");
    }
    // Construct the DOF index variable/context.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_context = var_db->getContext(object_name + "::CONTEXT");
    d_u_dof_index_var = new SideVariable<NDIM, int>(object_name + "::u_dof_index");
    if (var_db->checkVariableExists(d_u_dof_index_var->getName()))
    {
        d_u_dof_index_var = var_db->getVariable(d_u_dof_index_var->getName());
        d_u_dof_index_idx = var_db->mapVariableAndContextToIndex(d_u_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_u_dof_index_idx);
    }
    const int u_gcw = std::max(d_overlap_size.max(), SIDEG);
    d_u_dof_index_idx = var_db->registerVariableAndContext(d_u_dof_index_var, d_context, u_gcw);
    d_p_dof_index_var = new CellVariable<NDIM, int>(object_name + "::p_dof_index");
    if (var_db->checkVariableExists(d_p_dof_index_var->getName()))
    {
        d_p_dof_index_var = var_db->getVariable(d_p_dof_index_var->getName());
        d_p_dof_index_idx = var_db->mapVariableAndContextToIndex(d_p_dof_index_var, d_context);
        var_db->removePatchDataIndex(d_p_dof_index_idx);
    }
    const int p_gcw = std::max(d_overlap_size.max(), CELLG);
    d_p_dof_index_idx = var_db->registerVariableAndContext(d_p_dof_index_var, d_context, p_gcw);

    // Construct the nullspace variable/index.
    d_u_nullspace_var = new SideVariable<NDIM, double>(object_name + "::u_nullspace_var");
    if (var_db->checkVariableExists(d_u_nullspace_var->getName()))
    {
        d_u_nullspace_var = var_db->getVariable(d_u_nullspace_var->getName());
        d_u_nullspace_idx = var_db->mapVariableAndContextToIndex(d_u_nullspace_var, d_context);
        var_db->removePatchDataIndex(d_u_nullspace_idx);
    }
    d_u_nullspace_idx = var_db->registerVariableAndContext(d_u_nullspace_var, d_context, NOGHOST);
    d_p_nullspace_var = new CellVariable<NDIM, double>(object_name + "::p_nullspace_var");
    if (var_db->checkVariableExists(d_p_nullspace_var->getName()))
    {
        d_p_nullspace_var = var_db->getVariable(d_p_nullspace_var->getName());
        d_p_nullspace_idx = var_db->mapVariableAndContextToIndex(d_p_nullspace_var, d_context);
        var_db->removePatchDataIndex(d_p_nullspace_idx);
    }
    d_p_nullspace_idx = var_db->registerVariableAndContext(d_p_nullspace_var, d_context, NOGHOST);

    return;
} // StaggeredStokesPETScLevelSolver

StaggeredStokesPETScLevelSolver::~StaggeredStokesPETScLevelSolver()
{
    if (d_is_initialized) deallocateSolverState();
    return;
} // ~StaggeredStokesPETScLevelSolver

const std::vector<PetscInt>&
StaggeredStokesPETScLevelSolver::getCachedVelocityDOFs() const
{
    return d_velocity_dofs;
}

const std::vector<PetscInt>&
StaggeredStokesPETScLevelSolver::getCachedPressureDOFs() const
{
    return d_pressure_dofs;
}

const StaggeredStokesPETScMatUtilities::PatchLevelCellClosureMapData&
StaggeredStokesPETScLevelSolver::getCouplingAwareASMMapData() const
{
    return d_coupling_aware_asm_map_data;
}

const std::vector<int>&
StaggeredStokesPETScLevelSolver::getCouplingAwareASMSeedVelocityDOFs() const
{
    return d_coupling_aware_asm_seed_velocity_dofs;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

void
StaggeredStokesPETScLevelSolver::generateASMSubdomains(std::vector<std::set<int>>& overlap_is,
                                                       std::vector<std::set<int>>& nonoverlap_is)
{
    d_coupling_aware_asm_seed_velocity_dofs.clear();
    switch (d_asm_subdomain_construction_mode)
    {
    case ASMSubdomainConstructionMode::GEOMETRICAL:
        StaggeredStokesPETScMatUtilities::constructPatchLevelGeometricalASMSubdomains(overlap_is,
                                                                                      nonoverlap_is,
                                                                                      d_num_dofs_per_proc,
                                                                                      d_u_dof_index_idx,
                                                                                      d_level,
                                                                                      d_cf_boundary,
                                                                                      d_p_dof_index_idx,
                                                                                      d_box_size,
                                                                                      d_overlap_size);
        break;
    case ASMSubdomainConstructionMode::COUPLING_AWARE:
    {
        if (!d_petsc_mat)
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::generateASMSubdomains():\n"
                       << "  level matrix is not initialized for coupling-aware ASM subdomains.\n");
        }
        StaggeredStokesPETScMatUtilities::ensurePatchLevelCellClosureMapIsBuilt(
            d_coupling_aware_asm_map_data, d_u_dof_index_idx, d_p_dof_index_idx, d_level);
        StaggeredStokesPETScMatUtilities::computePatchLevelCouplingAwareASMSeedVelocityDofs(
            d_coupling_aware_asm_seed_velocity_dofs,
            d_u_dof_index_idx,
            d_level,
            d_coupling_aware_asm_map_data,
            d_coupling_aware_asm_seed_axis,
            d_coupling_aware_asm_seed_stride);
        if (d_coupling_aware_asm_closure_policy == CouplingAwareASMClosurePolicy::STRICT)
        {
            StaggeredStokesPETScMatUtilities::ensurePatchLevelVelocitySeedPairMapIsBuilt(
                d_coupling_aware_asm_map_data, d_u_dof_index_idx, d_level);
        }
        Mat A00_velocity_mat = nullptr;
        StaggeredStokesPETScMatUtilities::constructA00VelocitySubmatrix(
            A00_velocity_mat, d_petsc_mat, d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);
        StaggeredStokesPETScMatUtilities::constructPatchLevelCouplingAwareASMSubdomains(
            overlap_is,
            nonoverlap_is,
            d_num_dofs_per_proc,
            d_u_dof_index_idx,
            d_level,
            d_cf_boundary,
            A00_velocity_mat,
            d_coupling_aware_asm_map_data,
            d_coupling_aware_asm_seed_axis,
            d_coupling_aware_asm_seed_stride,
            d_coupling_aware_asm_closure_policy,
            d_coupling_aware_asm_relative_zero_tol);
        int ierr = MatDestroy(&A00_velocity_mat);
        IBTK_CHKERRQ(ierr);
        break;
    }
    case ASMSubdomainConstructionMode::UNKNOWN:
    default:
        TBOX_ERROR("StaggeredStokesPETScLevelSolver::generateASMSubdomains():\n"
                   << "  unsupported asm_subdomain_construction_mode = "
                   << enum_to_string(d_asm_subdomain_construction_mode) << ".\n");
    }

    if (d_log_asm_subdomains)
    {
        plog << d_object_name
             << "::generateASMSubdomains(): mode = " << enum_to_string(d_asm_subdomain_construction_mode) << "\n";
        std::size_t nonoverlap_sum = 0;
        std::size_t overlap_sum = 0;
        for (std::size_t k = 0; k < overlap_is.size(); ++k)
        {
            const std::size_t n_nonoverlap = k < nonoverlap_is.size() ? nonoverlap_is[k].size() : 0;
            const std::size_t n_overlap = overlap_is[k].size();
            nonoverlap_sum += n_nonoverlap;
            overlap_sum += n_overlap;
            plog << "  subdomain " << k << ": nonoverlap = " << n_nonoverlap << ", overlap = " << n_overlap
                 << ", delta = " << static_cast<long long>(n_overlap) - static_cast<long long>(n_nonoverlap) << "\n";
        }
        plog << "  totals: nonoverlap = " << nonoverlap_sum << ", overlap = " << overlap_sum
             << ", delta = " << static_cast<long long>(overlap_sum) - static_cast<long long>(nonoverlap_sum) << "\n";
    }

    return;
} // generateASMSubdomains

void
StaggeredStokesPETScLevelSolver::generateFieldSplitSubdomains(std::vector<std::string>& field_names,
                                                              std::vector<std::set<int>>& field_is)
{
    // Set IS'es for field split preconditioner.
    StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
        field_is, field_names, d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);

    return;
} // generateFieldSplitSubdomains

void
StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized(const SAMRAIVectorReal<NDIM, double>& x,
                                                                  const SAMRAIVectorReal<NDIM, double>& /*b*/)
{
    // Allocate DOF index data.
    if (!d_level->checkAllocated(d_u_dof_index_idx)) d_level->allocatePatchData(d_u_dof_index_idx);
    if (!d_level->checkAllocated(d_p_dof_index_idx)) d_level->allocatePatchData(d_p_dof_index_idx);
    StaggeredStokesPETScVecUtilities::constructPatchLevelDOFIndices(
        d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);

    // Setup PETSc objects.
    int ierr;
    const int mpi_rank = IBTK_MPI::getRank();
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_x);
    IBTK_CHKERRQ(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD, d_num_dofs_per_proc[mpi_rank], PETSC_DETERMINE, &d_petsc_b);
    IBTK_CHKERRQ(ierr);
    if (d_operator_mat)
    {
        ierr = MatDuplicate(d_operator_mat, MAT_COPY_VALUES, &d_petsc_mat);
        IBTK_CHKERRQ(ierr);
    }
    else
    {
        StaggeredStokesPETScMatUtilities::constructPatchLevelMACStokesOp(d_petsc_mat,
                                                                         d_U_problem_coefs,
                                                                         d_U_bc_coefs,
                                                                         d_new_time,
                                                                         d_num_dofs_per_proc,
                                                                         d_u_dof_index_idx,
                                                                         d_p_dof_index_idx,
                                                                         d_level);
    }

    std::vector<std::set<int>> field_is;
    std::vector<std::string> field_names;
    StaggeredStokesPETScMatUtilities::constructPatchLevelFields(
        field_is, field_names, d_num_dofs_per_proc, d_u_dof_index_idx, d_p_dof_index_idx, d_level);
    const auto get_field_dofs = [&field_is, &field_names](const std::string& field_name)
    {
        const auto field_name_it = std::find(field_names.begin(), field_names.end(), field_name);
        if (field_name_it == field_names.end())
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized():\n"
                       << "  unable to locate " << field_name << " field DOFs.\n");
        }
        const std::size_t field_idx = static_cast<std::size_t>(std::distance(field_names.begin(), field_name_it));
        return std::vector<PetscInt>(field_is[field_idx].begin(), field_is[field_idx].end());
    };
    d_velocity_dofs = get_field_dofs("velocity");
    d_pressure_dofs = get_field_dofs("pressure");
    d_velocity_dof_set.clear();
    d_pressure_dof_set.clear();
    d_velocity_dof_set.insert(d_velocity_dofs.begin(), d_velocity_dofs.end());
    d_pressure_dof_set.insert(d_pressure_dofs.begin(), d_pressure_dofs.end());
    d_custom_eigen_schur_subdomain_caches.clear();
    d_custom_eigen_a00_solver_storage.reset();

    if (d_augmented_operator_mat)
    {
        PetscInt full_m = 0, full_n = 0, aug_m = 0, aug_n = 0;
        ierr = MatGetSize(d_petsc_mat, &full_m, &full_n);
        IBTK_CHKERRQ(ierr);
        ierr = MatGetSize(d_augmented_operator_mat, &aug_m, &aug_n);
        IBTK_CHKERRQ(ierr);

        if (aug_m == full_m && aug_n == full_n)
        {
            ierr = MatAXPY(d_petsc_mat, 1.0, d_augmented_operator_mat, DIFFERENT_NONZERO_PATTERN);
            IBTK_CHKERRQ(ierr);
        }
        else
        {
            IS velocity_is = nullptr, velocity_is_all = nullptr;
            ierr = ISCreateGeneral(PETSC_COMM_WORLD,
                                   static_cast<PetscInt>(d_velocity_dofs.size()),
                                   d_velocity_dofs.empty() ? nullptr : d_velocity_dofs.data(),
                                   PETSC_COPY_VALUES,
                                   &velocity_is);
            IBTK_CHKERRQ(ierr);
            ierr = ISAllGather(velocity_is, &velocity_is_all);
            IBTK_CHKERRQ(ierr);

            PetscInt n_velocity_global = 0;
            ierr = ISGetSize(velocity_is_all, &n_velocity_global);
            IBTK_CHKERRQ(ierr);
            if (aug_m != n_velocity_global || aug_n != n_velocity_global)
            {
                ierr = ISDestroy(&velocity_is_all);
                IBTK_CHKERRQ(ierr);
                ierr = ISDestroy(&velocity_is);
                IBTK_CHKERRQ(ierr);
                TBOX_ERROR("StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized():\n"
                           << "  augmented operator has incompatible size: (" << aug_m << " x " << aug_n << ").\n"
                           << "  expected either full operator size (" << full_m << " x " << full_n
                           << ") or velocity block size (" << n_velocity_global << " x " << n_velocity_global
                           << ").\n");
            }

            const PetscInt* velocity_global_ids = nullptr;
            ierr = ISGetIndices(velocity_is_all, &velocity_global_ids);
            IBTK_CHKERRQ(ierr);

            PetscInt row_start = 0, row_end = 0;
            ierr = MatGetOwnershipRange(d_augmented_operator_mat, &row_start, &row_end);
            IBTK_CHKERRQ(ierr);
            ierr = MatSetOption(d_petsc_mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
            IBTK_CHKERRQ(ierr);
            std::vector<PetscInt> mapped_cols;
            for (PetscInt row = row_start; row < row_end; ++row)
            {
                PetscInt ncols = 0;
                const PetscInt* cols = nullptr;
                const PetscScalar* vals = nullptr;
                ierr = MatGetRow(d_augmented_operator_mat, row, &ncols, &cols, &vals);
                IBTK_CHKERRQ(ierr);
                if (row < 0 || row >= n_velocity_global)
                {
                    ierr = MatRestoreRow(d_augmented_operator_mat, row, &ncols, &cols, &vals);
                    IBTK_CHKERRQ(ierr);
                    TBOX_ERROR("StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized():\n"
                               << "  augmented velocity-block row index " << row << " is outside [0, "
                               << n_velocity_global - 1 << "].\n");
                }
                mapped_cols.resize(static_cast<std::size_t>(ncols));
                for (PetscInt k = 0; k < ncols; ++k)
                {
                    if (cols[k] < 0 || cols[k] >= n_velocity_global)
                    {
                        ierr = MatRestoreRow(d_augmented_operator_mat, row, &ncols, &cols, &vals);
                        IBTK_CHKERRQ(ierr);
                        TBOX_ERROR("StaggeredStokesPETScLevelSolver::initializeSolverStateSpecialized():\n"
                                   << "  augmented velocity-block column index " << cols[k] << " is outside [0, "
                                   << n_velocity_global - 1 << "] for row " << row << ".\n");
                    }
                    mapped_cols[static_cast<std::size_t>(k)] = velocity_global_ids[cols[k]];
                }
                const PetscInt full_row = velocity_global_ids[row];
                ierr = MatSetValues(d_petsc_mat, 1, &full_row, ncols, mapped_cols.data(), vals, ADD_VALUES);
                IBTK_CHKERRQ(ierr);
                ierr = MatRestoreRow(d_augmented_operator_mat, row, &ncols, &cols, &vals);
                IBTK_CHKERRQ(ierr);
            }
            ierr = MatAssemblyBegin(d_petsc_mat, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);
            ierr = MatAssemblyEnd(d_petsc_mat, MAT_FINAL_ASSEMBLY);
            IBTK_CHKERRQ(ierr);

            ierr = ISRestoreIndices(velocity_is_all, &velocity_global_ids);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&velocity_is_all);
            IBTK_CHKERRQ(ierr);
            ierr = ISDestroy(&velocity_is);
            IBTK_CHKERRQ(ierr);
        }
    }
    d_petsc_pc = d_petsc_mat;
    d_coupling_aware_asm_map_data.clear();
    d_coupling_aware_asm_seed_velocity_dofs.clear();
    // Set pressure nullspace if the level covers the entire domain.
    if (d_has_pressure_nullspace)
    {
        bool level_covers_entire_domain = d_level_num == 0;
        if (d_level_num > 0)
        {
            int local_cf_bdry_box_size = 0;
            for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
                const Array<BoundaryBox<NDIM>>& type_1_cf_bdry = d_cf_boundary->getBoundaries(patch->getPatchNumber(),
                                                                                              /* boundary type */ 1);
                local_cf_bdry_box_size += type_1_cf_bdry.size();
            }
            level_covers_entire_domain = IBTK_MPI::sumReduction(local_cf_bdry_box_size) == 0;
        }

        if (level_covers_entire_domain)
        {
            // Allocate pressure nullspace data.
            if (!d_level->checkAllocated(d_u_nullspace_idx)) d_level->allocatePatchData(d_u_nullspace_idx);
            if (!d_level->checkAllocated(d_p_nullspace_idx)) d_level->allocatePatchData(d_p_nullspace_idx);

            Pointer<SAMRAIVectorReal<NDIM, double>> nullspace_vec = new SAMRAIVectorReal<NDIM, double>(
                d_object_name + "nullspace_vec", d_hierarchy, d_level_num, d_level_num);
            nullspace_vec->addComponent(d_u_nullspace_var, d_u_nullspace_idx);
            nullspace_vec->addComponent(d_p_nullspace_var, d_p_nullspace_idx);
            for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
            {
                Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
                Pointer<SideData<NDIM, double>> u_patch_data = nullspace_vec->getComponentPatchData(0, *patch);
                u_patch_data->fill(0.0);
                Pointer<CellData<NDIM, double>> p_patch_data = nullspace_vec->getComponentPatchData(1, *patch);
                p_patch_data->fill(1.0);
            }

            LinearSolver::setNullSpace(
                /*const vec*/ false, std::vector<Pointer<SAMRAIVectorReal<NDIM, double>>>(1, nullspace_vec));
        }
    }

    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    d_data_synch_sched = StaggeredStokesPETScVecUtilities::constructDataSynchSchedule(u_idx, p_idx, d_level);
    d_ghost_fill_sched = StaggeredStokesPETScVecUtilities::constructGhostFillSchedule(u_idx, p_idx, d_level);
    return;
} // initializeSolverStateSpecialized

void
StaggeredStokesPETScLevelSolver::deallocateSolverStateSpecialized()
{
    // Deallocate DOF index data.
    if (d_level->checkAllocated(d_u_dof_index_idx)) d_level->deallocatePatchData(d_u_dof_index_idx);
    if (d_level->checkAllocated(d_p_dof_index_idx)) d_level->deallocatePatchData(d_p_dof_index_idx);
    d_coupling_aware_asm_map_data.clear();
    d_coupling_aware_asm_seed_velocity_dofs.clear();
    d_custom_eigen_schur_subdomain_caches.clear();
    d_custom_eigen_a00_solver_storage.reset();
    d_velocity_dof_set.clear();
    d_pressure_dof_set.clear();
    d_velocity_dofs.clear();
    d_pressure_dofs.clear();
    return;
} // deallocateSolverStateSpecialized

void
StaggeredStokesPETScLevelSolver::copyToPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level);
    return;
} // copyToPETScVec

void
StaggeredStokesPETScLevelSolver::copyFromPETScVec(Vec& petsc_x, SAMRAIVectorReal<NDIM, double>& x)
{
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    StaggeredStokesPETScVecUtilities::copyFromPatchLevelVec(
        petsc_x, u_idx, d_u_dof_index_idx, p_idx, d_p_dof_index_idx, d_level, d_data_synch_sched, d_ghost_fill_sched);
    return;
} // copyFromPETScVec

void
StaggeredStokesPETScLevelSolver::setupKSPVecs(Vec& petsc_x,
                                              Vec& petsc_b,
                                              SAMRAIVectorReal<NDIM, double>& x,
                                              SAMRAIVectorReal<NDIM, double>& b)
{
    if (d_initial_guess_nonzero) copyToPETScVec(petsc_x, x);
    const bool level_zero = (d_level_num == 0);
    const int u_idx = x.getComponentDescriptorIndex(0);
    const int p_idx = x.getComponentDescriptorIndex(1);
    const int f_idx = b.getComponentDescriptorIndex(0);
    const int h_idx = b.getComponentDescriptorIndex(1);
    const auto f_adj_idx = d_cached_eulerian_data.getCachedPatchDataIndex(f_idx);
    const auto h_adj_idx = d_cached_eulerian_data.getCachedPatchDataIndex(h_idx);
    for (PatchLevel<NDIM>::Iterator p(d_level); p; p++)
    {
        Pointer<Patch<NDIM>> patch = d_level->getPatch(p());
        Pointer<PatchGeometry<NDIM>> pgeom = patch->getPatchGeometry();
        Pointer<SideData<NDIM, double>> u_data = patch->getPatchData(u_idx);
        Pointer<SideData<NDIM, double>> f_data = patch->getPatchData(f_idx);
        Pointer<CellData<NDIM, double>> h_data = patch->getPatchData(h_idx);
        Pointer<SideData<NDIM, double>> f_adj_data = patch->getPatchData(f_adj_idx);
        Pointer<CellData<NDIM, double>> h_adj_data = patch->getPatchData(h_adj_idx);
        f_adj_data->copy(*f_data);
        h_adj_data->copy(*h_data);
        const bool at_physical_bdry = pgeom->intersectsPhysicalBoundary();
        // TODO: should we be using target data idx's here?
        StaggeredStokesPhysicalBoundaryHelper::setupBcCoefObjects(
            d_U_bc_coefs, d_P_bc_coef, u_idx, p_idx, d_homogeneous_bc);
        if (at_physical_bdry)
        {
            PoissonUtilities::adjustRHSAtPhysicalBoundary(
                *f_adj_data, patch, d_U_problem_coefs, d_U_bc_coefs, d_solution_time, d_homogeneous_bc);
            d_bc_helper->enforceNormalVelocityBoundaryConditions(
                f_adj_idx, h_adj_idx, d_U_bc_coefs, d_solution_time, d_homogeneous_bc, d_level_num, d_level_num);
        }
        const Array<BoundaryBox<NDIM>>& type_1_cf_bdry = level_zero ?
                                                             Array<BoundaryBox<NDIM>>() :
                                                             d_cf_boundary->getBoundaries(patch->getPatchNumber(),
                                                                                          /* boundary type */ 1);
        const bool at_cf_bdry = type_1_cf_bdry.size() > 0;
        if (at_cf_bdry)
        {
            PoissonUtilities::adjustRHSAtCoarseFineBoundary(
                *f_adj_data, *u_data, patch, d_U_problem_coefs, type_1_cf_bdry);
        }
    }

    StaggeredStokesPETScVecUtilities::copyToPatchLevelVec(
        petsc_b, f_adj_idx, d_u_dof_index_idx, h_adj_idx, d_p_dof_index_idx, d_level);

    copyToPETScVec(petsc_b, b);
    return;
} // setupKSPVecs

void
StaggeredStokesPETScLevelSolver::initializeEigenSubdomainSolver(const Eigen::MatrixXd& local_operator,
                                                                std::size_t subdomain_num)
{
    if (usingCustomEigenSolveSubdomainSolver())
    {
        std::vector<std::vector<int>>* overlap_subdomain_dofs = nullptr;
        getASMSubdomains(nullptr, &overlap_subdomain_dofs);
#if !defined(NDEBUG)
        TBOX_ASSERT(overlap_subdomain_dofs != nullptr);
        TBOX_ASSERT(subdomain_num < overlap_subdomain_dofs->size());
#endif
        if (d_custom_eigen_schur_subdomain_caches.size() != overlap_subdomain_dofs->size())
        {
            d_custom_eigen_schur_subdomain_caches.resize(overlap_subdomain_dofs->size());
        }
        if (!d_custom_eigen_a00_solver_storage)
        {
            dispatchEigenSolverType(d_custom_eigen_a00_solver_type,
                                    [this, n_subdomains = overlap_subdomain_dofs->size()](auto solver_tag)
                                    {
                                        using SolverType = typename decltype(solver_tag)::type;
                                        initializeCustomEigenA00SolveStorage<SolverType>(n_subdomains);
                                    });
        }

        const auto& overlap_dofs = (*overlap_subdomain_dofs)[subdomain_num];
        auto& cache = d_custom_eigen_schur_subdomain_caches[subdomain_num];
        cache = CustomEigenSchurSubdomainCache();
        cache.overlap_size = static_cast<int>(overlap_dofs.size());

        for (std::size_t local_pos = 0; local_pos < overlap_dofs.size(); ++local_pos)
        {
            const int dof = overlap_dofs[local_pos];
            if (d_velocity_dof_set.count(dof))
            {
                cache.velocity_positions.push_back(static_cast<int>(local_pos));
            }
            else if (d_pressure_dof_set.count(dof))
            {
                cache.pressure_positions.push_back(static_cast<int>(local_pos));
            }
            else
            {
                TBOX_ERROR("StaggeredStokesPETScLevelSolver::initializeEigenSubdomainSolver():\n"
                           << "  unable to classify local overlap DOF " << dof << " as velocity or pressure.\n");
            }
        }

        cache.A00 = extract_dense_block(local_operator, cache.velocity_positions, cache.velocity_positions);
        cache.A01 = extract_dense_block(local_operator, cache.velocity_positions, cache.pressure_positions);
        cache.A10 = extract_dense_block(local_operator, cache.pressure_positions, cache.velocity_positions);
        cache.A11 = extract_dense_block(local_operator, cache.pressure_positions, cache.pressure_positions);

        if (!cache.A00.size() && !cache.A11.size())
        {
            TBOX_ERROR("StaggeredStokesPETScLevelSolver::initializeEigenSubdomainSolver():\n"
                       << "  local custom Schur subdomain has no velocity or pressure DOFs.\n");
        }

        if (cache.A00.rows() > 0)
        {
            dispatchEigenSolverType(
                d_custom_eigen_a00_solver_type,
                [this, &cache, subdomain_num](auto solver_tag)
                {
                    using SolverType = typename decltype(solver_tag)::type;
                    auto& solver = getCustomEigenA00SolveStorage<SolverType>().solvers[subdomain_num];
                    initialize_custom_a00_solver(solver, cache.A00, d_custom_eigen_a00_solver_threshold, subdomain_num);
                });
        }

        if (cache.A11.rows() > 0)
        {
            if (cache.A00.rows() > 0)
            {
                dispatchEigenSolverType(d_custom_eigen_a00_solver_type,
                                        [this, &cache, subdomain_num](auto solver_tag)
                                        {
                                            using SolverType = typename decltype(solver_tag)::type;
                                            const auto& solver =
                                                getCustomEigenA00SolveStorage<SolverType>().solvers[subdomain_num];
                                            cache.A00_inv_A01 = solve_custom_a00(solver, cache.A01);
                                        });
                cache.schur = cache.A11 - cache.A10 * cache.A00_inv_A01;
            }
            else
            {
                cache.schur = cache.A11;
            }
            cache.schur_solve_matrix = buildCustomSchurSolveMatrix(cache.schur);
        }
        return;
    }

    PETScLevelSolver::initializeEigenSubdomainSolver(local_operator, subdomain_num);
    return;
} // initializeEigenSubdomainSolver

Eigen::VectorXd
StaggeredStokesPETScLevelSolver::solveEigenSubdomainSystem(const Eigen::VectorXd& rhs, std::size_t subdomain_num) const
{
    if (usingCustomEigenSolveSubdomainSolver())
    {
        const auto& cache = d_custom_eigen_schur_subdomain_caches[subdomain_num];
        Eigen::VectorXd solution = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(cache.overlap_size));

        const Eigen::VectorXd f0 = extract_subvector(rhs, cache.velocity_positions);
        const Eigen::VectorXd f1 = extract_subvector(rhs, cache.pressure_positions);

        if (cache.pressure_positions.empty())
        {
            Eigen::VectorXd u;
            dispatchEigenSolverType(d_custom_eigen_a00_solver_type,
                                    [this, subdomain_num, &f0, &u](auto solver_tag)
                                    {
                                        using SolverType = typename decltype(solver_tag)::type;
                                        const auto& solver =
                                            getCustomEigenA00SolveStorage<SolverType>().solvers[subdomain_num];
                                        u = solve_custom_a00(solver, f0).eval();
                                    });
            insert_subvector(solution, cache.velocity_positions, u);
            return solution;
        }

        if (cache.velocity_positions.empty())
        {
            const Eigen::VectorXd p = cache.schur_solve_matrix * f1;
            insert_subvector(solution, cache.pressure_positions, p);
            return solution;
        }

        Eigen::VectorXd y0;
        dispatchEigenSolverType(d_custom_eigen_a00_solver_type,
                                [this, subdomain_num, &f0, &y0](auto solver_tag)
                                {
                                    using SolverType = typename decltype(solver_tag)::type;
                                    const auto& solver =
                                        getCustomEigenA00SolveStorage<SolverType>().solvers[subdomain_num];
                                    y0 = solve_custom_a00(solver, f0).eval();
                                });
        const Eigen::VectorXd rhs_p = f1 - cache.A10 * y0;
        const Eigen::VectorXd p = cache.schur_solve_matrix * rhs_p;
        const Eigen::VectorXd u = y0 - cache.A00_inv_A01 * p;

        insert_subvector(solution, cache.velocity_positions, u);
        insert_subvector(solution, cache.pressure_positions, p);
        return solution;
    }

    return PETScLevelSolver::solveEigenSubdomainSystem(rhs, subdomain_num);
} // solveEigenSubdomainSystem

void
StaggeredStokesPETScLevelSolver::initializeCustomEigenShellData()
{
    const std::size_t n_subdomains = d_subdomain_dofs.size();
    if (d_custom_eigen_schur_subdomain_caches.size() != n_subdomains)
    {
        d_custom_eigen_schur_subdomain_caches.resize(n_subdomains);
    }

    auto initialize_impl = [this, n_subdomains](auto a00_solver_tag)
    {
        using SolverType = typename decltype(a00_solver_tag)::type;
        initializeCustomEigenA00SolveStorage<SolverType>(n_subdomains);
        initializeEigenShellDataWithLocalOperatorHook(
            [this](const Eigen::MatrixXd& local_operator, const std::size_t subdomain_num)
            {
                const auto& overlap_dofs = getEigenOverlapDofs(subdomain_num);
                auto& cache = d_custom_eigen_schur_subdomain_caches[subdomain_num];
                cache = CustomEigenSchurSubdomainCache();
                cache.overlap_size = static_cast<int>(overlap_dofs.size());
                cache.overlap_dofs = &overlap_dofs;

                for (std::size_t local_pos = 0; local_pos < overlap_dofs.size(); ++local_pos)
                {
                    const int dof = overlap_dofs[local_pos];
                    if (d_velocity_dof_set.count(dof))
                    {
                        cache.velocity_positions.push_back(static_cast<int>(local_pos));
                    }
                    else if (d_pressure_dof_set.count(dof))
                    {
                        cache.pressure_positions.push_back(static_cast<int>(local_pos));
                    }
                    else
                    {
                        TBOX_ERROR("StaggeredStokesPETScLevelSolver::initializeCustomEigenShellData():\n"
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
                    TBOX_ERROR("StaggeredStokesPETScLevelSolver::initializeCustomEigenShellData():\n"
                               << "  local custom Schur subdomain has no velocity or pressure DOFs.\n");
                }

                if (cache.A00.rows() > 0)
                {
                    auto& solver = getCustomEigenA00SolveStorage<SolverType>().solvers[subdomain_num];
                    initialize_custom_a00_solver(solver, cache.A00, d_custom_eigen_a00_solver_threshold, subdomain_num);
                }

                if (cache.A11.rows() > 0)
                {
                    if (cache.A00.rows() > 0)
                    {
                        const auto& solver = getCustomEigenA00SolveStorage<SolverType>().solvers[subdomain_num];
                        cache.A00_inv_A01 = solve_custom_a00(solver, cache.A01);
                        cache.schur = cache.A11 - cache.A10 * cache.A00_inv_A01;
                    }
                    else
                    {
                        cache.schur = cache.A11;
                    }
                    cache.schur_solve_matrix = buildCustomSchurSolveMatrix(cache.schur);
                }

                cache.rhs_workspace.resize(static_cast<Eigen::Index>(cache.overlap_size));
                cache.delta_workspace.resize(static_cast<Eigen::Index>(cache.overlap_size));
                cache.velocity_rhs_workspace.resize(static_cast<Eigen::Index>(cache.velocity_positions.size()));
                cache.pressure_rhs_workspace.resize(static_cast<Eigen::Index>(cache.pressure_positions.size()));
                cache.velocity_solution_workspace.resize(static_cast<Eigen::Index>(cache.velocity_positions.size()));
                cache.pressure_solution_workspace.resize(static_cast<Eigen::Index>(cache.pressure_positions.size()));
            });

        for (std::size_t subdomain_num = 0; subdomain_num < n_subdomains; ++subdomain_num)
        {
            auto& cache = d_custom_eigen_schur_subdomain_caches[subdomain_num];
            cache.update_dofs = &getEigenUpdateDofs(subdomain_num);
            cache.update_local_positions = &getEigenUpdateLocalPositions(subdomain_num);
            cache.active_residual_update_rows = &getEigenActiveResidualUpdateRows(subdomain_num);
            cache.active_residual_update_mat = &getEigenActiveResidualUpdateMat(subdomain_num);
            cache.residual_input_workspace.resize(static_cast<Eigen::Index>(cache.update_local_positions->size()));
            cache.residual_delta_workspace.resize(static_cast<Eigen::Index>(cache.active_residual_update_rows->size()));
        }
    };

    dispatchEigenSolverType(d_custom_eigen_a00_solver_type,
                            [&initialize_impl](auto solver_tag) { initialize_impl(solver_tag); });
    return;
} // initializeCustomEigenShellData

void
StaggeredStokesPETScLevelSolver::applyAdditiveCustomEigen(Vec x, Vec y)
{
    dispatchEigenSolverType(d_custom_eigen_a00_solver_type,
                            [this, x, y](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                applyAdditiveCustomEigenImpl<SolverType>(x, y);
                            });
    return;
} // applyAdditiveCustomEigen

void
StaggeredStokesPETScLevelSolver::applyMultiplicativeCustomEigen(Vec x, Vec y)
{
    dispatchEigenSolverType(d_custom_eigen_a00_solver_type,
                            [this, x, y](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                applyMultiplicativeCustomEigenImpl<SolverType>(x, y);
                            });
    return;
} // applyMultiplicativeCustomEigen

void
StaggeredStokesPETScLevelSolver::postprocessShellResult(Vec& y)
{
    const bool is_multiplicative_shell =
        d_shell_pc_type.compare(0, std::strlen("multiplicative"), "multiplicative") == 0;
    if (!is_multiplicative_shell || d_pressure_dofs.empty()) return;

    int ierr;
    PetscInt n_lo = 0, n_hi = 0;
    ierr = VecGetOwnershipRange(y, &n_lo, &n_hi);
    IBTK_CHKERRQ(ierr);

    PetscScalar* y_arr = nullptr;
    ierr = VecGetArray(y, &y_arr);
    IBTK_CHKERRQ(ierr);

    double local_sum = 0.0;
    int local_count = 0;
    for (const PetscInt dof : d_pressure_dofs)
    {
        if (dof < n_lo || dof >= n_hi) continue;
        const PetscInt local_idx = dof - n_lo;
        local_sum += PetscRealPart(y_arr[local_idx]);
        ++local_count;
    }
    const double global_sum = IBTK_MPI::sumReduction(local_sum);
    const int global_count = IBTK_MPI::sumReduction(local_count);
    if (global_count > 0)
    {
        const double pressure_mean = global_sum / static_cast<double>(global_count);
        for (const PetscInt dof : d_pressure_dofs)
        {
            if (dof < n_lo || dof >= n_hi) continue;
            const PetscInt local_idx = dof - n_lo;
            y_arr[local_idx] -= pressure_mean;
        }
    }
    ierr = VecRestoreArray(y, &y_arr);
    IBTK_CHKERRQ(ierr);
    return;
} // postprocessShellResult

void
StaggeredStokesPETScLevelSolver::setOperatorMat(Mat operator_mat)
{
    d_operator_mat = operator_mat;
    return;
} // setOperatorMat

void
StaggeredStokesPETScLevelSolver::setAugmentedOperatorMat(Mat augmented_operator_mat)
{
    d_augmented_operator_mat = augmented_operator_mat;
    return;
} // setAugmentedOperatorMat

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
