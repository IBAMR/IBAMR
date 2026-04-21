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

#include <ibtk/IBTK_MPI.h>
#include <ibtk/PETScLevelSolver.h>
#include <ibtk/private/PETScLevelSolverEigenReferenceShellBackend.h>

#include <tbox/Database.h>

namespace IBTK
{
namespace
{
std::unique_ptr<PETScLevelSolverShellBackend>
allocate_eigen_reference_shell_backend(PETScLevelSolver& /*solver*/)
{
    return std::make_unique<PETScLevelSolverEigenReferenceShellBackend>();
}

class PETScLevelSolverEigenReferenceShellBackendRegistrar
{
public:
    PETScLevelSolverEigenReferenceShellBackendRegistrar()
    {
        PETScLevelSolverShellBackendManager::getManager()->registerShellBackendFactoryFunction(
            "eigen-reference", allocate_eigen_reference_shell_backend);
    }
};

static PETScLevelSolverEigenReferenceShellBackendRegistrar eigen_reference_shell_backend_registrar;
} // namespace

const std::string PETScLevelSolverEigenReferenceShellBackend::s_backend_name = "EigenReference";

void
PETScLevelSolverEigenReferenceShellBackend::configureFromInputDatabase(
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    d_solver_type = EigenSubdomainSolverType::PARTIAL_PIV_LU;
    d_solver_threshold = -1.0;
    if (!input_db) return;

    if (input_db->keyExists("eigen_subdomain_solver_type"))
    {
        d_solver_type = parseEigenSubdomainSolverType(
            input_db->getString("eigen_subdomain_solver_type"),
            "PETScLevelSolverEigenReferenceShellBackend::PETScLevelSolverEigenReferenceShellBackend()");
    }
    if (input_db->keyExists("eigen_subdomain_solver_threshold"))
    {
        d_solver_threshold = input_db->getDouble("eigen_subdomain_solver_threshold");
    }
}

const std::string&
PETScLevelSolverEigenReferenceShellBackend::getName() const
{
    return s_backend_name;
}

void
PETScLevelSolverEigenReferenceShellBackend::initializeLocalSubdomainSolver(const Eigen::MatrixXd& local_operator,
                                                                           const std::size_t subdomain_num)
{
    dispatchEigenSolverType(getSolverType(),
                            [this, &local_operator, subdomain_num](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                initializeTypedSubdomainSolver<SolverType>(local_operator, subdomain_num);
                            });
}

void
PETScLevelSolverEigenReferenceShellBackend::initializeAdditionalSolverState()
{
    d_level_mat = copyPETScMatToEigenSparse(getSolverState().petsc_mat);
    dispatchEigenSolverType(getSolverType(),
                            [this](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                initializeSolveStorage<SolverType>();
                            });
}

void
PETScLevelSolverEigenReferenceShellBackend::deallocateAdditionalSolverState()
{
    d_solve_storage.reset();
    d_level_mat.resize(0, 0);
}

Eigen::VectorXd
PETScLevelSolverEigenReferenceShellBackend::solveLocalSubdomainSystem(const Eigen::VectorXd& rhs,
                                                                      const std::size_t subdomain_num) const
{
    Eigen::VectorXd solution;
    dispatchEigenSolverType(getSolverType(),
                            [this, &rhs, subdomain_num, &solution](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                solution = solveTypedSubdomainSystem<SolverType>(rhs, subdomain_num);
                            });
    return solution;
}

void
PETScLevelSolverEigenReferenceShellBackend::apply(Vec x, Vec y)
{
    if (!d_solver_state.use_multiplicative)
    {
        TBOX_ERROR("PETScLevelSolverEigenReferenceShellBackend::apply():\n"
                   << "  shell backend 'eigen-reference' only supports multiplicative mode.\n"
                   << "  choose shell_pc_type = multiplicative-eigen-reference (-basic/-restrict optional).\n");
    }
    TBOX_ASSERT(IBTK_MPI::getNodes() == 1);
    auto& common_subdomains = getCommonSubdomains();
    const Eigen::Index n = getNumDofs();
    TBOX_ASSERT(n > 0);

    {
        ConstPetscVecArrayMap x_array(x, n);
        PetscVecArrayMap y_array(y, n);
        const auto x_map = x_array.getMap();
        auto y_map = y_array.getMap();
        Eigen::VectorXd residual = Eigen::VectorXd::Zero(n);
        y_map.setZero();
        for (std::size_t subdomain_num = 0; subdomain_num < common_subdomains.size(); ++subdomain_num)
        {
            auto& cache = common_subdomains[subdomain_num];

            // Reference-style multiplicative sweep: recompute residual each subdomain.
            residual = x_map - d_level_mat * y_map;

            std::size_t rhs_idx = 0;
            for (const int dof : cache.overlap_dofs)
            {
                cache.rhs_workspace[static_cast<Eigen::Index>(rhs_idx++)] = residual[dof];
            }

            cache.delta_workspace = solveLocalSubdomainSystem(cache.rhs_workspace, subdomain_num);

            std::size_t update_idx = 0;
            for (const int dof : cache.update_dofs)
            {
                y_map[dof] +=
                    cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
            }
        }
    }
    d_solver_state.postprocess_result(y);
}

PETScLevelSolverEigenReferenceShellBackend::EigenSubdomainSolverType
PETScLevelSolverEigenReferenceShellBackend::getSolverType() const
{
    return d_solver_type;
}

double
PETScLevelSolverEigenReferenceShellBackend::getSolverThreshold() const
{
    return d_solver_threshold;
}
} // namespace IBTK
