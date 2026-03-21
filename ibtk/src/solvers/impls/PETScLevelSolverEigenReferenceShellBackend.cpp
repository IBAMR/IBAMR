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
PETScLevelSolverEigenReferenceShellBackend::PETScLevelSolverEigenReferenceShellBackend(PETScLevelSolver& solver)
    : PETScLevelSolverEigenShellBackendBase(solver)
{
}

const std::string&
PETScLevelSolverEigenReferenceShellBackend::getTypeKey() const
{
    return d_type_key;
}

const char*
PETScLevelSolverEigenReferenceShellBackend::getPCNameSuffixAdditive() const
{
    return "PC_AdditiveEigenReference";
}

const char*
PETScLevelSolverEigenReferenceShellBackend::getPCNameSuffixMultiplicative() const
{
    return "PC_MultiplicativeEigenReference";
}

void
PETScLevelSolverEigenReferenceShellBackend::configure(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db)
{
    d_solver_type = EigenSubdomainSolverType::PARTIAL_PIV_LU;
    d_solver_threshold = -1.0;
    if (!input_db) return;

    if (input_db->keyExists("eigen_subdomain_solver_type"))
    {
        d_solver_type = parseEigenSubdomainSolverType(input_db->getString("eigen_subdomain_solver_type"),
                                                      "PETScLevelSolverEigenReferenceShellBackend::configure()");
    }
    if (input_db->keyExists("eigen_subdomain_solver_threshold"))
    {
        d_solver_threshold = input_db->getDouble("eigen_subdomain_solver_threshold");
    }
}

void
PETScLevelSolverEigenReferenceShellBackend::initialize()
{
    d_level_mat = copyPETScMatToEigenSparse(d_context.getPETScMatForBackend());
    dispatchEigenSolverType(getSolverType(),
                            [this](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                initializeSolveStorage<SolverType>();
                            });
}

void
PETScLevelSolverEigenReferenceShellBackend::deallocate()
{
    clearCommonData();
    d_solve_storage.reset();
    d_level_mat.resize(0, 0);
}

void
PETScLevelSolverEigenReferenceShellBackend::initializeSubdomainSolver(const Eigen::MatrixXd& local_operator,
                                                                      const std::size_t subdomain_num)
{
    dispatchEigenSolverType(getSolverType(),
                            [this, &local_operator, subdomain_num](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                initializeSubdomainSolver<SolverType>(local_operator, subdomain_num);
                            });
}

Eigen::VectorXd
PETScLevelSolverEigenReferenceShellBackend::solveSubdomainSystem(const Eigen::VectorXd& rhs,
                                                                 const std::size_t subdomain_num) const
{
    Eigen::VectorXd solution;
    dispatchEigenSolverType(getSolverType(),
                            [this, &rhs, subdomain_num, &solution](auto solver_tag)
                            {
                                using SolverType = typename decltype(solver_tag)::type;
                                solution = solveSubdomainSystem<SolverType>(rhs, subdomain_num);
                            });
    return solution;
}

void
PETScLevelSolverEigenReferenceShellBackend::applyAdditive(Vec /*x*/, Vec /*y*/)
{
    TBOX_ERROR("PETScLevelSolverEigenReferenceShellBackend::applyAdditive():\n"
               << "  shell backend 'eigen-reference' only supports multiplicative mode.\n"
               << "  choose shell_pc_type = multiplicative-eigen-reference (-basic/-restrict optional).\n");
}

void
PETScLevelSolverEigenReferenceShellBackend::applyMultiplicative(Vec x, Vec y)
{
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

            cache.delta_workspace = solveSubdomainSystem(cache.rhs_workspace, subdomain_num);

            std::size_t update_idx = 0;
            for (const int dof : cache.update_dofs)
            {
                y_map[dof] +=
                    cache.delta_workspace[static_cast<Eigen::Index>(cache.update_local_positions[update_idx++])];
            }
        }
    }
    d_context.postprocessShellResultForBackend(y);
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
