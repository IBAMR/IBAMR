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

const std::string&
PETScLevelSolverEigenReferenceShellBackend::getName() const
{
    return s_backend_name;
}

void
PETScLevelSolverEigenReferenceShellBackend::initializeAdditionalSolverState()
{
    d_level_mat = copyPETScMatToEigenSparse(getSolverState().petsc_mat);
    PETScLevelSolverEigenFactorizedLocalSolveShellBackend::initializeAdditionalSolverState();
    return;
}

void
PETScLevelSolverEigenReferenceShellBackend::deallocateAdditionalSolverState()
{
    PETScLevelSolverEigenFactorizedLocalSolveShellBackend::deallocateAdditionalSolverState();
    d_level_mat.resize(0, 0);
    return;
}

void
PETScLevelSolverEigenReferenceShellBackend::initializeSubdomainSweep(Vec x, Vec y)
{
    if (!d_solver_state.use_multiplicative)
    {
        TBOX_ERROR("PETScLevelSolverEigenReferenceShellBackend::apply():\n"
                   << "  shell backend 'eigen-reference' only supports multiplicative mode.\n"
                   << "  choose shell_pc_type = multiplicative-eigen-reference (-basic/-restrict optional).\n");
    }
    PETScLevelSolverEigenShellBackendBase::initializeSubdomainSweep(x, y);
}

void
PETScLevelSolverEigenReferenceShellBackend::beginSubdomainRhs(const std::size_t subdomain_num, Vec x, Vec y)
{
    const auto x_map =
        Eigen::Map<const Eigen::VectorXd>(reinterpret_cast<const double*>(d_sweep_x_values), getNumDofs());
    const auto y_map =
        Eigen::Map<const Eigen::VectorXd>(reinterpret_cast<const double*>(d_sweep_y_values), getNumDofs());
    d_sweep_residual = x_map - d_level_mat * y_map;
    PETScLevelSolverEigenShellBackendBase::beginSubdomainRhs(subdomain_num, x, y);
}

void
PETScLevelSolverEigenReferenceShellBackend::updateSubdomainResidual(std::size_t /*subdomain_num*/, Vec /*x*/, Vec /*y*/)
{
    // The reference backend recomputes the original residual before each patch.
}

} // namespace IBTK
