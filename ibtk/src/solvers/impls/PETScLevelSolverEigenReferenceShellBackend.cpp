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

} // namespace IBTK
