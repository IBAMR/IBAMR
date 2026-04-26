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

#include <ibtk/PETScLevelSolver.h>
#include <ibtk/private/PETScLevelSolverEigenShellBackend.h>

namespace IBTK
{
namespace
{
std::unique_ptr<PETScLevelSolverShellBackend>
allocate_eigen_shell_backend(PETScLevelSolver& /*solver*/)
{
    return std::make_unique<PETScLevelSolverEigenShellBackend>();
}

class PETScLevelSolverEigenShellBackendRegistrar
{
public:
    PETScLevelSolverEigenShellBackendRegistrar()
    {
        PETScLevelSolverShellBackendManager::getManager()->registerShellBackendFactoryFunction(
            "eigen", allocate_eigen_shell_backend);
    }
};

static PETScLevelSolverEigenShellBackendRegistrar eigen_shell_backend_registrar;
} // namespace

const std::string PETScLevelSolverEigenShellBackend::s_backend_name = "Eigen";

const std::string&
PETScLevelSolverEigenShellBackend::getName() const
{
    return s_backend_name;
}
} // namespace IBTK
