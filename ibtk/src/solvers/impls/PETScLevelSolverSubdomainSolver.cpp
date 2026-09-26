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

#include <ibtk/PETScLevelSolverSubdomainSolver.h>

#include <tbox/Utilities.h>

namespace IBTK
{
PETScLevelSolverSubdomainSolver::PETScLevelSolverSubdomainSolver(PETScLevelSolverSubdomainSolver&& other) noexcept =
    default;

PETScLevelSolverSubdomainSolver&
PETScLevelSolverSubdomainSolver::operator=(PETScLevelSolverSubdomainSolver&& other) noexcept = default;

PETScLevelSolverSubdomainSolver::~PETScLevelSolverSubdomainSolver() = default;

PETScLevelSolverSubdomainSolver::operator bool() const
{
    return d_adapter != nullptr;
}

void
PETScLevelSolverSubdomainSolver::initializeSolverState(const std::vector<Mat>& matrices,
                                                       const std::vector<IS>& subdomains,
                                                       const std::string& options_prefix)
{
    TBOX_ASSERT(d_adapter);
    d_adapter->initializeSolverState(matrices, subdomains, options_prefix);
}

void
PETScLevelSolverSubdomainSolver::deallocateSolverState()
{
    TBOX_ASSERT(d_adapter);
    d_adapter->deallocateSolverState();
}

void
PETScLevelSolverSubdomainSolver::solve(const std::size_t first, const std::size_t last, Vec b, Vec x)
{
    TBOX_ASSERT(d_adapter);
    d_adapter->solve(first, last, b, x);
}
} // namespace IBTK
