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

#ifndef included_IBTK_PETScLevelSolverSubdomainSolver_inl
#define included_IBTK_PETScLevelSolverSubdomainSolver_inl

#include <ibtk/config.h>

#include <ibtk/PETScLevelSolverSubdomainSolver.h>

namespace IBTK
{
template <class Implementation>
class PETScLevelSolverSubdomainSolver::ImplementationAdapter final : public Adapter
{
public:
    template <class... Args>
    explicit ImplementationAdapter(Args&&... args) : d_implementation(std::forward<Args>(args)...)
    {
    }

    void initializeSolverState(const std::vector<Mat>& matrices,
                               const std::vector<IS>& subdomains,
                               const std::string& options_prefix) override
    {
        d_implementation.initializeSolverState(matrices, subdomains, options_prefix);
    }

    void deallocateSolverState() override
    {
        d_implementation.deallocateSolverState();
    }

    void solve(const std::size_t first, const std::size_t last, Vec b, Vec x) override
    {
        d_implementation.solve(first, last, b, x);
    }

private:
    Implementation d_implementation;
};

template <class Implementation, class... Args>
requires PETScLevelSolverSubdomainSolverImplementation<Implementation>&&
    std::constructible_from<Implementation, Args...>
    PETScLevelSolverSubdomainSolver::PETScLevelSolverSubdomainSolver(std::in_place_type_t<Implementation>,
                                                                     Args&&... args)
    : d_adapter(std::make_unique<ImplementationAdapter<Implementation>>(std::forward<Args>(args)...))
{
}
} // namespace IBTK

#endif
