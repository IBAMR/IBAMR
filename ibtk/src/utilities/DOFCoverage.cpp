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

#include <ibtk/DOFCoverage.h>
#include <ibtk/IBTK_CHKERRQ.h>
#include <ibtk/IBTK_MPI.h>

#include <tbox/Utilities.h>

#include <petscvec.h>

#include <algorithm>
#include <limits>
#include <string>

#include <ibtk/namespaces.h> // IWYU pragma: keep

namespace IBTK
{
namespace
{
// Count the occurrences of the DOFs in all subdomains, given a function that calls its argument with the DOFs of
// each subdomain, and check them.
template <class ForEachSubdomain>
void
check_counts(const std::string& context,
             const PetscInt n_local_dofs,
             const DOFCoverage coverage,
             ForEachSubdomain&& for_each_subdomain)
{
    Vec counts = nullptr;
    int ierr = VecCreateMPI(PETSC_COMM_WORLD, n_local_dofs, PETSC_DETERMINE, &counts);
    IBTK_CHKERRQ(ierr);
    ierr = VecSet(counts, 0.0);
    IBTK_CHKERRQ(ierr);
    std::vector<PetscScalar> ones;
    for_each_subdomain(
        [&](const PetscInt* dofs, const PetscInt n)
        {
            ones.resize(std::max<std::size_t>(ones.size(), n), 1.0);
            const int set_ierr = VecSetValues(counts, n, dofs, ones.data(), ADD_VALUES);
            IBTK_CHKERRQ(set_ierr);
        });
    ierr = VecAssemblyBegin(counts);
    IBTK_CHKERRQ(ierr);
    ierr = VecAssemblyEnd(counts);
    IBTK_CHKERRQ(ierr);

    PetscInt first = 0, last = 0;
    ierr = VecGetOwnershipRange(counts, &first, &last);
    IBTK_CHKERRQ(ierr);
    const PetscScalar* values = nullptr;
    ierr = VecGetArrayRead(counts, &values);
    IBTK_CHKERRQ(ierr);
    int n_missing = 0, n_repeated = 0;
    int first_bad = std::numeric_limits<int>::max();
    for (PetscInt i = 0; i < last - first; ++i)
    {
        const double count = PetscRealPart(values[i]);
        const bool missing = count < 0.5;
        const bool repeated = coverage == DOFCoverage::EXACTLY_ONCE && count > 1.5;
        n_missing += missing;
        n_repeated += repeated;
        if ((missing || repeated) && first_bad == std::numeric_limits<int>::max())
        {
            first_bad = static_cast<int>(first + i);
        }
    }
    ierr = VecRestoreArrayRead(counts, &values);
    IBTK_CHKERRQ(ierr);
    ierr = VecDestroy(&counts);
    IBTK_CHKERRQ(ierr);
    n_missing = IBTK_MPI::sumReduction(n_missing);
    n_repeated = IBTK_MPI::sumReduction(n_repeated);
    first_bad = IBTK_MPI::minReduction(first_bad);
    if (n_missing + n_repeated > 0)
    {
        const auto describe = [](const int n, const std::string& where)
        { return std::to_string(n) + (n == 1 ? " DOF is " : " DOFs are ") + where; };
        std::string failures;
        if (n_missing > 0) failures = describe(n_missing, "in no subdomain");
        if (n_repeated > 0)
        {
            failures += (failures.empty() ? "" : " and ") + describe(n_repeated, "in more than one subdomain");
        }
        TBOX_ERROR(context << ":\n"
                           << "  the subdomains do not cover the DOFs as required: " << failures
                           << "; the lowest DOF that fails is " << first_bad << ".\n");
    }
}
} // namespace

bool
default_check_dof_coverage()
{
#if defined(NDEBUG)
    return false;
#else
    return true;
#endif
}

void
check_dof_coverage(const std::string& context,
                   const std::vector<IS>& subdomains,
                   const PetscInt n_local_dofs,
                   const DOFCoverage coverage)
{
    check_counts(context,
                 n_local_dofs,
                 coverage,
                 [&](auto&& count)
                 {
                     for (IS subdomain : subdomains)
                     {
                         PetscInt n = 0;
                         const PetscInt* dofs = nullptr;
                         int ierr = ISGetLocalSize(subdomain, &n);
                         IBTK_CHKERRQ(ierr);
                         ierr = ISGetIndices(subdomain, &dofs);
                         IBTK_CHKERRQ(ierr);
                         count(dofs, n);
                         ierr = ISRestoreIndices(subdomain, &dofs);
                         IBTK_CHKERRQ(ierr);
                     }
                 });
}

void
check_dof_coverage(const std::string& context,
                   const std::vector<std::set<int>>& subdomains,
                   const PetscInt n_local_dofs,
                   const DOFCoverage coverage)
{
    check_counts(context,
                 n_local_dofs,
                 coverage,
                 [&](auto&& count)
                 {
                     std::vector<PetscInt> dofs;
                     for (const std::set<int>& subdomain : subdomains)
                     {
                         dofs.assign(subdomain.begin(), subdomain.end());
                         count(dofs.data(), static_cast<PetscInt>(dofs.size()));
                     }
                 });
}
} // namespace IBTK
