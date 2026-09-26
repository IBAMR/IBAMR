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

#ifndef included_IBTK_DOFCoverage
#define included_IBTK_DOFCoverage

#include <ibtk/config.h>

#include <petscis.h>

#include <set>
#include <string>
#include <vector>

namespace IBTK
{
/*!
 * \brief How often each DOF must occur in a collection of subdomains.
 */
enum class DOFCoverage
{
    //! Each DOF occurs in at least one subdomain.
    AT_LEAST_ONCE,
    //! Each DOF occurs in exactly one subdomain.
    EXACTLY_ONCE
};

/*!
 * \brief Whether solvers check the coverage of their subdomains, when the input database does not say. The
 * default is to check in debug builds only: a solver rebuilds its subdomains, and so would repeat this check,
 * every time its state is initialized, which happens whenever the operator's boundary configuration changes, so
 * a collective check that touches every DOF on every rank is not paid on every rebuild in a release build unless
 * requested.
 */
bool default_check_dof_coverage();

/*!
 * \brief Check that a collection of subdomains covers the DOFs as required.
 *
 * Only some solvers and preconditioners need their subdomains to cover the DOFs, so this check is for the
 * implementations that do, and they should run it only when requested (see default_check_dof_coverage()). The
 * DOFs are the global indices 0, ..., N - 1 of a vector whose numbers of entries on the ranks are n_local_dofs,
 * and every rank passes its own subdomains. A subdomain may contain DOFs of other ranks. The check is
 * collective; it counts the occurrences of each DOF in the subdomains of all ranks and, if the requirement
 * fails, reports on every rank the number of DOFs that do not satisfy it and the lowest such DOF.
 *
 * \param context Text that starts the error message, such as the name of the solver and the method.
 */
void check_dof_coverage(const std::string& context,
                        const std::vector<IS>& subdomains,
                        PetscInt n_local_dofs,
                        DOFCoverage coverage);

/*!
 * \brief Check a collection of subdomains given as sets of DOFs; see the overload for index sets.
 */
void check_dof_coverage(const std::string& context,
                        const std::vector<std::set<int>>& subdomains,
                        PetscInt n_local_dofs,
                        DOFCoverage coverage);
} // namespace IBTK

#endif
