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
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_CHKERRQ.h>

#include <tbox/Logger.h>
#include <tbox/PIO.h>

#include <mpi.h>

#include <set>
#include <string>
#include <vector>

#include "../tests.h"

using namespace IBTK;

// The scenario is the part of the input file name that follows "dof_coverage.". Each rank owns 6 DOFs.
int
main(int argc, char* argv[])
{
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);
    SAMRAI::tbox::Pointer<SAMRAI::tbox::Logger::Appender> abort_appender = new TestAppender();
    SAMRAI::tbox::Logger::getInstance()->setAbortAppender(abort_appender);
    SAMRAI::tbox::PIO::logOnlyNodeZero("output");

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const std::string path = argc > 1 ? argv[1] : "";
    const std::string name = path.substr(path.rfind('/') + 1);
    std::string scenario = name.substr(name.find('.') + 1, name.find('.', name.find('.') + 1) - name.find('.') - 1);
    // A scenario with the suffix "_is" checks the index sets only, so that a violation is reported by that overload.
    const bool index_sets_only = scenario.size() > 3 && scenario.compare(scenario.size() - 3, 3, "_is") == 0;
    if (index_sets_only) scenario.resize(scenario.size() - 3);
    const int n = 6;
    auto own = [&](const int owner, const int first, const int last)
    {
        std::set<int> dofs;
        for (int i = first; i < last; ++i) dofs.insert(n * owner + i);
        return dofs;
    };
    const int next = (rank + 1) % size;

    std::vector<std::set<int>> subdomains;
    DOFCoverage coverage = DOFCoverage::EXACTLY_ONCE;
    if (scenario == "exact")
    {
        subdomains = { own(rank, 0, 3), own(rank, 3, 6) };
    }
    else if (scenario == "cover_by_overlap")
    {
        coverage = DOFCoverage::AT_LEAST_ONCE;
        subdomains = { own(rank, 0, 4), own(rank, 2, 6) };
    }
    else if (scenario == "repeated")
    {
        subdomains = { own(rank, 0, 4), own(rank, 2, 6) };
    }
    else if (scenario == "gap")
    {
        coverage = DOFCoverage::AT_LEAST_ONCE;
        subdomains = { own(rank, 0, 3) };
    }
    else if (scenario == "remote")
    {
        // Each rank's subdomain holds the DOFs of the next rank.
        subdomains = { own(next, 0, 6) };
    }
    else if (scenario == "remote_gap")
    {
        subdomains = { own(next, 0, rank == 0 ? 6 : 5) };
    }
    else
    {
        TBOX_ERROR("unknown scenario " << scenario << "\n");
    }

    std::vector<IS> index_sets;
    for (const std::set<int>& subdomain : subdomains)
    {
        std::vector<PetscInt> dofs(subdomain.begin(), subdomain.end());
        IS is = nullptr;
        int ierr =
            ISCreateGeneral(PETSC_COMM_SELF, static_cast<PetscInt>(dofs.size()), dofs.data(), PETSC_COPY_VALUES, &is);
        IBTK_CHKERRQ(ierr);
        index_sets.push_back(is);
    }
    if (!index_sets_only) check_dof_coverage("check_sets", subdomains, n, coverage);
    check_dof_coverage("check_index_sets", index_sets, n, coverage);
    for (IS& is : index_sets)
    {
        const int ierr = ISDestroy(&is);
        IBTK_CHKERRQ(ierr);
    }
    SAMRAI::tbox::plog << scenario << ": covered as required\n";
    return 0;
}
