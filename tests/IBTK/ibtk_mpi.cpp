// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

// Config files

// Set up application namespace declarations
#include <ibtk/IBTKInit.h>
#include <ibtk/IBTK_MPI.h>

#include <fstream>
#include <limits>

#include <ibtk/app_namespaces.h>

template <typename T>
bool minReduction(T x, T exact, int proc);

template <typename T>
bool maxReduction(T x, T exact, int proc);

template <typename T>
bool sumReduction(T x, T exact);

template <typename T>
bool allToOneSumReduction(T x, T exact, int root);

template <typename T>
bool bcast(T x, T exact, int root);

template <typename T>
bool sendAndRecv(T x, T exact);

template <typename T>
bool allGather(T x);

/*******************************************************************************
 * For each run, the input filename must be given on the command line.  In all *
 * cases, the command line is:                                                 *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize IBTK
    IBTK::IBTKInit ibtkInit(argc, argv, MPI_COMM_WORLD);

    int num_nodes = IBTK::IBTK_MPI::getNodes();
    int rank = IBTK::IBTK_MPI::getRank();

    std::ofstream output_file;
    if (!rank) output_file.open("output");

    // integers
    int x = IBTK::IBTK_MPI::getRank();
    bool passed = minReduction(x, 0, 0);

    passed = IBTK::IBTK_MPI::maxReduction(passed ? 1 : 0);
    if (!rank) output_file << "Min test " << (passed ? "passed" : "failed") << ".\n";

    passed = maxReduction(x, IBTK::IBTK_MPI::getNodes() - 1, IBTK::IBTK_MPI::getNodes() - 1);

    passed = IBTK::IBTK_MPI::maxReduction(passed ? 1 : 0);
    if (!rank) output_file << "Max test " << (passed ? "passed" : "failed") << ".\n";

    passed = sumReduction(x, num_nodes * (num_nodes - 1) / 2);

    passed = IBTK::IBTK_MPI::maxReduction(passed ? 1 : 0);
    if (!rank) output_file << "sumReduction test " << (passed ? "passed" : "failed") << ".\n";

    passed = allToOneSumReduction(x, num_nodes * (num_nodes - 1) / 2, 0);

    passed = IBTK::IBTK_MPI::maxReduction(passed ? 1 : 0);
    if (!rank) output_file << "allToOneSum test " << (passed ? "passed" : "failed") << ".\n";

    passed = bcast(x, num_nodes - 1, num_nodes - 1);

    passed = IBTK::IBTK_MPI::maxReduction(passed ? 1 : 0);
    if (!rank) output_file << "bcast test " << (passed ? "passed" : "failed") << ".\n";

    passed = sendAndRecv(x, (x - 1 + num_nodes) % num_nodes);

    passed = IBTK::IBTK_MPI::maxReduction(passed ? 1 : 0);
    if (!rank) output_file << "send and recv test " << (passed ? "passed" : "failed") << ".\n";
    passed = allGather(x);

    passed = IBTK::IBTK_MPI::maxReduction(passed ? 1 : 0);
    if (!rank) output_file << "all gather test " << (passed ? "passed" : "failed") << ".\n";

    double x_min = rank, x_max = rank;
    IBTK::IBTK_MPI::minMaxReduction(&x_min, &x_max);
    if (!rank) output_file << "minMaxReduction: " << x_min << " " << x_max << "\n";

    // Integer reductions must be correct for the extreme values of the type.
    // In the last entry the largest value is small and one rank holds the
    // lowest value.
    const int i_lowest = std::numeric_limits<int>::lowest(), i_highest = std::numeric_limits<int>::max();
    int i_min[3] = { rank == num_nodes - 1 ? i_lowest : rank, i_highest, rank };
    int i_max[3] = { rank == 0 ? i_highest : -rank, i_lowest, rank == 0 ? i_lowest : -rank };
    IBTK::IBTK_MPI::minMaxReduction(i_min, i_max, 3);
    if (!rank)
    {
        output_file << "minMaxReduction of ints: " << i_min[0] << " " << i_min[1] << " " << i_min[2] << " " << i_max[0]
                    << " " << i_max[1] << " " << i_max[2] << "\n";
    }

    const unsigned int u_highest = std::numeric_limits<unsigned int>::max();
    unsigned int u_min[2] = { static_cast<unsigned int>(rank) + 1u, u_highest };
    unsigned int u_max[2] = { rank == 0 ? u_highest : static_cast<unsigned int>(rank),
                              static_cast<unsigned int>(rank) };
    IBTK::IBTK_MPI::minMaxReduction(u_min, u_max, 2);
    if (!rank)
    {
        output_file << "minMaxReduction of unsigned ints: " << u_min[0] << " " << u_min[1] << " " << u_max[0] << " "
                    << u_max[1] << "\n";
    }

    int a[2] = { rank, 1 }, b = 2 * rank;
    IBTK::IBTK_MPI::sumReduction<int>({ { a, 2 }, { &b, 1 } });
    if (!rank) output_file << "sumReduction of several arrays: " << a[0] << " " << a[1] << " " << b << "\n";

    if (!rank) output_file.close();
} // main

template <typename T>
bool
minReduction(T x, T exact, int proc)
{
    int min_proc = std::numeric_limits<int>::max();
    // perform a min reduction
    T val = IBTK::IBTK_MPI::minReduction(x, &min_proc);
    if (val == exact && proc == min_proc)
        return true;
    else
        return false;
}

template <typename T>
bool
maxReduction(T x, T exact, int proc)
{
    int max_proc = std::numeric_limits<int>::min();
    // perform a max reduction
    T val = IBTK::IBTK_MPI::maxReduction(x, &max_proc);
    if (val == exact && proc == max_proc)
        return true;
    else
        return false;
}

template <typename T>
bool
sumReduction(T x, T exact)
{
    T val = IBTK::IBTK_MPI::sumReduction(x);
    if (val == exact)
        return true;
    else
        return false;
}

template <typename T>
bool
allToOneSumReduction(T x, T exact, int root)
{
    IBTK::IBTK_MPI::allToOneSumReduction(&x, 1, root);
    if (IBTK::IBTK_MPI::getRank() == root)
    {
        if (x == exact)
            return true;
        else
            return false;
    }
    else
    {
        return true;
    }
}

template <typename T>
bool
bcast(T x, T exact, int root)
{
    T other = IBTK::IBTK_MPI::bcast(x, root);
    if (other == exact)
        return true;
    else
        return false;
}

template <typename T>
bool
sendAndRecv(T x, T exact)
{
    T other;
    int num_nodes = IBTK::IBTK_MPI::getNodes();
    int size = 1;
    IBTK::IBTK_MPI::send(&x, size, (IBTK::IBTK_MPI::getRank() + 1) % num_nodes, false);
    IBTK::IBTK_MPI::recv(&other, size, (IBTK::IBTK_MPI::getRank() - 1 + num_nodes) % num_nodes, false);

    if (other == exact)
        return true;
    else
        return false;
}

template <typename T>
bool
allGather(T x)
{
    std::vector<T> other(IBTK::IBTK_MPI::getNodes());
    IBTK::IBTK_MPI::allGather(x, other.data());

    bool passed = true;
    for (int i = 0; i < IBTK::IBTK_MPI::getNodes(); ++i)
        if (other[i] != i) passed = false;
    return passed;
}
