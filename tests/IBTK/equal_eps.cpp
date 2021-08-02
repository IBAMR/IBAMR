// ---------------------------------------------------------------------
//
// Copyright (c) 2020 - 2021 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <ibtk/AppInitializer.h>
#include <ibtk/IBTKInit.h>
#include <ibtk/ibtk_utilities.h>

#include <tbox/MathUtilities.h>

#include <ibamr/app_namespaces.h>

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    IBTKInit ibtk_init(argc, argv, MPI_COMM_WORLD);

    std::ofstream out("output");

    auto test = [&out](const double a, const double b, std::function<bool(double, double)> compare, std::string type)
    {
        out.precision(std::numeric_limits<double>::max_digits10);
        out << type << " test:\n";
        out << "Testing a = " << a << " and b = " << b << "\n";
        if (compare(a, b))
        {
            out << "EQUAL!\n";
        }
        else
        {
            out << "NOT EQUAL!\n";
        }
        out << "\n";
    };
    std::vector<std::function<bool(double, double)> > fcns; // = {IBTK::rel_equal_eps, IBTK::abs_equal_eps,
                                                            // MathUtilities<double>::equalEps};
    fcns.push_back([](double a, double b) -> bool { return IBTK::rel_equal_eps(a, b); });
    fcns.push_back([](double a, double b) -> bool { return IBTK::abs_equal_eps(a, b); });
    fcns.push_back(MathUtilities<double>::equalEps);
    std::vector<std::string> strs = { "Relative", "Absolute", "SAMRAI" };

    auto loop_test = [&](const double a, const double b)
    {
        for (unsigned int i = 0; i < fcns.size(); ++i)
        {
            test(a, b, fcns[i], strs[i]);
        }
    };

    std::vector<double> add_vec = { 0.0, 1.0 };
    for (const auto& add : add_vec)
    {
        const double min = add + std::numeric_limits<double>::min();
        const double zero = add + 0.0;
        const double eps = add + std::numeric_limits<double>::epsilon();
        const double eps_1 = add + 1.0e-8;

        loop_test(min, zero);
        loop_test(zero, zero);
        loop_test(eps, zero);
        loop_test(eps, min);
        loop_test(eps_1, min);
        loop_test(eps_1, zero);
        loop_test(eps_1, eps);
    }

    return 0;
} // main
