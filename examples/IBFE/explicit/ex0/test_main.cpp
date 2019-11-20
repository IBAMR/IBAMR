// ---------------------------------------------------------------------
//
// Copyright (c) 2017 - 2019 by the IBAMR developers
// All rights reserved.
//
// This file is part of IBAMR.
//
// IBAMR is free software and is distributed under the 3-clause BSD
// license. The full text of the license can be found in the file
// COPYRIGHT at the top level directory of IBAMR.
//
// ---------------------------------------------------------------------

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "example.cpp"

int ex_argc;
char** ex_argv;
bool ex_runs;
bool run_example(int, char**);
std::vector<double> bench_u_err, bench_p_err;
std::vector<double> u_err, p_err;
static const double REL_ERROR = 1.0e-8;
double bench;
double actual;
static const int L1_IDX = 0;
static const int L2_IDX = 1;
static const int MAX_IDX = 2;

// Set names of test based on if compiled with 2D or 3D libraries
#if (NDIM == 2)
#define TEST_CASE_NAME IBFE_explicit_ex0_2d
#endif
#if (NDIM == 3)
#define TEST_CASE_NAME IBFE_explicit_ex0_3d
#endif

TEST(TEST_CASE_NAME, example_runs)
{
    EXPECT_TRUE(ex_runs);
}

TEST(TEST_CASE_NAME, u_L1Norm)
{
    bench = bench_u_err[L1_IDX];
    actual = u_err[L1_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

TEST(TEST_CASE_NAME, u_L2Norm)
{
    bench = bench_u_err[L2_IDX];
    actual = u_err[L2_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

TEST(TEST_CASE_NAME, u_MaxNorm)
{
    bench = bench_u_err[MAX_IDX];
    actual = u_err[MAX_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

TEST(TEST_CASE_NAME, p_L1Norm)
{
    bench = bench_p_err[L1_IDX];
    actual = p_err[L1_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

TEST(TEST_CASE_NAME, p_L2Norm)
{
    bench = bench_p_err[L2_IDX];
    actual = p_err[L2_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

TEST(TEST_CASE_NAME, p_MaxNorm)
{
    bench = bench_p_err[MAX_IDX];
    actual = p_err[MAX_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

int
main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);

    // Error recorded from main2d running
    //     ./test2d input2d.test -ksp_rtol 1e-16 -ksp_max_it 1000
    // on commit 883709e1a9b2db99b98bcd89bcdf2eb7a97d84ef in November 2019.
    //
    // benchmark error in Q

    // 2d
    // Error in u at time 0.01953125:

    bench_u_err.resize(3);
    if (NDIM == 2)
    {
        bench_u_err[L1_IDX] = 4.5994873134723505926e-05;  // 2d L1Norm
        bench_u_err[L2_IDX] = 6.7446238822015566988e-05;  // 2d L2Norm
        bench_u_err[MAX_IDX] = 0.00033814404916269108275; // 2d maxNorm
    }

    // 2d
    // Error in p at time 0.0185546875:

    bench_p_err.resize(3);
    if (NDIM == 2)
    {
        bench_p_err[L1_IDX] = 0.22936754783016169; // 2d L1Norm
        bench_p_err[L2_IDX] = 0.97160102092371714; // 2d L2Norm
        bench_p_err[MAX_IDX] = 7.7487656327819296; // 2d maxNorm
    }

    ex_argc = argc;
    ex_argv = argv;
    ex_runs = run_example(ex_argc, ex_argv, u_err, p_err);
    return RUN_ALL_TESTS();
}
