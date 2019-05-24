#include "example.cpp"
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

int ex_argc;
char** ex_argv;
bool ex_runs;
static const double REL_ERROR = 1.0e-8;
std::vector<double> bench_U_err;
std::vector<double> U_err;
double bench;
double actual;
static const int L1_IDX = 0;
static const int L2_IDX = 1;
static const int MAX_IDX = 2;

// Set names of test based on if compiled with 2D or 3D libraries
#if (NDIM == 2)
#define TEST_CASE_NAME adv_diff_ex1_2d
#endif
#if (NDIM == 3)
#define TEST_CASE_NAME adv_diff_ex1_3d
#endif

TEST(TEST_CASE_NAME, example_runs)
{
    EXPECT_TRUE(ex_runs);
}

TEST(TEST_CASE_NAME, U_L1Norm)
{
    bench = bench_U_err[L1_IDX];
    actual = U_err[L1_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

TEST(TEST_CASE_NAME, U_L2Norm)
{
    bench = bench_U_err[L2_IDX];
    actual = U_err[L2_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

TEST(TEST_CASE_NAME, U_MaxNorm)
{
    bench = bench_U_err[MAX_IDX];
    actual = U_err[MAX_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

int
main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);

    // error recorded from main2d running input2d.test Oct 28, 2016
    // benchmark error in Q

    // 2d
    // Error in U at time 0.125:
    // L1-norm:  0.001721208912
    // L2-norm:  0.001361852495
    // max-norm: 0.001586238609

    bench_U_err.resize(3);
    if (NDIM == 2)
    {
        bench_U_err[L1_IDX] = 0.001721208912;  // 2d L1Norm
        bench_U_err[L2_IDX] = 0.001361852495;  // 2d L2Norm
        bench_U_err[MAX_IDX] = 0.001586238609; // 2d maxNorm
    }

    // 3d
    // Error in U at time 0.5:
    // L1-norm:  0.001088142652
    // L2-norm:  0.0007832940189
    // max-norm: 0.001055833193

    else if (NDIM == 3)
    {
        bench_U_err[L1_IDX] = 0.001088142652;  // 3d L1Norm
        bench_U_err[L2_IDX] = 0.0007832940189; // 3d L2Norm
        bench_U_err[MAX_IDX] = 0.001055833193; // 3d maxNorm
    }

    ex_argc = argc;
    ex_argv = argv;
    ex_runs = run_example(ex_argc, ex_argv, U_err);
    return RUN_ALL_TESTS();
}
