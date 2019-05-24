#include "example.cpp"
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

int ex_argc;
char** ex_argv;
bool ex_runs;
static const double REL_ERROR = 1.0e-8;
std::vector<double> bench_Q_err;
std::vector<double> Q_err;
double bench;
double actual;
static const int L1_IDX = 0;
static const int L2_IDX = 1;
static const int MAX_IDX = 2;

// Set names of test based on if compiled with 2D or 3D libraries
#if (NDIM == 2)
#define TEST_CASE_NAME adv_diff_ex0_2d
#endif
#if (NDIM == 3)
#define TEST_CASE_NAME adv_diff_ex0_3d
#endif

TEST(TEST_CASE_NAME, example_runs)
{
    EXPECT_TRUE(ex_runs);
}

TEST(TEST_CASE_NAME, Q_L1Norm)
{
    bench = bench_Q_err[L1_IDX];
    actual = Q_err[L1_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

TEST(TEST_CASE_NAME, Q_L2Norm)
{
    bench = bench_Q_err[L2_IDX];
    actual = Q_err[L2_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

TEST(TEST_CASE_NAME, Q_MaxNorm)
{
    bench = bench_Q_err[MAX_IDX];
    actual = Q_err[MAX_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

int
main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);

    // error recorded from main2d running input2d.test Oct 28, 2016
    // benchmark error in Q

    // 2d
    // Error in Q at time 0.2:
    // L1-norm:  0.0002722711746
    // L2-norm:  0.0004827779198
    // max-norm: 0.002898502357

    bench_Q_err.resize(3);
    if (NDIM == 2)
    {
        bench_Q_err[L1_IDX] = 0.0002722711746; // 2d L1Norm
        bench_Q_err[L2_IDX] = 0.0004827779198; // 2d L2Norm
        bench_Q_err[MAX_IDX] = 0.002898502357; // 2d maxNorm
    }

    // 3d
    // Error in Q at time 0.2:
    // L1-norm:  0.1668961774
    // L2-norm:  0.3181721473
    // max-norm: 2.450737448

    else if (NDIM == 3)
    {
        bench_Q_err[L1_IDX] = 0.1668961774; // 3d L1Norm
        bench_Q_err[L2_IDX] = 0.3181721473; // 3d L2Norm
        bench_Q_err[MAX_IDX] = 2.450737448; // 3d maxNorm
    }

    ex_argc = argc;
    ex_argv = argv;
    ex_runs = run_example(ex_argc, ex_argv, Q_err);
    return RUN_ALL_TESTS();
}
