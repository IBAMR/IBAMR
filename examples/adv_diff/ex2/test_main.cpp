#include "example.cpp"
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

int ex_argc;
char** ex_argv;
bool ex_runs;
static const double REL_ERROR = 1.0e-8;
std::vector<double> bench_C_err;
std::vector<double> C_err;
double bench;
double actual;
static const int L1_IDX = 0;
static const int L2_IDX = 1;
static const int MAX_IDX = 2;

// Set names of test based on if compiled with 2D or 3D libraries
#if (NDIM == 2)
#define TEST_CASE_NAME adv_diff_ex2_2d
#endif
#if (NDIM == 3)
#define TEST_CASE_NAME adv_diff_ex2_3d
#endif

TEST(TEST_CASE_NAME, example_runs)
{
    EXPECT_TRUE(ex_runs);
}

TEST(TEST_CASE_NAME, C_L1Norm)
{
    bench = bench_C_err[L1_IDX];
    actual = C_err[L1_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

TEST(TEST_CASE_NAME, C_L2Norm)
{
    bench = bench_C_err[L2_IDX];
    actual = C_err[L2_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

TEST(TEST_CASE_NAME, C_MaxNorm)
{
    bench = bench_C_err[MAX_IDX];
    actual = C_err[MAX_IDX];
    EXPECT_LE(std::abs((actual - bench)), (bench * REL_ERROR));
}

int
main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);

    // error recorded from main2d running input2d.test Oct 31, 2016
    // benchmark error in U

    // 2d
    // Error in U at time 2:
    // L1-norm:  8.678686881e-05
    // L2-norm:  0.0001199845941
    // max-norm: 0.0007675396914

    bench_C_err.resize(3);
    if (NDIM == 2)
    {
        bench_C_err[L1_IDX] = 8.678686881e-05;  // 2d L1Norm
        bench_C_err[L2_IDX] = 0.0001199845941;  // 2d L2Norm
        bench_C_err[MAX_IDX] = 0.0007675396914; // 2d maxNorm
    }

    // 3d
    // No benchmark data at this time

    //    else if (NDIM == 3) {
    //        bench_C_err[L1_IDX] = 0; //3d L1Norm
    //        bench_C_err[L2_IDX] = 0; //3d L2Norm
    //        bench_C_err[MAX_IDX] = 0;  //3d maxNorm
    //    }

    ex_argc = argc;
    ex_argv = argv;
    ex_runs = run_example(ex_argc, ex_argv, C_err);
    return RUN_ALL_TESTS();
}
