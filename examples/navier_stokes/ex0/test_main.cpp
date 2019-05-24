#include "example.cpp"
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

int ex_argc;
char** ex_argv;
bool ex_runs;
static const double REL_ERROR = 1.0e-5;
std::vector<double> bench_u_err, bench_p_err;
std::vector<double> u_err, p_err;
double bench;
double actual;
static const int L1_IDX = 0;
static const int L2_IDX = 1;
static const int MAX_IDX = 2;

// Set names of test based on if compiled with 2D or 3D libraries
#if (NDIM == 2)
#define TEST_CASE_NAME navier_stokes_ex0_2d
#endif
#if (NDIM == 3)
#define TEST_CASE_NAME navier_stokes_ex0_3d
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

// tests for p error

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

    // error recorded from main2d running input2d.test Oct 4, 2016
    // benchmark error in u

    bench_u_err.resize(3);
    if (NDIM == 2)
    {
        bench_u_err[L1_IDX] = 0.00357601; // 2d L1Norm
        bench_u_err[L2_IDX] = 0.00439633; // 2d L2Norm
        bench_u_err[MAX_IDX] = 0.0417876; // 2d maxNorm
    }
    else if (NDIM == 3)
    {
        bench_u_err[L1_IDX] = 0.000510578; // 3d L1Norm
        bench_u_err[L2_IDX] = 0.000477978; // 3d L2Norm
        bench_u_err[MAX_IDX] = 0.0055429;  // 3d maxNorm
    }

    // benchmark error in p
    bench_p_err.resize(3);
    if (NDIM == 2)
    {
        bench_p_err[L1_IDX] = 0.0219484; // 2d L1Norm
        bench_p_err[L2_IDX] = 0.0295763; // 2d L2Norm
        bench_p_err[MAX_IDX] = 0.220349; // 2d maxNorm
    }
    else if (NDIM == 3)
    {
        bench_p_err[L1_IDX] = 0.0021559;   // 3d L1Norm
        bench_p_err[L2_IDX] = 0.00249808;  // 3d L2Norm
        bench_p_err[MAX_IDX] = 0.00843954; // 3d maxNorm
    }

    ex_argc = argc;
    ex_argv = argv;
    ex_runs = run_example(ex_argc, ex_argv, u_err, p_err);
    return RUN_ALL_TESTS();
}
