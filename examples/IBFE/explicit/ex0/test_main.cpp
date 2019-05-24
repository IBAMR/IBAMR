#include "example.cpp"
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

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
    // on commit a205f3676cd in August 2018.
    //
    // benchmark error in Q

    // 2d
    // Error in u at time 0.01953125:
    // L1-norm:  4.63166742e-05
    // L2-norm:  6.774968741e-05
    // max-norm: 0.0003396844457

    bench_u_err.resize(3);
    if (NDIM == 2)
    {
        bench_u_err[L1_IDX] = 4.6316950814468583e-05;  // 2d L1Norm
        bench_u_err[L2_IDX] = 6.7750371581914943e-05;  // 2d L2Norm
        bench_u_err[MAX_IDX] = 0.00033968718606779092; // 2d maxNorm
    }

    // 2d
    // Error in p at time 0.0185546875:
    // L1-norm:  0.2293675067
    // L2-norm:  0.9715761264
    // max-norm: 7.749929666

    bench_p_err.resize(3);
    if (NDIM == 2)
    {
        bench_p_err[L1_IDX] = 0.22936750682483956; // 2d L1Norm
        bench_p_err[L2_IDX] = 0.97157612634508372; // 2d L2Norm
        bench_p_err[MAX_IDX] = 7.7499296947099507; // 2d maxNorm
    }

    ex_argc = argc;
    ex_argv = argv;
    ex_runs = run_example(ex_argc, ex_argv, u_err, p_err);
    return RUN_ALL_TESTS();
}
