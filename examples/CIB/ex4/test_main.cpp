#include "example.cpp"
#include <gtest/gtest.h>

int ex_argc;
char** ex_argv;
bool ex_runs;
double analytic_end_u;
static const double SHORT_END_TIME = 0.015;
static const double LONG_END_TIME = 1.798348053;
static const double REL_ERROR = 1.0e-1;

double acceptable_error;
double ex_end_time, ex_end_u;

#if (NDIM == 2)
#define TEST_CASE_NAME CIB_ex4_2d
#endif
#if (NDIM == 3)
#define TEST_CASE_NAME CIB_ex4_3d
#endif

TEST(TEST_CASE_NAME, example_runs)
{
    EXPECT_TRUE(ex_runs);
}

TEST(TEST_CASE_NAME, u_convergence)
{
    EXPECT_LE(std::abs(analytic_end_u - ex_end_u), acceptable_error);
}

int
main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    ex_argc = argc;
    ex_argv = argv;
    ex_runs = run_example(ex_argc, ex_argv, ex_end_time, ex_end_u);

    if (std::abs(ex_end_time - SHORT_END_TIME) <= REL_ERROR)
    {
        analytic_end_u = 0.0171;
        acceptable_error = analytic_end_u;
    }
    else if (std::abs(ex_end_time - LONG_END_TIME) <= REL_ERROR)
    {
        analytic_end_u = 0.9552;
        acceptable_error = analytic_end_u * 0.25;
    }

    return RUN_ALL_TESTS();
}
