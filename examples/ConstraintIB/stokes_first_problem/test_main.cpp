#include <gtest/gtest.h>

#include "example.cpp"

int ex_argc;
char** ex_argv;
bool ex_runs;
bool run_example(int, char**);

#if (NDIM == 2)
#define TEST_CASE_NAME ConstraintIB_stokes_first_problem
#endif
#if (NDIM == 3)
#define TEST_CASE_NAME ConstraintIB_stokes_first_problem
#endif

TEST(TEST_CASE_NAME, example_runs)
{
    EXPECT_TRUE(ex_runs);
}

int
main(int argc, char** argv)
{
    testing::InitGoogleTest(&argc, argv);
    ex_argc = argc;
    ex_argv = argv;
    ex_runs = run_example(ex_argc, ex_argv);
    return RUN_ALL_TESTS();
}
