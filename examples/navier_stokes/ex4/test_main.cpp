#include "example.cpp"
#include <cmath>
#include <gtest/gtest.h>

int ex_argc;
char** ex_argv;
bool ex_runs;
bool run_example(int, char**);

// Set names of test based on if compiled with 2D or 3D libraries
#if (NDIM == 2)
#define TEST_CASE_NAME navier_stokes_ex4_2d
#endif
#if (NDIM == 3)
#define TEST_CASE_NAME navier_stokes_ex4_3d
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
