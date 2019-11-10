#include "example.cpp"
#include <gtest/gtest.h>

int ex_argc;
char** ex_argv;
std::vector<double> Q_err;
bool ex_runs;
bool run_example(int, char**, std::vector<double>&);

// Set names of test based on if compiled with 2D or 3D libraries
#if (NDIM == 2)
#define TEST_CASE_NAME level_set_ex1_2d
#endif
#if (NDIM == 3)
#define TEST_CASE_NAME level_set_ex1_3d
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
    ex_runs = run_example(ex_argc, ex_argv, Q_err);
    return RUN_ALL_TESTS();
}
