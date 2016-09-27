#include <gtest/gtest.h>
#include "example.cpp"
#include <cmath>

int ex_argc;
char** ex_argv;
bool result;
double run_example(int, char**);
double EPSILON = 0.01;
double std_2d_u_error =  0.0147675; // recorded from main2d running same number of timesteps September 26, 2016
double std_3d_u_error =  0.0055429; // recorded from main3d running same number of timesteps September 26, 2016
double u_error;

TEST(navier_stokes_ex0, 2d_uMaxNorm) {
    EXPECT_LE(std::abs(u_error - std_2d_u_error), EPSILON);
}

TEST(navier_stokes_ex0, 3d_uMaxNorm) {
  EXPECT_LE(std::abs(u_error - std_3d_u_error), EPSILON);
}

int main( int argc, char** argv ) {
    testing::InitGoogleTest( &argc, argv );
    ex_argc = argc;
    ex_argv = argv;
    u_error = run_example(ex_argc, ex_argv);
//    test_errors[1] = run_example(ex_argc, ex_argv)[1];
    return RUN_ALL_TESTS( );
}
