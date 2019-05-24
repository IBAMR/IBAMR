#include "example.cpp"
#include <vector>

int
main(int argc, char** argv)
{
    std::vector<double> u_err, p_err;
    run_example(argc, argv, u_err, p_err);
    return 0;
}
