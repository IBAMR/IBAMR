#include <vector>

#include "example.cpp"

int
main(int argc, char** argv)
{
    std::vector<double> C_err;
    run_example(argc, argv, C_err);
    return 0;
}
