#include <vector>

#include "example.cpp"

int
main(int argc, char** argv)
{
    std::vector<double> U_err;
    run_example(argc, argv, U_err);
    return 0;
}
