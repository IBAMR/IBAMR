#include "ex5.cpp"

int example_argc;
char** example_argv;
bool runExample(int, char**);

int main(int argc, char* argv[])
{
    example_argc = argc;
    example_argv = argv;
    runExample(example_argc, example_argv);
    return 0;
}
