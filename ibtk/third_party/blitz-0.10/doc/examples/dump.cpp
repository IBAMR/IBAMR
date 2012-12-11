#include <blitz/array.h>

using namespace blitz;

int main()
{
    Array<float,4> A(3,7,8,2,FortranArray<4>());
    A.dumpStructureInformation(cerr);
    return 0;
}

