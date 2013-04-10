#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<float,1> A(4), B(4);
    firstIndex i;
    A = i;
    B = 0, 1, 2, 3;
    BZTEST(count(A==B)==4);
    return 0;
}

