#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<double,2> A;
    A.resizeAndPreserve(2,2);
    BZTEST(A.numElements() == 4);

    Array<double,5> B;
    B.resizeAndPreserve(2,2,2,2,2);
    BZTEST(B.numElements() == 32);

    return 0;
}

