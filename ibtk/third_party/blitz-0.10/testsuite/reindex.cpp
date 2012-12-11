#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    BZ_USING_NAMESPACE(blitz::tensor)

    Array<int,2> A(3,3);
    A = i*3 + j;

    // 0, 1, 2
    // 3, 4, 5
    // 6, 7, 8

    A.reindexSelf(TinyVector<int,2>(3,4));
    BZTEST(A(3,4) == 0);
    BZTEST(A(5,6) == 8);

    Array<int,2> B = A.reindex(TinyVector<int,2>(0,0));
    BZTEST(B(0,0) == 0);
    BZTEST(B(2,2) == 8);

    return 0;
}

