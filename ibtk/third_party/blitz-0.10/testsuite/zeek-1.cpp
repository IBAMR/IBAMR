#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    BZ_USING_NAMESPACE(tensor)

    Array<double, 2> A(4, 4);
    A = 1,  2,  3,  4,
        5,  6,  7,  8,
        9,  10, 11, 12,
        13, 14, 15, 16;

    Array<double, 1> a(4);
    a = sum(A, j);
    Array<double, 1> b(sum(A, j));

    BZTEST(all(a == b));

    return 0;
}

