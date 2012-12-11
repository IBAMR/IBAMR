// Creating arrays from array-expressions
#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<float,2> A(4,3,fortranArray);
    Array<float,2> B(4,3,fortranArray);

    A = 2.0;
    B = 1.0;

    Array<float,2> C(A+B);

    BZTEST(C.ordering(0) == A.ordering(0));
    BZTEST(C.ordering(1) == A.ordering(1));

    Array<float,2> D(4,3);

    beginCheckAssert();
    Array<float,2> C3(A+D);   // conflicting storage formats
    endCheckAssert();

    return 0;
}

