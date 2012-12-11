#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)
BZ_USING_NAMESPACE(blitz::tensor)

int main()
{
    Array<int,2> A(2,3), B(2,3);

    A = 0, 3, 5, 
        1, 6, 9;
    B = 2, 5, 1,
        9, 3, 4;

    Array<int,2> C(A+2*B);
    Array<int,2> C2(2,3);
    C2 =  4, 13,  7,
         19, 12, 17;
    BZTEST(count(C2 == C) == 6); 

    beginCheckAssert();
    Array<int,2> D(i*10+j);
    endCheckAssert();
}

