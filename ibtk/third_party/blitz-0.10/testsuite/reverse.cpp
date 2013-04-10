#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<int,2> A(3,4), B(3,4);
    A = 0, 1, 2, 3,
        4, 5, 6, 7,
        8, 9, 10, 11;
   
    B = 3, 2, 1, 0,
        7, 6, 5, 4,
        11, 10, 9, 8;
    A.reverseSelf(secondDim);

    BZTEST(count(A==B) == 12);

    A.reverseSelf(secondDim);
    B = 0, 1, 2, 3,
        4, 5, 6, 7,
        8, 9, 10, 11;
    BZTEST(count(A==B) == 12);

    A.reverseSelf(firstDim);
    B = 8, 9, 10, 11,
        4, 5, 6, 7,
        0, 1, 2, 3;
    BZTEST(count(A==B) == 12);

    A.reverseSelf(secondDim);
    B = 11, 10, 9, 8,
        7, 6, 5, 4,
        3, 2, 1, 0;
    BZTEST(count(A==B) == 12);

    Array<int,2> E(3,4);
    E = 0, 1, 2, 3,
        4, 5, 6, 7,
        8, 9, 10, 11;

    Array<int,2> C = E(Range(1,2),Range(2,3));
    Array<int,2> D(2,2);
    D = 6, 7,
        10, 11;
    BZTEST(count(C==D) == 4);

    C.reverseSelf(firstDim);
    D = 10, 11,
        6, 7;
    BZTEST(count(C==D) == 4);

    C.reverseSelf(firstDim);
    D = 6, 7, 
       10, 11;
    BZTEST(count(C==D) == 4);

    C.reverseSelf(secondDim);
    D = 7, 6,
       11, 10;

    BZTEST(count(C==D) == 4);
    C.reverseSelf(firstDim);
    D = 11, 10,
         7, 6;

    BZTEST(count(C==D) == 4);
    return 0;
}

