#include "testsuite.h"

#include<blitz/array.h>
#include <blitz/tinymat2.h>
#include <blitz/tinymat2.cc>

BZ_USING_NAMESPACE(blitz)

int main()
{
  TinyMatrix<int,2,3> A,B,C;
    A = 7;
    BZTEST(A(1,2) == 7);
    BZTEST(A(0,0) == 7);

    B = 0, 1, 2,
        2, 3, 4;
    BZTEST(B(0,0) == 0);
    BZTEST(B(1,2) == 4);

    C = A*B;

    BZTEST(C(0,0) == 0);
    BZTEST(C(1,2) == 28);

    return 0;
}
