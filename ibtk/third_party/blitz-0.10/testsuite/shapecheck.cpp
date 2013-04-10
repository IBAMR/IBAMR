#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<int,2> A(4,4), B(3,4), C(4,4), D(3,4);
    Array<int,3> E(4,4,4);

    A = 0;
    B = 0;
    C = 0;
    D = 0;
    E = 0;

    beginCheckAssert();
    A = D;
    endCheckAssert();
 
    A = C;
 
    beginCheckAssert();
    A = B + D;    
    endCheckAssert();

    return 0;
}


