#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<int,1> d(5);
    d=1,-1,2,-2,0;
    BZTEST(d(0) == 1);
    BZTEST(d(1) == -1);
    BZTEST(d(2) == 2);
    BZTEST(d(3) == -2);
    BZTEST(d(4) == 0);

    Array<double,1> c(5);
    c=1.0,-1.0,2.0,-2.0,0.0;
    BZTEST(c(0) == 1);
    BZTEST(c(1) == -1);
    BZTEST(c(2) == 2);
    BZTEST(c(3) == -2);
    BZTEST(c(4) == 0);
    
    // and finally test initializations where implicit conversions to
    // T_numtype are necessary

    Array<double,1> b(5);
    b=1,-1,2,-2,0;
    BZTEST(b(0) == 1);
    BZTEST(b(1) == -1);
    BZTEST(b(2) == 2);
    BZTEST(b(3) == -2);
    BZTEST(b(4) == 0);

    return 0;
}

