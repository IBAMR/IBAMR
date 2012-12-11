#include "testsuite.h"
#include <blitz/tinyvec2.h>
#include <blitz/tinyvec2.cc>
#include <blitz/array.h>
#include <blitz/array.cc>

BZ_USING_NAMESPACE(blitz)

int main()
{
  TinyVector<int,5> d;
    d=1,-1,2,-2,0;
    BZTEST(d(0) == 1);
    BZTEST(d(1) == -1);
    BZTEST(d(2) == 2);
    BZTEST(d(3) == -2);
    BZTEST(d(4) == 0);

  TinyVector<double,5> c;
    c=1.0,-1.0,2.0,-2.0,0.0;
    BZTEST(c(0) == 1);
    BZTEST(c(1) == -1);
    BZTEST(c(2) == 2);
    BZTEST(c(3) == -2);
    BZTEST(c(4) == 0);
    
    // and finally test initializations where implicit conversions to
    // T_numtype are necessary

  TinyVector<double,5> b;
    b=1,-1,2,-2,0;
    BZTEST(b(0) == 1);
    BZTEST(b(1) == -1);
    BZTEST(b(2) == 2);
    BZTEST(b(3) == -2);
    BZTEST(b(4) == 0);

    return 0;
}

