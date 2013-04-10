#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

typedef TinyVector<uint8_t,4> tv;

int main()
{

  {
    // First test bitwise update operators taking a scalar
    const tv a(1,2,4,8);
    tv b(a);
    b>>=1;
    BZTEST(all(b==tv(0,1,2,4)));

    b&=3;
    BZTEST(all(b==tv(0,1,2,0)));

    b|=4;
    BZTEST(all(b==tv(4,5,6,4)));

    b^=3;
    BZTEST(all(b==tv(7,6,5,7)));

    b= ~b;
    BZTEST(all(b==tv(~7,~6,~5,~7)));
  }
  
  {
    // now bitwise update operators taking a tv
    const tv a(1,2,4,8);
    tv b(a);
    b>>= tv(0,1,2,3);
    BZTEST(all(b==tv(1,1,1,1)));

    b<<= tv(3,2,1,0);
    BZTEST(all(b==tv(8,4,2,1)));

    b=tv(0xff, 0xff, 0xff, 0xff);
    b&=tv(0,1,3,5);
    BZTEST(all(b==tv(0,1,3,5)));

    b|=tv(5,3,1,0);
    BZTEST(all(b==tv(5,3,3,5)));

    b^=tv(1,3,5,0);
    BZTEST(all(b==tv(4,0,6,5)));
  }

  {
    // bitwise operators taking a scalar
    const tv a(1,2,4,8);

    BZTEST(all((a>>1)==tv(0,1,2,4)));
    BZTEST(all((a<<1)==tv(2,4,8,16)));

    BZTEST(all((a&3)==tv(1,2,0,0)));

    BZTEST(all((a|4)==tv(5,6,4,12)));

    BZTEST(all((a^3)==tv(2,1,7,11)));
  }

  {
    // bitwise operators taking two tvs
    const tv a(1,2,4,8);
    const tv b(0,1,2,3);
    const tv c(3,1,7,1);

    BZTEST(all((a>>b)==tv(1,1,1,1)));

    BZTEST(all((a<<b)==tv(1,4,16,64)));

    BZTEST(all((a&c)==tv(1,0,4,0)));

    BZTEST(all((a|c)==tv(3,3,7,9)));

    BZTEST(all((a^c)==tv(2,3,3,9)));

  }
  
    return 0;
}

