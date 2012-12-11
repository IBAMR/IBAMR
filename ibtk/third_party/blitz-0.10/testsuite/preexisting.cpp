#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

// tests that preexisting memory constructors work with padding
int main()
{
  Array<int,2> b(3,3);
  b = 1,2,3,4,5,6,7,8,9;

  {
    int p[] = {1,2,3,4,5,6,7,8,9};
    Array<int,2> a(p, shape(3,3), neverDeleteData);
    BZTEST(all(a==b));
  }

  {
    int* pp;
    pp = new int[3*3];
    Array<int,2> a(pp, shape(3,3), deleteDataWhenDone);
    a = 1,2,3,4,5,6,7,8,9;
    BZTEST(all(a==b));
  }


  const int w =simdTypes<int>::vecWidth;
  if(w>1) {
    int p[] = {1,2,3,4,5,6,7,8,9};
    beginCheckAssert();
    Array<int,2> a(p, shape(3,3), neverDeleteData, FortranArray<2>(paddedData));
    endCheckAssert();
  }

    return 0;
}

