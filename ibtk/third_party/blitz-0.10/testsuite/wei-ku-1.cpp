#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
  BZ_USING_NAMESPACE(blitz::tensor)
  Range all = Range::all();

  {
    Array<int,5> A(6,5,4,3,2);
    int ld_begin = 1;
    int ld_end = 3;
    Array<int,5> B(6,5,2,3,2);

    A = 11111 + 10000*i + 1000*j + 100*k + 10*l + m;

    B = A(all, all, Range(ld_begin, ld_end-1), all, all);
    // cout << B << endl;    

    BZTEST(B(0,0,0,0,0) == 11211);
    BZTEST(B(1,1,1,1,1) == 22322);
    BZTEST(B(5,4,1,2,1) == 65332);
  }

  {
    Array<int,5> B(10,1,25,6,6);
    Array<int,5> A(10,1,125,6,6);
    int ld_begin = 1, ld_end = 26;
    A = 11111 + 10000*i + 1000*j + 100*k + 10*l + m;
    B = A(all, all, Range(ld_begin, ld_end-1), all, all);

    BZTEST(B(0,0,0,0,0) == 11211);
    BZTEST(B(1,0,1,1,1) == 21322);
    BZTEST(B(5,0,3,4,5) == 61556);
  }

  return 0;
}

