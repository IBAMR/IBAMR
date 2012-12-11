#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

// Test index operators with different bases. Tests bug XXX.

int main()
{
  Array<int,1> bar(3); // c-Array
  Array<int,1> qux(3,fortranArray); // fortran
  Array<int,1> foo(Range(5,7)); // custom
  Array<int,1> base(Range(5,15));
  Array<int,1> subbase(base(Range(7,11,2)));

  blitz::firstIndex ii;
  // Indexed assignment from 0
  bar = ii;
  BZTEST(bar(0)==0);
  BZTEST(bar(1)==1);
  BZTEST(bar(2)==2);
  cout << bar;
  // FortranArray
  qux = ii;
  BZTEST(qux(1)==1);
  BZTEST(qux(2)==2);
  BZTEST(qux(3)==3);
  cout << qux;
  // Indexed assignment on non-base-0
  foo = ii;
  BZTEST(foo(5)==5);
  BZTEST(foo(6)==6);
  BZTEST(foo(7)==7);
  cout << foo;
  subbase=ii;
  BZTEST(subbase(5)==5);
  BZTEST(subbase(6)==6);
  BZTEST(subbase(7)==7);
  cout << subbase;
  return 0;
}

