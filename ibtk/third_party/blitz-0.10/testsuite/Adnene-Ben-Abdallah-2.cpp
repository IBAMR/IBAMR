#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main() 
{
  BZ_USING_NAMESPACE(tensor)

  Array<double,1> X(4), Y1(4), Y2(4);
  X = 1, 2, 3, 4;
  Array<double,2> A(4,4);
  A = 1,  2,  3,  4,
      5,  6,  7,  8,
      9,  10, 11, 12,
      13, 14, 15, 16;

  Y1 = 5, 6, 7, 8;
  Y2 = 5, 6, 7, 8;
// this does not work under IBM xlC compiler!
#if !defined(_IBMR2)
  Y1 += sum(A(i,j) * X(j), j);
  Y2 = Y2 + sum(A(i,j)* X(j), j);
#endif
  BZTEST(count(Y1==Y2)==4);

  return 0;
}

