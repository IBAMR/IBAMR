#include "testsuite.h"
#include <blitz/array.h>


BZ_USING_NAMESPACE(blitz)

int main() {
  Array<int, 2> A(4);
  Array<int, 1> B(4), C(4);

  firstIndex i;
  secondIndex j;

  A = i+j;

  B = max(A,j);

  C = 3, 4, 5, 6;
  BZTEST(count(B == C) == 4);

  return 0;
}


