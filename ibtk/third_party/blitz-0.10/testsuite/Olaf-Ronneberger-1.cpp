#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main( void)
{
  Array<float,3> A(1,10,10);
  Array<float,2> B(10,10);

  A = 0.0f;
  B = 0.0f;

  beginCheckAssert();  
  A(1, Range::all(), Range::all()) = B;  // first index should have been 0!
  endCheckAssert();

  beginCheckAssert();
  A(1, Range::all(), Range::all()) = B;  // not okay: shape error
  endCheckAssert();

  beginCheckAssert();
  A(Range::all(), 10, Range::all()) = 0;  // not okay: second index out of range
  endCheckAssert();

  A(0, Range::all(), Range::all()) = B;  // okay

  exit(0);
}

