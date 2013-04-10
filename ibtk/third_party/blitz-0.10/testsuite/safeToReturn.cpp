#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

const int N=1000;

// To return an expression that references local arrays, the
// expression has to be wrapped in safeToReturn(), which turns
// FastArrayIterators into FastArrayCopyIterators which keep a copy of
// the arrays instead of referencing them.

// this is safe
BzBinaryExprResult<Multiply, Array<double,1>, Array<double,1> >::T_result::T_range_result
good()
{
  Array<double,1> a(N);
  a=tensor::i;
  Array<double,1> b(N);
  b=sin(a);
  return safeToReturn(a*b);
}

// this is unsafe
BzBinaryExprResult<Multiply, Array<double,1>, Array<double,1> >::T_result
bad()
{
  Array<double,1> a(N);
  a=tensor::i;
  Array<double,1> b(N);
  b=sin(a);
  return a*b;
}

int main()
{
  // switch to expr=bad() and it will crash.
  BzBinaryExprResult<Multiply, Array<double,1>, Array<double,1> >::T_result::T_range_result expr=good();

  //BzBinaryExprResult<Multiply, Array<double,1>, Array<double,1> >::T_result expr=bad();
  Array<double,1> result(N);
  result=expr;
  BZTEST(all(abs(result-tensor::i*sin(1.0*tensor::i))<1e-5));

    return 0;
}

