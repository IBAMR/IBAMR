#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

// test that weak references work ok for various types of arrays
// including nonunit strides and strange orderings resulting from
// slicing.

#define test_expr(d1,d2) BZTEST(all((d1)==(d2)));

int main()
{
  Array<double, 2> a(10,10);
  Array<double, 3> aa(10,10,10);
  Array<double,2> b;
  Array<double,1> c;

  a=tensor::i+10*tensor::j;
  aa=tensor::i+10*tensor::j+100*tensor::k;

  b.weakReference(a);
  test_expr(a,b);
  c.weakReference(a(1,Range::all()));
  test_expr(a(1,Range::all()), c);
  c.weakReference(a(Range::all(),2));
  test_expr(a(Range::all(),2), c);
  a.weakReference(aa(Range(2,6), Range::all(), 7));
  test_expr(aa(Range(2,6), Range::all(), 7), a);
  
  return 0;
}

