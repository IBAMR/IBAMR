#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

// test that math functions operating on integers return doubles and not ints

int main()
{
  Array<double,1> a(10);
  Array<int,1> b(10);
  b=tensor::i;
  a=sqrt(b);

  cout << a ; 
  BZTEST(a(2)==sqrt(2));
  BZTEST(a(5)==sqrt(5));
  
  a=sqrt(tensor::i); 
  cout << a ; 
  BZTEST(a(2)==sqrt(2));
  BZTEST(a(5)==sqrt(5));

  // we can't really test that the powN and abs are not cast to
  // doubles because it's not really visible from the outside.

  return 0;
}

