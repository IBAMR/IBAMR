#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
  Array<float,1> f(3);
  Array<double,1> d(3),e(3);
  f=1,2,4;
  e=1,2,4;

  d=cast<double>(f);
  BZTEST(all(f==e));

    return 0;
}

