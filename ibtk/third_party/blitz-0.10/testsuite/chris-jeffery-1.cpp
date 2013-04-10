#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int zippy()
{ return 0; }

int main()
{
    if (zippy())
    {
#ifdef BZ_HAVE_COMPLEX
#ifdef BZ_HAVE_COMPLEX_MATH1

{
    Array<complex<float>,1> test(5);
    test = sqrt(test);

    Array<double,1> test2(5);
    test2 = sqrt(test2);
}

#endif  // BZ_HAVE_COMPLEX_MATH1

{
  Array<complex<float>,1> test(5),test2(5) ;
  Array<float,1> test3(5) ;
  test3 = 0.64f ;
  test2 = complex<float>(3,3) ;
  test = 2*test3*test2 ;
  test = test3*test2;
}

#endif  // BZ_HAVE_COMPLEX

    }

    return 0;
}

