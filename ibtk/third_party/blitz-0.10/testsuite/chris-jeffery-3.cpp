#include "testsuite.h"
#include <blitz/array.h>


BZ_USING_NAMESPACE(blitz)

int main()
{
Array<float,2> test(8,8), test2(5,5) ;

test = 5;

Range I(2,6) ;
Range J(3,7) ;

// Koenig lookup hack
#if defined(__GNUC__) && (__GNUC__ < 3)
test2 = where(blitz::operator> (test(I,J), test(I-1,J)), 0, test(I,J));
#else
test2 = where(test(I,J) > test(I-1,J), 0, test(I,J));
#endif

BZTEST(test2(3,3) == 5);

cout << test2 << endl ;

}


