// Check blitz/minmax.h, which provides a templated, type-promoting
// version of min(a,b) and max(a,b).

// Some MS Windows and boost header files provide macro definitions of min/max
// that can conflict with our templated functions in the blitz::extrema 
// namespace.  Placing these macros here as a test that we can work around
// this problem by enclosing our min/max invocations in a set of parentheses.
#define min(a,b) (a < b ? a : b)
#define max(a,b) (a > b ? a : b)

#include "testsuite.h"
#include <blitz/minmax.h>

int main()
{
    double a = 1.0;
    float b = 2.0;
    int c = 3;

    BZTEST(static_cast<int>((blitz::extrema::min)(a,b)) == 1);
    BZTEST(static_cast<int>((blitz::extrema::max)(a,b)) == 2);
    BZTEST(static_cast<int>((blitz::extrema::min)(a,c)) == 1);
    BZTEST(static_cast<int>((blitz::extrema::min)(b,c)) == 2);
    BZTEST(static_cast<int>((blitz::extrema::max)(b,c)) == 3);
    BZTEST(static_cast<int>((blitz::extrema::max)(a,c)) == 3);

    return 0;
}

