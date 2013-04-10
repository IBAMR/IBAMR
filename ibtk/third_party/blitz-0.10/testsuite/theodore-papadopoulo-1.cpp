#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<double,1> F(4);

    F = -1.0, -2.0, -3.0, -4.0;

    TinyVector<int,1> i = maxIndex(F);
    double f = max(F);

    BZTEST(i[0] == 0);
    BZTEST(f == -1.0);
    return 0;
}

