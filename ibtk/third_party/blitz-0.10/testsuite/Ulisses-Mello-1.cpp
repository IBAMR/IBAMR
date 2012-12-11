#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<double,3> geo(5,10,15);
    Array<double,3> bio;
        
    geo = -1.0;
        
    bio.resize(geo.shape());
    BZTEST(geo.extent(0) == bio.extent(0)
        && geo.extent(1) == bio.extent(1)
        && geo.extent(2) == bio.extent(2));

    bio = geo;
    return 0;
}

