#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main(int argc, char* argv[]) {
    Array<float,2> image(100, 200);
    for (int x=0;x<100;x++)
        for (int y=0;y<200;y++)
            image(x,y) = 1000*y+x;

    RectDomain<2>  rect_domain(shape(0,0),shape(99,0));
    Array<float,2> slice = image(rect_domain);

    unsigned k = 0;
    for (Array<float,2>::iterator i=slice.begin();i!=slice.end();++i)
        BZTEST(*i==k++);

    slice.reverseSelf(0);
    for (Array<float,2>::iterator i=slice.begin();i!=slice.end();++i)
        BZTEST(*i==--k);

    return 0;
} 
