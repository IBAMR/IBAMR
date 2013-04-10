#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

template <typename T>
void subsample(Array<T, 2>& in, Array<T, 2>& sub1)
{
    Range all = Range::all();

    sub1.resize((in.extent(0)-3)/2+1, in.extent(1));

    Range R0(1, in.extent(0)-2, 2);

    sub1(all, all) = (2*in(R0, all) + in(R0+1, all) + in(R0-1, all))/4.0;
}

int main()
{
    Array<float,2> A(9,5), B;

    BZ_USING_NAMESPACE(blitz::tensor)

    A = 20*i+j;

    subsample(A,B);

    Array<float,2> C(4,5);
    C =  20, 21, 22, 23, 24,
         60, 61, 62, 63, 64,
        100, 101, 102, 103, 104,
        140, 141, 142, 143, 144;
    BZTEST(count(B==C) == 20);

    return 0;
}

