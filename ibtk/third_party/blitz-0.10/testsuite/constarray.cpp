#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)
BZ_USING_NAMESPACE(blitz::tensor)

// Test that subarrays may be taken of const arrays.

void foo(const Array<int,2>& A)
{
    int n1 = sum(A(1,Range::all()));
    BZTEST(n1 == 46);

    Array<int,1> B(4), C(4);
    B = sum(A(j,i),j);
    C = 60, 64, 68, 72;
    BZTEST(count(B==C)==4);

    n1 = sum(A(Range::all(), Range::all()));
    BZTEST(n1 == 264);

    n1 = sum(2*A);
    BZTEST(n1 == 264*2);
}

// 0 1 2 3
// 10 11 12 13
// 20 21 22 23
// 30 31 32 33

int main()
{
    Array<int,2> A(4,4);

    A = 10*i+j;

    foo(A);

    return 0;
}

