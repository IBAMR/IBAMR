#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
}

void foo()
{
    const int N = 10;
    Array<float,1> A(N), B(N), C(N), D(N), E(N);
    A = B*C + D*E;
}

