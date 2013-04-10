#define BZ_TAU_PROFILING

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    TAU_PROFILE("main()", "int ()", TAU_DEFAULT);

    const int N = 32;

    Array<float,2> A(N,N), B(N,N), C(N,N), D(N,N), E(N,N);
    A = 5.0;
    B = 0.0;
    C = 0.0;

    for (int i=0; i < 20; ++i)
    {
        D = A + B + C;
        D /= sum(pow2(D));
        A = B * cos(D) + C * sin(D);
        B += exp(-D);

        float x = sum(A);
        float y = sum(A+B);
        float z = sum(sqr(A)+sqr(B));
        C = x*A+y*B+z*C;

        D = exp(-sqr(A)-sqr(B));
        E = A + B + C + D;
        float q = min(A);
        float r = max(B);
    }

    return 0;
}

