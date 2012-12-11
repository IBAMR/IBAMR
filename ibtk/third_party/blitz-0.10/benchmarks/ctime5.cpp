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
    A = B+C + D*cos(E);
    A = B*sin(C) + D*sin(E);
    B = C + D;
    A = A+B+C+D;
    A = B*C*C + D*E;
    A = B-C + log(D)*cos(E);
    A = B*sin(C) + tan(D)/E;
    B = C - D;
    A = A+B*C+D;
    C = A+B+C+D+E;
    E = (A+B)*(C-D);
    D = A/B*C/D;
    B = (A*B) + (A/B);
    D = sin(A) + sin(E);
}

