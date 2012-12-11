#include <math.h>

int main()
{
}

int N;
float* A, * B, * C, * D, * E;

void foo()
{
    for (int i=0; i < N; ++i)
        A[i] = B[i]*C[i] + D[i]*E[i];

    for (int i=0; i < N; ++i)
        A[i] = B[i]+C[i] + D[i]*cos(E[i]);

    for (int i=0; i < N; ++i)
        A[i] = B[i]*sin(C[i]) + D[i]*sin(E[i]);

    for (int i=0; i < N; ++i)
        B[i] = C[i] + D[i];

    for (int i=0; i < N; ++i)
        A[i] = A[i] + B[i] + C[i] + D[i];

    for (int i=0; i < N; ++i)
        A[i] = B[i]*C[i]*C[i] + D[i]*E[i];

    for (int i=0; i < N; ++i)
        A[i] = B[i]-C[i] + log(D[i])*cos(E[i]);

    for (int i=0; i < N; ++i)
        A[i] = B[i]*sin(C[i]) + tan(D[i])/E[i];

    for (int i=0; i < N; ++i)
        B[i] = C[i] - D[i];

    for (int i=0; i < N; ++i)
        A[i] = A[i]+B[i]*C[i]+D[i];

    for (int i=0; i < N; ++i)
        C[i] = A[i]+B[i]+C[i]+D[i]+E[i];

    for (int i=0; i < N; ++i)
        E[i] = (A[i]+B[i])*(C[i]-D[i]);

    for (int i=0; i < N; ++i)
        D[i] = A[i]/B[i]*C[i]/D[i];

    for (int i=0; i < N; ++i)
        B[i] = (A[i]*B[i]) + (A[i]/B[i]);

    for (int i=0; i < N; ++i)
        D[i] = sin(A[i]) + sin(E[i]);
}

