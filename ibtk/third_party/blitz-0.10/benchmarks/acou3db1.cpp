#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

void setupInitialConditions(Array<float,3>& P1, Array<float,3>& P2,
    Array<float,3>& P3, Array<float,3>& c, int N);

float acoustic3D_BlitzRaw(int N, int niters)
{
    // Allocate the arrays as a group.  Blitz++ will interlace them in
    // memory, improving data locality.

    Array<float,3> P1(N,N,N), P2(N,N,N), P3(N,N,N), c(N,N,N);
    Range I(1,N-2), J(1,N-2), K(1,N-2);

    setupInitialConditions(P1, P2, P3, c, N);

    for (int iter=0; iter < niters; ++iter)
    {
        P3(I,J,K) = (2-6*c(I,J,K)) * P2(I,J,K)
          + c(I,J,K)*(P2(I-1,J,K) + P2(I+1,J,K) + P2(I,J-1,K) + P2(I,J+1,K)
          + P2(I,J,K-1) + P2(I,J,K+1)) - P1(I,J,K);

        P1 = P2;
        P2 = P3;
    }

    return P1(N/2,N/2,N/2);
}

