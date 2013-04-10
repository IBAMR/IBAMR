#include <blitz/array.h>
#include <blitz/array/stencil-et.h>

BZ_USING_NAMESPACE(blitz)

void setupInitialConditions(Array<float,3>& P1, Array<float,3>& P2,
    Array<float,3>& P3, Array<float,3>& c, int N);

BZ_DECLARE_STENCIL4(acoustic3D, P1, P2, P3, c)
    P3 = 2 * P2 + c * Laplacian3D_stencilop(P2) - P1;
BZ_END_STENCIL

float acoustic3D_BlitzStencil(int N, int niters)
{
    Array<float,3> P1, P2, P3, c;
    allocateArrays(shape(N,N,N), P1, P2, P3, c);

    setupInitialConditions(P1, P2, P3, c, N);

    for (int iter=0; iter < niters; ++iter)
    {
        applyStencil(acoustic3D(), P1, P2, P3, c);
        cycleArrays(P1, P2, P3);
    }

    return P1(N/2,N/2,N/2);
}

