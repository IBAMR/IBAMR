#include <blitz/array.h>
#include <blitz/timer.h>

BZ_USING_NAMESPACE(blitz)

BZ_DECLARE_STENCIL4(acoustic2D_stencil,P1,P2,P3,c)
  P3 = 2 * P2 + c * Laplacian2D(P2) - P1;
BZ_END_STENCIL

int benchmark(int N, int nIterations, int blockSize)
{
    Array<float,2> P1, P2, P3, c;
    allocateArrays(shape(N,N), P1, P2, P3, c);

    // Initial conditions: obviously in a real application these
    // wouldn't be zeroed...
    Range I(0,blockSize-1), J(0,blockSize-1);

    P1(I,J) = 0;
    P2(I,J) = 0;
    P3(I,J) = 0;
    c(I,J) = 0;

    for (int i=0; i < nIterations; ++i)
    {
        // Apply the stencil object to the arrays
        applyStencil(acoustic2D_stencil(), P1(I,J), P2(I,J), P3(I,J), c(I,J));

        // Set [P1,P2,P3] <- [P2,P3,P1] to set up for the next
        // time step
        cycleArrays(P1,P2,P3);
    }

    return 0;
}

int main()
{
    Timer timer;

    cout << "N\tMflops" << endl;

    const int blockSize = 27;

    for (int N=2000; N < 2100; ++N)
    {
        double stencilPoints = pow(blockSize-2,2.0);
        int nIterations = 5000;

        timer.start();
        benchmark(N, nIterations, blockSize);
        timer.stop();

        double flops = (4 + 7) * stencilPoints * nIterations;
        double Mflops = flops / timer.elapsedSeconds() / 1.0E+6;
        cout << N << "\t" << Mflops << endl;
    }

    return 0;
}

