/*****************************************************************************
 * array.cpp        Blitz++ Array stencilling example
 *****************************************************************************/

#include <blitz/array.h>   

BZ_USING_NAMESPACE(blitz)

int main()
{
    int N = 64;

    // Create three-dimensional arrays of float
    Array<float,3> A(N,N,N), B(N,N,N);

    // Set up initial conditions: +30 C over an interior block,
    // and +22 C elsewhere
    A = 22.0;

    Range interior(N/4,3*N/4);
    A(interior,interior,interior) = 30.0;

    int numIters = 301;

    Range I(1,N-2), J(1,N-2), K(1,N-2);

#ifdef BZ_HAVE_STD
#ifdef BZ_ARRAY_SPACE_FILLING_TRAVERSAL
    generateFastTraversalOrder(TinyVector<int,2>(N-2,N-2));
#endif
#endif

    for (int i=0; i < numIters; ++i)
    {
        double c = 1/6.5;
 
        B(I,J,K) = c * (.5 * A(I,J,K) + A(I+1,J,K) + A(I-1,J,K)
            + A(I,J+1,K) + A(I,J-1,K) + A(I,J,K+1) + A(I,J,K-1));

        A(I,J,K) = c * (.5 * B(I,J,K) + B(I+1,J,K) + B(I-1,J,K)
            + B(I,J+1,K) + B(I,J-1,K) + B(I,J,K+1) + B(I,J,K-1));

        // Output the result along a line through the centre
        for (int j=0; j < 8; ++j)
            cout << setprecision(2) << A(N/2,N/2,j*N/8) << " ";

        cout << endl;
        cout.flush();
    }

    return 0;
}

