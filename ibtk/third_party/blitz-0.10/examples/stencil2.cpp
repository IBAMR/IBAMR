#include <blitz/array.h>
#include <blitz/array/stencil-et.h>

BZ_USING_NAMESPACE(blitz)

// This example illustrates a simple use of stencils.
// It applies a Laplacian operator to an array which is
// all zero, except for a single 1 in the centre.
// This results in the stencil pattern (although reversed
// and upside-down).

int main()
{
    const int N = 9;

    // Create some arrays
    Array<float,2> A(N,N), B(N,N), C(N,N), D(N,N);

    // Initialize A to 1 in the centre, 0 elsewhere
    A = 0;
    A(N/2,N/2) = 1.0;

    // Get the interior of the arrays.  Since central42 is the
    // stencil [ 1 -4 6 -4 1 ], we have to leave a padding of 2
    // cells around the interior.  The range 2..N-3 does this.
    Range I(2,N-3);
    Array<float,2> Ai = A(I,I), Bi = B(I,I), Ci = C(I,I), Di = D(I,I);

    // These declarations let us call the first dimension "x"
    // and the second dimension "y".
    const int x = firstDim, y = secondDim;

    // Apply the stencils.  Note that Laplacian2D is identical to
    // central22(Ai,x) + central22(Ai,y).

    Bi = Laplacian2D(Ai);
    Ci = central42(Ai,x) + central42(Ai,y);
    Di = central22(Ai,x) + central22(Ai,y);

    // Display some of the results

    Range middle(N/4,3*N/4);

    cout << "Laplacian2D: " << B(middle,middle) << endl
         << "central42: " << C(middle,middle) << endl
         << "central22: " << D(middle,middle) << endl;

    return 0;
}

