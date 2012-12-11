// Nested heterogeneous arrays

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<Array<int,1>,1> A(3);

    A(0).resize(3);
    A(0) = 0, 1, 2;

    A(1).resize(5);
    A(1) = 5, 7, 18, 2, 1;

    BZ_USING_NAMESPACE(blitz::tensor);

    A(2).resize(4);
    A(2) = pow2(i+1);

    cout << "A = " << A << endl;

    // A = [ [ 0 1 2 ] [ 5 7 18 2 1 ] [ 1 4 9 16 ] ]

    return 0;
}

