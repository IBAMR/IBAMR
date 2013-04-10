#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    // 3x3 C-style row major storage, base zero
    Array<int,2> A(3, 3);

    // 3x3 Fortran-style column major storage, base one
    Array<int,2> B(3, 3, ColumnMajorArray<2>());

    // A custom storage format: 
    // Indices have range 0..3, 0..3
    // Column major ordering
    // Rows are stored ascending, columns stored descending
    GeneralArrayStorage<2> storage;
    storage.ordering() = firstRank, secondRank;
    storage.base() = 0, 0;
    storage.ascendingFlag() = true, false;

    //    Array<int,2> C(3, 3, storage);
    Array<int,2> C;
    C.setStorage(storage);
    C.resize(3,3);

    // Set each array equal to
    // [ 1 2 3 ]
    // [ 4 5 6 ]
    // [ 7 8 9 ]

    A = 1, 2, 3,
        4, 5, 6, 
        7, 8, 9;

    // Comma-delimited lists initialize in memory-storage order only.
    // Hence we list the values in column-major order to initialize B:

    B = 1, 4, 7, 2, 5, 8, 3, 6, 9;

    BZTEST(all(A == B));
    BZTEST(B(2,1) == 8);
    BZTEST(B(0,2) == 3);

    // Array C is stored in column major, plus the columns are stored
    // in descending order

    C = 3, 6, 9, 2, 5, 8, 1, 4, 7;
    BZTEST(all(A == C));
    BZTEST(count(A==C) == 9);

    Array<int,2> D(3,3);
    D = A + B + C;

    Array<int,2> E(3,3);
    E = 3 * A;
    BZTEST(all(D == E));
    BZTEST(!all(D == A));

    return 0;
}

