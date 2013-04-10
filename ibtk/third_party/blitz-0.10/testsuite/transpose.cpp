#include "testsuite.h"
#include <blitz/array.h>
#include <random/discrete-uniform.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
  // try with all kinds of different sizes to make sure it works for
  // aligned and nonaligned lengths.
  for(int i=0;i<20;++i) {
    const int sx=ranlib::DiscreteUniform<int>(10).random();
    const int sy=ranlib::DiscreteUniform<int>(10).random();
    const int sz=ranlib::DiscreteUniform<int>(10).random();
    Array<int,3> A(sx,sy,sz, contiguousArray);

//    A.dumpStructureInformation();

    A.transposeSelf(secondDim, thirdDim, firstDim);

    BZTEST(A.ordering(0) == 1 && A.ordering(1) == 0 && A.ordering(2) == 2);
    BZTEST(A.length(0) == sy && A.length(1) == sz && A.length(2) == sx);
    BZTEST(A.stride(0) == sz && A.stride(1) == 1 && A.stride(2) == sz*sy);
  }

//    A.dumpStructureInformation();

    Array<int,2> B(3,3);
    B = 0, 1, 2,
        3, 4, 5,
        6, 7, 8;

//    cout << B << endl;

    B.transposeSelf(secondDim, firstDim);
    Array<int,2> C(3,3);
    C = 0, 3, 6,
        1, 4, 7,
        2, 5, 8;
    BZTEST(count(B==C) == 9);

//    cout << B << endl;

    Array<int,2> D(3,3,FortranArray<2>());
    D = 0, 1, 2,
        3, 4, 5,
        6, 7, 8;
    D.transposeSelf(secondDim, firstDim);
    Array<int,2> E(3,3,FortranArray<2>());
    E = 0, 3, 6,
        1, 4, 7,
        2, 5, 8;
    BZTEST(count(D==E) == 9);

    return 0;
}

