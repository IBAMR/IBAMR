#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    // NEEDS_WORK: resize

    // One-dimensional resizeAndPreserve
    Array<int,1> A(5);
    A = 0, 1, 2, 3, 4;
    A.resizeAndPreserve(10);
    A(Range(5,toEnd)) = 0;
  
    BZTEST(A.lbound(0) == 0);
    BZTEST(A.ubound(0) == 9);
 
    Array<int,1> B(5);
    B = 0, 1, 2, 3, 4;
     
    BZTEST(all(A(Range(0,4)) == B));
    BZTEST(all(A(Range(5,9)) == 0));

    A.resizeAndPreserve(5);
    BZTEST(all(A == B));
    A.resizeAndPreserve(2);
    BZTEST((A(0) == 0) && (A(1) == 1));
    BZTEST(A.lbound(0) == 0);
    BZTEST(A.ubound(0) == 1);

    // Two-dimensional resizeAndPreserve, with Fortran-style arrays
    Array<int,2> C(3,3,FortranArray<2>());
    C = 0, 1, 2, 
        3, 4, 5,
        6, 7, 8;
    C.resizeAndPreserve(2,4);
    // valgrind here correctly complains about use of uninitialized
    // values, because the (i,4) elements are uninitialized. harmless.
    BZTEST(all(C(Range(1,2),Range(1,4)) == C(Range(1,2),Range(1,4))));
    BZTEST(C.lbound(0) == 1);
    BZTEST(C.lbound(1) == 1);
    BZTEST(C.ubound(0) == 2);
    BZTEST(C.ubound(1) == 4);

    // Three dimensional resize
    Array<int,3> G(4,5,6);
    G.resize(7,8,9);
    BZTEST(G.extent(0) == 7);
    BZTEST(G.extent(1) == 8);
    BZTEST(G.extent(2) == 9);
    return 0;
}

