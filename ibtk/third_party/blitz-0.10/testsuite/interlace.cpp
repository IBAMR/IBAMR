// Intended coverage: blitz/arrayinter.h (Array interlacing)

#include "testsuite.h"

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<int,2> A, B;
    interlaceArrays(shape(32,33), A, B);
    BZTEST(A.stride(1) == 2);
    BZTEST(B.stride(1) == 2);
    BZTEST(&A(0,0) + 1 == &B(0,0));
    BZTEST(&A(5,0) + 1 == &B(5,0));
    A = 0;
    BZTEST(A.rows() == 32);
    BZTEST(A.cols() == 33);
    BZTEST(count(A == 0) == 32*33);
    B = 1;
    BZTEST(B.rows() == 32);
    BZTEST(B.cols() == 33);
    BZTEST(count(B == 1) == 32*33);
    BZTEST(count(A == 0) == 32*33);


    Array<int,2> C;
    interlaceArrays(shape(13,17), A, B, C);
    BZTEST(&A(0,0) + 1 == &B(0,0));
    BZTEST(&A(0,0) + 2 == &C(0,0));
    A = 0;
    B = 1;
    C = 2;
    B = 1;
    A = 0;
    BZTEST(count(A == 0) == 13*17);
    BZTEST(count(B == 1) == 13*17);
    BZTEST(count(C == 2) == 13*17);

    Array<int,2> D, E, F, G, H;
    interlaceArrays(shape(11,19), A, B, C, D, E, F, G, H);
    BZTEST(A.stride(1) == 8);
    BZTEST(E.stride(1) == 8);
    BZTEST(H.stride(1) == 8);
    BZTEST(&A(0,0)+1 == &B(0,0));
    BZTEST(&A(0,0)+2 == &C(0,0));
    BZTEST(&A(0,0)+3 == &D(0,0));
    BZTEST(&A(0,0)+4 == &E(0,0));
    BZTEST(&A(0,0)+5 == &F(0,0));
    BZTEST(&A(0,0)+6 == &G(0,0));
    BZTEST(&A(0,0)+7 == &H(0,0));
    A = 0;
    B = 1;
    C = 2;
    D = 3;
    E = 4;
    F = 5;
    G = 6;
    H = 7;
    G = 6;
    F = 5;
    E = 4;
    D = 3;
    C = 2;
    B = 1;
    A = 0;
    BZTEST(count(A == 0) == 11*19);
    BZTEST(count(B == 1) == 11*19);
    BZTEST(count(C == 2) == 11*19);
    BZTEST(count(D == 3) == 11*19);
    BZTEST(count(E == 4) == 11*19);
    BZTEST(count(F == 5) == 11*19);
    BZTEST(count(G == 6) == 11*19);
    BZTEST(count(H == 7) == 11*19);
    BZTEST(count(A != 1) == 11*19);
}

