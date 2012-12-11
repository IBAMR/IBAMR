#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

const int w=simdTypes<int>::vecWidth;

void check1DArrays();
void check3DArrays();
void checkEmpty1DArray(Array<int,1>& A);
void checkFull1DArray(Array<int,1>& A, int length);
void checkFull3DArray(Array<int,3>& A, int L, int N, int M);
void checkFull3DFortranArray(Array<int,3>& A, int L, int N, int M);

int main()
{
    check1DArrays();
    check3DArrays();
}

void check1DArrays()
{
    Array<int,1> A;

    checkEmpty1DArray(A);
    Array<int,1> B = A;
    checkEmpty1DArray(B);
    Array<int,1> C(A.copy());
    checkEmpty1DArray(C);

    const int N = 10;
    A.resize(N);
    checkFull1DArray(A, N);
    B.resizeAndPreserve(2*N);
    checkFull1DArray(B,2*N);
    Array<int,1> D(N);
    checkFull1DArray(D, N);
}

void checkEmpty1DArray(Array<int,1>& A)
{
    BZTEST(A.base(0) == 0);
    beginCheckAssert();  A.columns();  endCheckAssert();
    BZTEST(A.data() == 0);
    BZTEST(A.dataZero() == 0);
    BZTEST(A.dataFirst() == 0);
    beginCheckAssert();  A.depth();  endCheckAssert();
    BZTEST(A.dimensions() == 1);
    BZTEST(A.extent(0) == 0);
    BZTEST(A.isMajorRank(0));
    BZTEST(!A.isMinorRank(0));
    BZTEST(A.isRankStoredAscending(0));
    BZTEST(A.lbound(0) == 0);
    BZTEST(A.length(0) == 0);
    BZTEST(A.numElements() == 0);
    BZTEST(A.ordering(0) == 0);
    BZTEST(A.rank() == 1);
    BZTEST(A.rows() == 0);
    BZTEST(A.size() == 0);
    BZTEST(A.stride(0) == 0);
    BZTEST(A.ubound(0) == -1);
    BZTEST(A.zeroOffset() == 0);
}

void checkFull1DArray(Array<int,1>& A, int length)
{
    BZTEST(A.base(0) == 0);
    beginCheckAssert();  A.columns();  endCheckAssert();
    BZTEST(A.data() != 0);
    BZTEST(A.dataZero() == A.data());
    BZTEST(A.dataFirst() == A.data());
    beginCheckAssert();  A.depth();  endCheckAssert();
    BZTEST(A.dimensions() == 1);
    BZTEST(A.extent(0) == length);
    BZTEST(A.isMajorRank(0));
    BZTEST(!A.isMinorRank(0));
    BZTEST(A.isRankStoredAscending(0));
    BZTEST(A.lbound(0) == 0);
    BZTEST(A.length(0) == length);
    BZTEST(A.numElements() == length);
    BZTEST(A.ordering(0) == 0);
    BZTEST(A.rank() == 1);
    BZTEST(A.rows() == length);
    BZTEST(A.size() == length);
    BZTEST(A.stride(0) == 1);
    BZTEST(A.ubound(0) == length-1);
    BZTEST(A.zeroOffset() == 0);
}

void check3DArrays()
{
    Array<int,3> AA(2*w,3*w,4*w);
    checkFull3DArray(AA,2*w,3*w,4*w);
    Array<int,3> A(2,3,4);
    checkFull3DArray(A,2,3,4);

    Array<int,3> BB(2*w,3*w,4*w,FortranArray<3>());
    checkFull3DFortranArray(BB,2*w,3*w,4*w);
    Array<int,3> B(2,3,4,FortranArray<3>());
    checkFull3DFortranArray(B,2,3,4);
}

void checkFull3DArray(Array<int,3>& A, int L, int M, int N)
{
    BZTEST(A.base(0) == 0);
    BZTEST(A.columns() == M);
    BZTEST(A.data() != 0);
    BZTEST(A.dataZero() == A.data());
    BZTEST(A.dataFirst() == A.data());
    BZTEST(A.depth() == N);
    BZTEST(A.dimensions() == 3);
    BZTEST(A.extent(0) == L);
    BZTEST(A.extent(1) == M);
    BZTEST(A.extent(2) == N);
    beginCheckAssert(); A.extent(3); endCheckAssert();    
    BZTEST(A.isMajorRank(0));
    BZTEST(!A.isMajorRank(1));
    BZTEST(!A.isMajorRank(2));
    BZTEST(!A.isMinorRank(0));
    BZTEST(A.isMinorRank(1));
    BZTEST(A.isMinorRank(2));
    BZTEST(A.isRankStoredAscending(0));
    BZTEST(A.isRankStoredAscending(1));
    BZTEST(A.isRankStoredAscending(2));
    BZTEST(A.lbound(0) == 0);
    BZTEST(A.lbound(1) == 0);
    BZTEST(A.lbound(2) == 0);
    BZTEST(A.numElements() == L*N*M);
    BZTEST(A.ordering(0) == 2);
    BZTEST(A.ordering(1) == 1);
    BZTEST(A.ordering(2) == 0);
    BZTEST(A.rank() == 3);
    BZTEST(A.rows() == L);
    BZTEST(A.size() == L*N*M);
#if BZ_PAD_ARRAYS
    if(N%w) {
      beginCheckAssert();
      BZASSERT(A.stride(0) == N*M);
      BZASSERT(A.stride(1) == N);
      endCheckAssert();
    }
    else 
#endif
      {
	BZTEST(A.stride(0) == N*M);
	BZTEST(A.stride(1) == N);
      }
    BZTEST(A.stride(2) == 1);
    BZTEST(A.ubound(0) == L-1);
    BZTEST(A.ubound(1) == M-1);
    BZTEST(A.ubound(2) == N-1);
    BZTEST(A.zeroOffset() == 0);
}

void checkFull3DFortranArray(Array<int,3>& A, int L, int M, int N)
{
    BZTEST(A.base(0) == 1);
    BZTEST(A.columns() == M);
    BZTEST(A.data() != 0);
    BZTEST(A.dataZero() == A.data() - A.stride(0) - A.stride(1) - A.stride(2));
    BZTEST(A.dataFirst() == A.data());
    BZTEST(A.depth() == N);
    BZTEST(A.dimensions() == 3);
    BZTEST(A.extent(0) == L);
    BZTEST(A.extent(1) == M);
    BZTEST(A.extent(2) == N);
    beginCheckAssert(); A.extent(3); endCheckAssert();
    BZTEST(!A.isMajorRank(0));
    BZTEST(!A.isMajorRank(1));
    BZTEST(A.isMajorRank(2));
    BZTEST(A.isMinorRank(0));
    BZTEST(A.isMinorRank(1));
    BZTEST(!A.isMinorRank(2));
    BZTEST(A.isRankStoredAscending(0));
    BZTEST(A.isRankStoredAscending(1));
    BZTEST(A.isRankStoredAscending(2));
    BZTEST(A.lbound(0) == 1);
    BZTEST(A.lbound(1) == 1);
    BZTEST(A.lbound(2) == 1);
    BZTEST(A.numElements() == L*N*M);
    BZTEST(A.ordering(0) == 0);
    BZTEST(A.ordering(1) == 1);
    BZTEST(A.ordering(2) == 2);
    BZTEST(A.rank() == 3);
    BZTEST(A.rows() == L);
    BZTEST(A.size() == L*N*M);
    BZTEST(A.stride(0) == 1);
#if BZ_PAD_ARRAYS
    if(L%w) {
      beginCheckAssert();
      BZASSERT(A.stride(1) == L);
      BZASSERT(A.stride(2) == L*M);
      endCheckAssert();
    }
    else
#endif
      {
	BZTEST(A.stride(1) == L);
	BZTEST(A.stride(2) == L*M);
      }
    BZTEST(A.ubound(0) == L);
    BZTEST(A.ubound(1) == M);
    BZTEST(A.ubound(2) == N);
    BZTEST(- A.zeroOffset() == A.stride(0) + A.stride(1) + A.stride(2));
}

