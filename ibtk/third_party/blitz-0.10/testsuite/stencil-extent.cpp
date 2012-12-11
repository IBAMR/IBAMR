#include <testsuite.h>
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

// test that stencil extent detection works

BZ_DECLARE_STENCIL2(test, P1, P2)
P2 = forward11_stencilop(P1,0)+forward21_stencilop(P1,1);
BZ_END_STENCIL

typedef Array<double,2> array2d;

int main()
{
const int N=5;
array2d A(N,N);
array2d B(N,N);
// Fill a three-dimensional array with a Gaussian function
firstIndex i;
secondIndex j;
thirdIndex k;
float midpoint = 3.;
float c = - 1.;

 A = exp(c * (sqr(i-midpoint) + sqr(j-midpoint)));
 B = 42;
 //cout << A << endl;
 //cout << B << endl;

applyStencil(test(),A,B);

 BZTEST(B(0,0)!=42);
 BZTEST(B(4,1)==42);
 BZTEST(B(2,2)!=42);
 BZTEST(B(1,3)==42);
 BZTEST(B(1,4)==42);
 //cout << B << endl;

 
}

