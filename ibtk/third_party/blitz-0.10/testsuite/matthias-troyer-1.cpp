#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

#if 0

BZ_DECLARE_STENCIL2(kinEnergy,A,B)
    B=Laplacian3D_stencilop(A);
BZ_END_STENCIL_WITH_SHAPE(shape(-1,-1,-1),shape(1,1,1))

typedef complex<double> T_num;

typedef Array<T_num,3> array3d;

int main()
{
const int N=5;
array3d A(N,N,N);
array3d B(N,N,N);

// Fill a three-dimensional array with a Gaussian function
firstIndex i;
secondIndex j;
thirdIndex k;
float midpoint = 15/2.;
float c = - 1/3.0;
A = exp(c * (sqr(i-midpoint) + sqr(j-midpoint)
    + sqr(k-midpoint)));

applyStencil(kinEnergy(), A, B);

Array<T_num,1> out_view(B.data(),shape(N*N*N));
cout << out_view;
}
#endif

BZ_DECLARE_STENCIL2(footprint,A,B)
    B = Laplacian2D4_stencilop(A);
BZ_END_STENCIL_WITH_SHAPE(shape(-2,-2),shape(+2,+2))

int main()
{
    int N = 9;
    Array<double,2> A(N,N), B(N,N);
    A = 0;
    A(4,4) = 1;
    applyStencil(footprint(), A, B);
    cout << B(Range(2,6),Range(2,6)) << endl;
}

