#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

//declaration of the stencil kinEnergy
BZ_DECLARE_STENCIL2(kinEnergy,A,B)
B=Laplacian3D_stencilop(A);
BZ_END_STENCIL_WITH_SHAPE(shape(-1,-1,-1),shape(+1,+1,+1))

typedef Array<complex<double>,3> array3d;

int main()
{
const int N=5;
array3d A(N,N,N);
array3d B(N,N,N);
// Fill a three-dimensional array with a Gaussian function
firstIndex i;
secondIndex j;
thirdIndex k;
float midpoint = 3.;
float c = - 1.;
//A = exp(c * (sqr(i-midpoint) + sqr(j-midpoint)
//    + sqr(k-midpoint)));
A = zip( exp(c * (sqr(i-midpoint) + sqr(j-midpoint)
    + sqr(k-midpoint))), 0.0, complex<double>());

applyStencil(kinEnergy(),A,B);

Array<complex<double>,1> a_view(A.data(),shape(N*N*N));
cout << a_view;
Array<complex<double>,1> out_view(B.data(),shape(N*N*N));
cout << out_view<<endl;
}

