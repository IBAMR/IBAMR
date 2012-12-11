// Test real, imag, conj, abs, arg, polar
// Added some tests of complex<long double>

#include "testsuite.h"

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

#ifndef BZ_HAVE_COMPLEX_FCNS


int main()
{
    cout << "(Cannot test complex math functions-- not supported on "
            "this platform)" << endl;

    return 0;
}

#else
typedef Array<complex<float>, 1> CArray;
static const double eps = 0.0001;
int main()
{
  const int sz=7;
    CArray A(sz);
    A(0) = complex<float>(1,2);
    A(1) = complex<float>(3,4);

    Array<float,1> Ar = real(A);
    BZTEST(int(Ar(0)) == 1 && int(Ar(1)) == 3);

    Array<float,1> Ai = imag(A);
    BZTEST(int(Ai(0)) == 2 && int(Ai(1)) == 4);

    CArray Ac(sz);
    Ac = conj(A);

    BZTEST(Ac(0) == complex<float>(1,-2));
    BZTEST(Ac(1) == complex<float>(3,-4));

    Array<float,1> Ab(sz);
    Ab = abs(A);
    BZTEST(fabs(Ab(0) - 2.236068) < eps);
    BZTEST(fabs(Ab(1) - 5.0) < eps);

    Ab = arg(A);
    BZTEST(fabs(Ab(0) - atan(2.0)) < eps);
    BZTEST(fabs(Ab(1) - atan(4.0/3.0)) < eps);

    Array<float,1> r(sz), theta(sz);
    r(0) = 4.0f;
    r(1) = 15.0f;
    theta(0) = float(3.141592/3.0);
    theta(1) = float(3.0*3.141592/2.0);
    Ac = blitz::polar(r,theta);
    BZTEST(fabs(real(Ac(0)) - 2) < eps);
    BZTEST(fabs(imag(Ac(0)) - 3.4641012) < eps);
    BZTEST(fabs(real(Ac(1)) - 0.0) < eps);
    BZTEST(fabs(imag(Ac(1)) + 15.0) < eps);

    Array<complex<long double>,1> A11(5),B11(5),C11(5);
    A11=1,2,3,4,5;
    B11=1,2,3,4,5;
    C11=A11+B11;
    BZTEST(fabs(real(C11(0)) - 2.) < eps);
    C11=A11/B11;
    BZTEST(fabs(real(C11(1)) - 1.) < eps);
    C11=1.0l/A11;
    BZTEST(fabs(real(C11(2)) - 1/3.) < eps);
    C11=A11/1.0l;
    BZTEST(fabs(real(C11(3)) - 4.) < eps);
    C11=complex<long double>(0,1)/A11;
    BZTEST(fabs(imag(C11(4)) - 1/5.) < eps);
    C11=A11/complex<long double>(0,1);
    BZTEST(fabs(imag(C11(0)) - -1.) < eps);

    return 0;
}
#endif

