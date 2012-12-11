/*
 * Complex arrays example
 */

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

const double pi = 3.14159265358979323846264338327950288;

#ifndef BZ_HAVE_COMPLEX
int main()
{
    cout << "Complex number support is required to compile this example."
         << endl;
    return 0;
}

#else

int main()
{
    const int N = 16;

    Array<complex<double>,1> A(N);
    Array<double,1> theta(N);

    BZ_USING_NAMESPACE(blitz::tensor);

    // Fill the theta array with angles from 0..2 Pi, evenly spaced
    theta = (2 * pi * i) / N;

    // Set A[i] = cos(theta[i]) + _I * sin(theta[i])
    A = zip(cos(theta), sin(theta), complex<double>());

    cout << A << endl;

#ifdef BZ_HAVE_COMPLEX_MATH1
    // Here's another way of doing it, which eliminates the need for
    // the theta array:
    // Set A[i] = exp(Pi i _I / N)
    A = exp(zip(0.0, (2 * pi * i) / N, complex<double>()));

    cout << A << endl;
#endif

    return 0;
}
#endif // BZ_HAVE_COMPLEX
