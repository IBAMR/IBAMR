/*
 * This benchmark times the performance of B += sqr(A), where A and B
 * are complex<double> arrays.
 *
 * Note: need to use -mv8 for SPARC v8.
 */

#include <blitz/array.h>
#include <blitz/timer.h>

BZ_USING_NAMESPACE(blitz)

#include <blitz/vector.h>

typedef Array<complex<double>,1> CArray;

void setup(Array<complex<double>,1>& A, Array<complex<double>,1>& B)
{
    int n = A.extent(firstDim);

    for (int i=0; i < n; ++i)
    {
        double x=-10.+20./(n-1.0)*i;
        A(i)=sin(x);
        B(i)=sin(x);
    }
}

void version1(CArray& A, CArray& B, int nIters)
{
    Timer timer;

    // Array notation
    setup(A, B);
    timer.start();
    for (int i=0; i < nIters; ++i)
    {
        B += A*A;
    }
    timer.stop();
    cout << "Time using array notation b += a*a: " << timer.elapsedSeconds()
         << endl;
}

void version2(CArray& A, CArray& B, int nIters)
{
    // Array notation, using sqr(a)
    Timer timer;

    setup(A,B);
    timer.start();
    for (int i=0; i < nIters; ++i)
    {
        B += sqr(A);
    }
    timer.stop();
    cout << "Time using array notation b += sqr(a): " << timer.elapsedSeconds()
         << endl;
}

void version2c(CArray& A, CArray& B)
{
    B += sqr(A);
}

void version2b(CArray& A, CArray& B, int nIters)
{
    // Array notation, using sqr(a)
    Timer timer;

    setup(A,B);
    timer.start();
    for (int i=0; i < nIters; ++i)
    {
        version2c(A,B);
    }
    timer.stop();
    cout << "Time using array notation b += sqr(a): " << timer.elapsedSeconds()
         << endl;
}

void version3(CArray& A, CArray& B, int nIters)
{
    Timer timer;

    int N = A.extent(firstDim);

    // Low-level implementation
    setup(A,B);
    timer.start();
    for (int i=0; i < nIters; ++i)
    {
        for (int j=0; j < N; ++j)
            B(j) += A(j) * A(j);
    }
    timer.stop();
    cout << "Time using low-level version: " << timer.elapsedSeconds()
         << endl;
}

void version4(CArray& A, CArray& B, int nIters)
{
    Timer timer;

    struct cmplx {
        double re, im;
    };
    cmplx* a = (cmplx*)A.data();
    cmplx* b = (cmplx*)B.data();
    setup(A,B);
    int N = A.extent(firstDim);

    timer.start();
    for (int i=0; i < nIters; ++i)
    {
        for (int j=0; j < N; ++j)
        {
            double ar = a[j].re;
            double ai = a[j].im;
            b[j].re += ar*ar - ai*ai;
            b[j].im += 2 * ar * ai;
        }
    }
    timer.stop();
    cout << "Time using really low-level version: " << timer.elapsedSeconds()
         << endl;
}
           
int run(int N, int nIters)
{
    Array<complex<double>,1> A(N), B(N);

    version1(A,B,nIters);
    version2(A,B,nIters);
    version2b(A,B,nIters);
    version3(A,B,nIters);
    version4(A,B,nIters);

    return 0;
}

int main()
{
    cout << "In-cache:" << endl;
    run(256,39063);

    cout << endl << "Out-of-cache:" << endl;
    run(1000000,10);
}

