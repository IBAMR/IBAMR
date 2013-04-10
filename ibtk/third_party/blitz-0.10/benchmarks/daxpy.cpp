#include <blitz/vector2.h>
#include <blitz/array.h>
#include <random/uniform.h>
#include <blitz/benchext.h>

#ifdef BZ_HAVE_VALARRAY
 #define BENCHMARK_VALARRAY
#endif

#ifdef BENCHMARK_VALARRAY
#include <valarray>
#endif

BZ_NAMESPACE(blitz)
extern void sink();
BZ_NAMESPACE_END

BZ_USING_NAMESPACE(blitz)

#if defined(BZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES)
 #define fdaxpy   fdaxpy_
 #define daxpy    daxpy_
 #define f90daxpy f90daxpy_
 #define fidaxpy  fidaxpy_
 #define fidaxpyo fidaxpyo_
#elif defined(BZ_FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES)
 #define fdaxpy   fdaxpy__
 #define daxpy    daxpy__
 #define f90daxpy f90daxpy__
 #define fidaxpy  fidaxpy__
 #define fidaxpyo fidaxpyo__
#elif defined(BZ_FORTRAN_SYMBOLS_CAPS)
 #define fdaxpy   FDAXPY
 #define daxpy    DAXPY
 #define f90daxpy F90DAXPY
 #define fidaxpy  FIDAXPY
 #define fidaxpyo FIDAXPYO
#endif

extern "C" {
  void fdaxpy(const int& N, const double& da, double* x,
    const int& xstride, const double* y, const int& ystride);

  void daxpy(const int& N, const double& da, double* x,
    const int& xstride, const double* y, const int& ystride);

  void f90daxpy(const double& a, double* x, 
    const double* y, const int& length, const int& iters);

  void fidaxpy(const double& a, double* x, const double* y,
    const int& length, const int& iters);

  void fidaxpyo(const double& a, double* x, const double* y,
    const int& length, const int& iters);
}

void daxpyVectorVersion(BenchmarkExt<int>& bench, double a, double b);
void daxpyArrayVersion(BenchmarkExt<int>& bench, double a);
void daxpyF77Version(BenchmarkExt<int>& bench, double a);
void daxpyBLASVersion(BenchmarkExt<int>& bench, double a);
#ifdef FORTRAN_90
void daxpyF90Version(BenchmarkExt<int>& bench, double a);
#endif
#ifdef BENCHMARK_VALARRAY
void daxpyValarrayVersion(BenchmarkExt<int>& bench, double a);
#endif

int main()
{
    int numBenchmarks = 6;
#ifndef BENCHMARK_VALARRAY
		numBenchmarks--;   // No  valarray
#endif
#ifndef FORTRAN_90
		numBenchmarks--;   // No fortran 90
#endif

    BenchmarkExt<int> bench("DAXPY Benchmark", numBenchmarks);

    const int numSizes = 19;
    bench.setNumParameters(numSizes);

    Vector<int> parameters(numSizes);
    Vector<long> iters(numSizes);
    Vector<double> flops(numSizes);

    for (int i=0; i < numSizes; ++i)
    {
        parameters(i) = static_cast<int>(pow(10.0, 0.25*(i+1)));
        iters(i) = 50000000L / parameters(i);
        if (iters(i) < 2)
            iters(i) = 2;
        flops(i) = 2 * parameters(i) * 2;
    }

    bench.setParameterVector(parameters);
    bench.setIterations(iters);
    bench.setOpsPerIteration(flops);

    bench.beginBenchmarking();

    float a = .398498293819823;

    daxpyVectorVersion(bench, a, -a);
    daxpyArrayVersion(bench, a);
    daxpyF77Version(bench, a);
    daxpyBLASVersion(bench, a);
#ifdef FORTRAN_90
    daxpyF90Version(bench, a);
#endif
#ifdef BENCHMARK_VALARRAY
    daxpyValarrayVersion(bench, a);
#endif

    bench.endBenchmarking();

    bench.saveMatlabGraph("daxpy.m");

    return 0;
}


template<class T>
void initializeRandomDouble(T* data, int numElements, int stride = 1)
{
    ranlib::Uniform<T> rnd;

    for (int i=0; i < numElements; ++i)
        data[size_t(i*stride)] = rnd.random();
}

template<class T>
void initializeRandomDouble(valarray<T>& data, int numElements, int stride = 1)
{
    ranlib::Uniform<T> rnd;

    for (int i=0; i < numElements; ++i)
        data[size_t(i*stride)] = rnd.random();
}

void daxpyVectorVersion(BenchmarkExt<int>& bench, double a, double b)
{
    bench.beginImplementation("Vector<T>");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Vector<T>: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Vector<double> x(N), y(N);
        initializeRandomDouble(x.data(), N);
        initializeRandomDouble(y.data(), N);

        bench.start();
        for (long i=0; i < iters; ++i)
        { 
            y += a * x;
            y += b * x;
        }
        bench.stop();
    }

    bench.endImplementation();
}


void daxpyArrayVersion(BenchmarkExt<int>& bench, double a)
{
    bench.beginImplementation("Array<T,1>");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Array<T,1>: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Array<double,1> x(N), y(N);
        initializeRandomDouble(x.data(), N);
        initializeRandomDouble(y.data(), N);

        double b = - a;

        bench.start();
        for (long i=0; i < iters; ++i)
        {
            y = y + a * x;
            y = y + b * x;
	    sink();
        }
        bench.stop();

        bench.startOverhead();
        for (long i=0; i < iters; ++i) {
            sink();
	}
        bench.stopOverhead();
    }

    bench.endImplementation();
}

void daxpyF77Version(BenchmarkExt<int>& bench, double a)
{
    bench.beginImplementation("Fortran 77");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Fortran 77: N = " << N << endl;
        cout.flush();

        int iters = bench.getIterations();

        double* x = new double[N];
        double* y = new double[N];
        initializeRandomDouble(x, N);
        initializeRandomDouble(y, N);

        bench.start();
        fidaxpy(a, x, y, N, iters);
        bench.stop();

        delete [] x;
        delete [] y;
    }

    bench.endImplementation();
}


void daxpyBLASVersion(BenchmarkExt<int>& bench, double a)
{
#ifdef USE_LIBBLAS
    bench.beginImplementation("Platform BLAS");
#else
    bench.beginImplementation("Fortran BLAS");
#endif

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Fortran BLAS: N = " << N << endl;
        cout.flush();

        int iters = bench.getIterations();

        double* x = new double[N];
        double* y = new double[N];
        initializeRandomDouble(x, N);
        initializeRandomDouble(y, N);

        int xstride = 1, ystride = 1;
        double b = - a;

        bench.start();
        for (long i=0; i < iters; ++i)
        {
#ifdef USE_LIBBLAS
            daxpy(N, a, x, xstride, y, ystride);
            daxpy(N, b, x, xstride, y, ystride);
#else
            fdaxpy(N, a, x, xstride, y, ystride);
            fdaxpy(N, b, x, xstride, y, ystride);
#endif
        }
        bench.stop();

        delete [] x;
        delete [] y;
    }

    bench.endImplementation();
}

#ifdef FORTRAN_90
void daxpyF90Version(BenchmarkExt<int>& bench, double a)
{
    bench.beginImplementation("Fortran 90");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Fortran 90: N = " << N << endl;
        cout.flush();

        int iters = bench.getIterations();

        double* x = new double[N];
        double* y = new double[N];
        initializeRandomDouble(x, N);
        initializeRandomDouble(y, N);

        bench.start();
        f90daxpy(a, x, y, N, iters);
        bench.stop();

        delete [] x;
        delete [] y;
    }

    bench.endImplementation();
}
#endif

#ifdef BENCHMARK_VALARRAY
void daxpyValarrayVersion(BenchmarkExt<int>& bench, double a)
{
    bench.beginImplementation("valarray<T>");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "valarray<T>: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        valarray<double> x(N), y(N);
        initializeRandomDouble(x, N);
        initializeRandomDouble(y, N);

        double b = - a;

        bench.start();
        for (long i=0; i < iters; ++i)
        {
            y += a * x;
            y += b * x;
        }
        bench.stop();
    }

    bench.endImplementation();
}
#endif
