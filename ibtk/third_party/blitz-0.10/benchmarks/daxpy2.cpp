// In KAI C++ 3.2c, restrict causes problems for copy propagation.
// Temporary kludge is to disable use of the restrict keyword.

#define BZ_DISABLE_RESTRICT

#include <blitz/vector.h>
#include <blitz/array.h>
#include <blitz/rand-uniform.h>
#include <blitz/benchext.h>

#ifdef BENCHMARK_VALARRAY
#include <valarray>
#endif

BZ_USING_NAMESPACE(blitz)

#ifdef BZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES
 #define fdaxpy   fdaxpy_
 #define f90daxpy f90daxpy_
 #define fidaxpy  fidaxpy_
 #define fidaxpyo fidaxpyo_
#endif

#ifdef BZ_FORTRAN_SYMBOLS_CAPS
 #define fdaxpy   FDAXPY
 #define f90daxpy F90DAXPY
 #define fidaxpy  FIDAXPY
 #define fidaxpyo FIDAXPYO
#endif

extern "C" {
  void fdaxpy(const int& N, const double& da, double* x,
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
void daxpyF90Version(BenchmarkExt<int>& bench, double a);

#ifdef BENCHMARK_VALARRAY
void daxpyValarrayVersion(BenchmarkExt<int>& bench, double a);
#endif

int main()
{

#ifdef BENCHMARK_VALARRAY
    int numBenchmarks = 6;
#else
    int numBenchmarks = 5;
#endif

    BenchmarkExt<int> bench("DAXPY Benchmark", numBenchmarks);

    const int numSizes = 19;
    bench.setNumParameters(numSizes);
    bench.setRateDescription("Mflops/s");

    Vector<int> parameters(numSizes);
    Vector<long> iters(numSizes);
    Vector<double> flops(numSizes);

    for (int i=0; i < numSizes; ++i)
    {
        parameters[i] = pow(10.0, (i+1)/4.0);
        iters[i] = 50000000L / parameters[i];
        if (iters[i] < 2)
            iters[i] = 2;
        flops[i] = 2 * parameters[i] * 2;
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
    daxpyF90Version(bench, a);

#ifdef BENCHMARK_VALARRAY
    daxpyValarrayVersion(bench, a);
#endif

    bench.endBenchmarking();

    bench.saveMatlabGraph("daxpy2.m");

    return 0;
}

void initializeRandomDouble(double* data, int numElements, int stride = 1)
{
    static Random<Uniform> rnd;

    for (int i=0; i < numElements; ++i)
        data[i*stride] = rnd.random();
}

template<class T>
void initializeArray(T& array, int numElements)
{
    static Random<Uniform> rnd;

    for (size_t i=0; i < numElements; ++i)
        array[i] = rnd.random();
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
            y += a * x;
            y += b * x;
        }
        bench.stop();
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
    bench.beginImplementation("Fortran BLAS");

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
            fdaxpy(N, a, x, xstride, y, ystride);
            fdaxpy(N, b, x, xstride, y, ystride);
        }
        bench.stop();

        delete [] x;
        delete [] y;
    }

    bench.endImplementation();
}

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
        initializeArray(x, N);
        initializeArray(y, N);

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
