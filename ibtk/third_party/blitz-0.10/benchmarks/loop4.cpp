// In KAI C++ 3.2, restrict causes problems for copy propagation.
// Temporary fix: disable restrict

#define BZ_DISABLE_RESTRICT

#include <blitz/vector.h>
#include <blitz/array.h>
#include <blitz/rand-uniform.h>
#include <blitz/benchext.h>

#ifdef BZ_HAVE_VALARRAY
	#define BENCHMARK_VALARRAY
#endif

#ifdef BENCHMARK_VALARRAY
#include <valarray>
#endif

BZ_USING_NAMESPACE(blitz)

#if defined(BZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES)
 #define loop4_f77 loop4_f77_
 #define loop4_f90 loop4_f90_
#elif defined(BZ_FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES)
 #define loop4_f77 loop4_f77__
 #define loop4_f90 loop4_f90__
#elif defined(BZ_FORTRAN_SYMBOLS_CAPS)
 #define loop4_f77 LOOP4_F77
 #define loop4_f90 LOOP4_F90
#endif

extern "C" {
  void loop4_f77(const int& niters, const int& N, double* x, double* a, double* b);

  void loop4_f90(const int& niters, const int& N, double* x, double* a, double* b);

}

void VectorVersion(BenchmarkExt<int>& bench);
void ArrayVersion(BenchmarkExt<int>& bench);
void F77Version(BenchmarkExt<int>& bench);
void F90Version(BenchmarkExt<int>& bench);

#ifdef BENCHMARK_VALARRAY
void ValarrayVersion(BenchmarkExt<int>& bench);
#endif

int main()
{
#ifdef BENCHMARK_VALARRAY
    int numBenchmarks = 5;
#else
    int numBenchmarks = 4;
#endif

    BenchmarkExt<int> bench("$x=$a+$b", numBenchmarks);

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
        flops[i] = 1 * parameters[i];
    }

    bench.setParameterVector(parameters);
    bench.setIterations(iters);
    bench.setOpsPerIteration(flops);

    bench.beginBenchmarking();


    VectorVersion(bench);
    ArrayVersion(bench);
    F77Version(bench);
#ifdef FORTRAN_90
    F90Version(bench);
#endif
#ifdef BENCHMARK_VALARRAY
    ValarrayVersion(bench);
#endif

    bench.endBenchmarking();

    bench.saveMatlabGraph("loop4.m");

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

void VectorVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Vector<T>");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();
        cout << "Vector<T>: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Vector<double> x(N);
        initializeArray(x);
        Vector<double> a(N);
        initializeArray(a);
        Vector<double> b(N);
        initializeArray(b);

        bench.start();
        for (long i=0; i < iters; ++i)
        {
            x=a+b;
        }
        bench.stop();
    }

    bench.endImplementation();
}

void ArrayVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Array<T,1>");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();
        cout << "Array<T,1>: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Array<double> x(N);
        initializeArray(x);
        Array<double> a(N);
        initializeArray(a);
        Array<double> b(N);
        initializeArray(b);

        bench.start();
        for (long i=0; i < iters; ++i)
        {
            x=a+b;
        }
        bench.stop();
    }

    bench.endImplementation();
}

void F77Version(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Fortran 77");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Fortran 77: N = " << N << endl;
        cout.flush();

        int iters = bench.getIterations();

        double* x = new double[N];
        initializeRandomDouble(x, N);
        double* a = new double[N];
        initializeRandomDouble(a, N);
        double* b = new double[N];
        initializeRandomDouble(b, N);

        bench.start();
        loop4_f77(iters, N, x, a, b);
        bench.stop();

        delete [] x;
        delete [] a;
        delete [] b;
    }

    bench.endImplementation();
}

#ifdef FORTRAN_90
void F90Version(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Fortran 90");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Fortran 90: N = " << N << endl;
        cout.flush();

        int iters = bench.getIterations();

        double* x = new double[N];
        initializeRandomDouble(x, N);
        double* a = new double[N];
        initializeRandomDouble(a, N);
        double* b = new double[N];
        initializeRandomDouble(b, N);

        bench.start();
        loop4_f90(iters, N, x, a, b);
        bench.stop();

        delete [] x;
        delete [] a;
        delete [] b;
    }

    bench.endImplementation();
}
#endif
