// Array DAXPY benchmark

#include <blitz/benchext.h>
#include <blitz/array.h>
#include <blitz/vector2.h>
#include <random/uniform.h>

BZ_NAMESPACE(blitz)
extern void sink();
BZ_NAMESPACE_END

BZ_USING_NAMESPACE(blitz)

#ifdef BZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES
  #define arrdaxpyf arrdaxpyf_
#elif defined(BZ_FORTRAN_SYMBOLS_WITH_DOUBLE_TRAILING_UNDERSCORES)
  #define arrdaxpyf arrdaxpyf__
#endif

extern "C" {
    void arrdaxpyf(double* A, double* B, int& N, double& a);
}

void arrdaxpyFortran77Version(BenchmarkExt<int>& bench);
void arrdaxpyBlitzVersion(BenchmarkExt<int>& bench);

int main()
{
    BenchmarkExt<int> bench("Array DAXPY", 2);

    const int numSizes = 8;

    bench.setNumParameters(numSizes);
    bench.setDependentVariable("flops");

    Vector<int> parameters(numSizes);
    Vector<long> iters(numSizes);
    Vector<double> flops(numSizes);

    parameters = pow(2.,tensor::i);
    cout << parameters;
    iters = 100*16*32*8*8*8/pow3(parameters);
    cout << iters;
    flops = pow3(parameters) * 2 * 2;
    cout << flops;

    bench.setParameterVector(parameters);
    bench.setParameterDescription("3D Array size");
    bench.setIterations(iters);
    bench.setOpsPerIteration(flops);

    bench.beginBenchmarking();
    arrdaxpyBlitzVersion(bench);
    arrdaxpyFortran77Version(bench);
    bench.endBenchmarking();

    bench.saveMatlabGraph("arrdaxpy.m");

    return 0;
}

void initializeRandomDouble(double* data, int numElements)
{
  ranlib::Uniform<double> rnd;

    for (int i=0; i < numElements; ++i)
        data[i] = rnd.random();
}

void arrdaxpyBlitzVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Blitz++");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Blitz++: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Array<double,3> A(N,N,N), B(N,N,N);
        initializeRandomDouble(A.data(), N*N*N);
        initializeRandomDouble(B.data(), N*N*N);
        TinyVector<int,2> size = N-2;
        double a = 0.34928313;
        double b = - a; 

        bench.start();
        for (long i=0; i < iters; ++i)
        {
            A += a * B;
            A += b * B;
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

void arrdaxpyFortran77Version(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Fortran 77");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Fortran 77: N = " << N << endl;
        cout.flush();

        int iters = (int)bench.getIterations();

        size_t arraySize = size_t(N) * size_t(N) * N;
       
        double* A = new double[arraySize];
        double* B = new double[arraySize];

        initializeRandomDouble(A, arraySize);
        initializeRandomDouble(B, arraySize);

        double a = 0.34928313;

        for (long i=0; i < iters; ++i)
        {
	  arrdaxpyf(A,B,N,a);
	  sink();
        }
        bench.stop();

        bench.startOverhead();
        for (long i=0; i < iters; ++i) {
            sink();
	}
        bench.stopOverhead();

        delete [] A;
        delete [] B;
    }

    bench.endImplementation();
}
