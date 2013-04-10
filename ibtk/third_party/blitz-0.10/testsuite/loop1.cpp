/*
 * This benchmark causes some versions of gcc to crash.
 */

#include <blitz/vector2.h>
#include <blitz/array.h>
#include <random/uniform.h>
#include <blitz/benchext.h>

BZ_NAMESPACE(blitz)
extern void sink();
BZ_NAMESPACE_END

BZ_USING_NAMESPACE(blitz)

void VectorVersion(BenchmarkExt<int>& bench);
void ArrayVersion(BenchmarkExt<int>& bench);

int dontActuallyRunBenchmark()
{
    return 1;
}

int main()
{
    if (dontActuallyRunBenchmark())
        return 0;

    int numBenchmarks = 2;

    BenchmarkExt<int> bench("loop1: $x=sqrt($y)", numBenchmarks);

    const int numSizes = 26;
    bench.setNumParameters(numSizes);
    //bench.setRateDescription("Mflops/s"); removed

    Vector<int> parameters(numSizes);
    Vector<long> iters(numSizes);
    Vector<double> flops(numSizes);

    for (int i=0; i < numSizes; ++i)
    {
      parameters(i) = (int)pow(10.0, (i+1)/4.0);
      iters[i] = 50000000L / parameters(i);
        if (iters(i) < 2)
	  iters(i) = 2;
        flops(i) = 1 * parameters(i);
    }

    bench.setParameterVector(parameters);
    bench.setIterations(iters);
    bench.setOpsPerIteration(flops);

    bench.beginBenchmarking();


    VectorVersion(bench);
    ArrayVersion(bench);

    bench.endBenchmarking();

    bench.saveMatlabGraph("loop1.m");

    return 0;
}

template<typename T>
void initializeRandomDouble(T* data, int numElements, int stride = 1)
{
  ranlib::Uniform<T> rnd;

    for (int i=0; i < numElements; ++i)
        data[size_t(i*stride)] = rnd.random();
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
        initializeRandomDouble(x.data(), N);
        Vector<double> y(N);
        initializeRandomDouble(y.data(), N);

        bench.start();
        long i;
        for (i=0; i < iters; ++i)
        {
            x=sqrt(y);
            sink();
        }
        bench.stop();

        bench.startOverhead();
        for (i=0; i < iters; ++i)
            sink();
        bench.stopOverhead();
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

        Array<double, 1> x(N);
        initializeRandomDouble(x.dataFirst(), N);
        Array<double, 1> y(N);
        initializeRandomDouble(y.dataFirst(), N);

        bench.start();
        long i;
        for (i=0; i < iters; ++i)
        {
            x=sqrt(y);
            sink();
        }
        bench.stop();

        bench.startOverhead();
        for (i=0; i < iters; ++i)
            sink();
        bench.stopOverhead();
    }

    bench.endImplementation();
}

