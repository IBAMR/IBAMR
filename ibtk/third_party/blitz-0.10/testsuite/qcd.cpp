// This program causes/caused some versions of egcs to crash.

#include <blitz/benchext.h>
#include <random/uniform.h>
#include <blitz/array.h>

#ifdef BZ_HAVE_COMPLEX

void initializeRandomDouble(double* data, int numElements);

BZ_USING_NAMESPACE(blitz)

int QCDBlitzVersion(BenchmarkExt<int>& bench);
int QCDBlitzTunedVersion(BenchmarkExt<int>& bench);

int dontActuallyRunBenchmark()
{
    return 1;
}

int main()
{
    if (dontActuallyRunBenchmark())
        return 0;

    cout << "Blitz++ QCD Benchmark" << endl
         << "Working... (this may take a while)  ";
    cout.flush();

    BenchmarkExt<int> bench("Lattice QCD Benchmark", 2);

    //bench.setRateDescription("Millions of operations/s"); \todo removed for now
    bench.beginBenchmarking();

    QCDBlitzVersion(bench);
    QCDBlitzTunedVersion(bench);

    bench.endBenchmarking();

    bench.saveMatlabGraph("qcd.m");

    cout << "Done." << endl;

    return 0;
}

int QCDBlitzVersion(BenchmarkExt<int>& bench)
{
    typedef TinyMatrix<complex<double>, 3, 2> spinor;
    typedef TinyMatrix<complex<double>, 3, 3> SU3Gauge;

    bench.beginImplementation("Blitz++");

    while (!bench.doneImplementationBenchmark())
    {
        int length = bench.getParameter();
        int iters = (int)bench.getIterations();

        Array<spinor,1> res(length), src(length);
        Array<SU3Gauge,1> M(length);

        initializeRandomDouble((double*)src.data(), 
            length * sizeof(spinor) / sizeof(double));
        initializeRandomDouble((double*)M.data(),
            length * sizeof(SU3Gauge) / sizeof(double));

        bench.start();
        long i;
        for (i=0; i < iters; ++i)
        {
            for (int i=0; i < length; ++i)
	      res(i) = M(i)*src(i);
        }
        bench.stop();

        // Time overhead
        bench.startOverhead();
        for (i=0; i < iters; ++i)
        {
        }
        bench.stopOverhead();
    }

    bench.endImplementation();

    return 0;
}

    typedef TinyMatrix<complex<double>, 3, 2> spinor;
    typedef TinyMatrix<complex<double>, 3, 3> gaugeFieldElement;

    struct latticeUnit {
        spinor one;
        gaugeFieldElement gauge;
        spinor two;
    };

int QCDBlitzTunedVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Blitz++ (tuned)");

    while (!bench.doneImplementationBenchmark())
    {
        int length = bench.getParameter();
        int iters = (int)bench.getIterations();

        Array<latticeUnit,1> lattice(length);

        initializeRandomDouble((double*)lattice.data(),
            length * sizeof(latticeUnit) / sizeof(double));

        bench.start();
        long i;
        for (i=0; i < iters; ++i)
        {
            for (int i=0; i < length; ++i)
	      lattice(i).two = lattice(i).gauge * lattice(i).two;
        }
        bench.stop();

        // Time overhead
        bench.startOverhead();
        for (i=0; i < iters; ++i)
        {
        }
        bench.stopOverhead();
    }

    bench.endImplementation();

    return 0;
}

void initializeRandomDouble(double* data, int numElements)
{
    // This is a temporary kludge until I implement random complex
    // numbers.

  ranlib::Uniform<double> rnd;

    for (int i=0; i < numElements; ++i)
        data[i] = rnd.random();
}

#else // BZ_HAVE_COMPLEX

#include <blitz/blitz.h>

int main()
{
    cout << "This benchmark requires <complex> from the ISO/ANSI C++ standard."
         << endl;
    return 0;
}

#endif // BZ_HAVE_COMPLEX
