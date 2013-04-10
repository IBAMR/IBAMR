#include <blitz/tinymat.h>
#include <blitz/vector.h>
#include <blitz/benchext.h>
#include <blitz/rand-uniform.h>

#ifdef BZ_HAVE_COMPLEX

BZ_USING_NAMESPACE(blitz)

#if defined(BZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES)
  #define qcdf  qcdf_
  #define qcdf2 qcdf2_
#elif defined( BZ_FORTRAN_SYMBOLS_CAPS)
  #define qcdf  QCDF
  #define qcdf2 QCDF2
#endif

extern "C" {
  void qcdf(const void* M, void* res, const void* src, const int& N,
    const int& iters);
  void qcdf2(const void* M, void* res, const void* src, const int& N,
    const int& iters);
}

int QCDBlitzVersion(BenchmarkExt<int>& bench);
int QCDBlitzTunedVersion(BenchmarkExt<int>& bench);
int QCDFortran77Version(BenchmarkExt<int>& bench);
int QCDFortran77TunedVersion(BenchmarkExt<int>& bench);

void initializeRandomDouble(double* data, int numElements);

int main()
{
    cout << "Blitz++ QCD Benchmark" << endl
         << "Working... (this may take a while)  ";
    cout.flush();

    BenchmarkExt<int> bench("Lattice QCD Benchmark", 4);

    bench.setRateDescription("Millions of operations/s");
    bench.beginBenchmarking();

    QCDBlitzVersion(bench);
    QCDBlitzTunedVersion(bench);
    QCDFortran77Version(bench);
    QCDFortran77TunedVersion(bench);

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

cout << "length = " << length << endl;

        Vector<spinor> res(length), src(length);
        Vector<SU3Gauge> M(length);

        initializeRandomDouble((double*)src.data(), 
            length * sizeof(spinor) / sizeof(double));
        initializeRandomDouble((double*)M.data(),
            length * sizeof(SU3Gauge) / sizeof(double));

        bench.start();
        for (long i=0; i < iters; ++i)
        {
            for (int i=0; i < length; ++i)
                res[i] = product(M[i], src[i]);
        }
        bench.stop();

        // Time overhead
        bench.startOverhead();
        for (long i=0; i < iters; ++i)
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

        Vector<latticeUnit> lattice(length);

        initializeRandomDouble((double*)lattice.data(),
            length * sizeof(latticeUnit) / sizeof(double));

        bench.start();
        for (long i=0; i < iters; ++i)
        {
            for (int i=0; i < length; ++i)
                lattice[i].two = product(lattice[i].gauge, lattice[i].one);
        }
        bench.stop();

        // Time overhead
        bench.startOverhead();
        for (long i=0; i < iters; ++i)
        {
        }
        bench.stopOverhead();
    }

    bench.endImplementation();

    return 0;
}

int QCDFortran77Version(BenchmarkExt<int>& bench)
{
    // Use Blitz++ library only to allocate space for the
    // arrays. 
    typedef TinyMatrix<complex<double>, 3, 2> spinor;
    typedef TinyMatrix<complex<double>, 3, 3> SU3Gauge;

    bench.beginImplementation("Fortran 77");

    while (!bench.doneImplementationBenchmark())
    {
        int length = bench.getParameter();
        int iters = (int)bench.getIterations();

        Vector<spinor> res(length), src(length);
        Vector<SU3Gauge> M(length);

        initializeRandomDouble((double*)src.data(),
            length * sizeof(spinor) / sizeof(double));
        initializeRandomDouble((double*)M.data(),
            length * sizeof(SU3Gauge) / sizeof(double));

        bench.start();
        qcdf(M.data(), res.data(), src.data(), length, iters);
        bench.stop();

        // Time overhead
        bench.startOverhead();
        for (long i=0; i < iters; ++i)
        {
        }
        bench.stopOverhead();
    }

    bench.endImplementation();

    return 0;
}

int QCDFortran77TunedVersion(BenchmarkExt<int>& bench)
{
    // Use Blitz++ library only to allocate space for the
    // arrays.
    typedef TinyMatrix<complex<double>, 3, 2> spinor;
    typedef TinyMatrix<complex<double>, 3, 3> SU3Gauge;

    bench.beginImplementation("Fortran 77 Hand-tuned");

    while (!bench.doneImplementationBenchmark())
    {
        int length = bench.getParameter();
        int iters = (int)bench.getIterations();

        Vector<spinor> res(length), src(length);
        Vector<SU3Gauge> M(length);

        initializeRandomDouble((double*)src.data(),
            length * sizeof(spinor) / sizeof(double));
        initializeRandomDouble((double*)M.data(),
            length * sizeof(SU3Gauge) / sizeof(double));

        bench.start();
        qcdf2(M.data(), res.data(), src.data(), length, iters);
        bench.stop();

        // Time overhead
        bench.startOverhead();
        for (long i=0; i < iters; ++i)
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

    static Random<Uniform> rnd;

    for (int i=0; i < numElements; ++i)
        data[i] = rnd.random();
}

#else // BZ_HAVE_COMPLEX

#include <iostream.h>

int main()
{
    cout << "This benchmark requires <complex> from the ISO/ANSI C++ standard."
         << endl;
    return 0;
}

#endif // BZ_HAVE_COMPLEX
