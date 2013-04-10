// Array expression benchmark

#include <blitz/array.h>
#include <blitz/benchext.h>

BZ_USING_NAMESPACE(blitz)

void blitzVersion(BenchmarkExt<int>& bench);
void CVersion(BenchmarkExt<int>& bench);

int main()
{
    BenchmarkExt<int> bench("Array expression", 2);

    bench.beginBenchmarking();
    blitzVersion(bench);
    CVersion(bench);
    bench.endBenchmarking();

    bench.saveMatlabGraph("arrexpr1.m");

    return 0;
}

void blitzVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Blitz++");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Blitz++: N = " << N << endl;

        long iters = bench.getIterations();

        Array<double,1> x(N);

        // Tickle
        x = 0.;
        firstIndex i;

        bench.start();
        for (long it=0; it < iters; ++it)
        {
            x = i * i;
        }
        bench.stop();
    }

    bench.endImplementation();
}

void CVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("C");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "C: N = " << N << endl;       

        long iters = bench.getIterations();

        double* x = new double[N];

        // Tickle
        for (int i=0; i < N; ++i)
            x[i] = 0;

        bench.start();
        for (long it=0; it < iters; ++it)
        {
            for (int i=0; i < N; ++i)
                x[i] = i * i;
        }
        bench.stop();

        delete [] x;
    }

    bench.endImplementation();
}

