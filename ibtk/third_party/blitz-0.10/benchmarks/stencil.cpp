// Array stencil benchmark

#include <blitz/array.h>
#include <blitz/traversal.h>
#include <blitz/benchext.h>
#include <blitz/rand-uniform.h>
#include <blitz/array/stencil-et.h>

BZ_USING_NAMESPACE(blitz)

#if defined(BZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES)
  #define stencilf stencilf_
  #define stencilftiled stencilftiled_
  #define stencilf90 stencilf90_
#elif defined(BZ_FORTRAN_SYMBOLS_CAPS)
  #define stencilf       STENCILF
  #define stencilftiled  STENCILFTILED
  #define stencilf90     STENCILF90
#endif

extern "C" {
    void stencilf(double* A, double* B, int& N, int& iters);
    void stencilftiled(double* A, double* B, int& N, int& iters);
    void stencilf90(double* A, double* B, int& N, int& iters);
}

#ifdef FORTRAN_90
void stencilFortran90Version(BenchmarkExt<int>& bench);
#endif
void stencilFortran77Version(BenchmarkExt<int>& bench);
void stencilFortran77VersionTiled(BenchmarkExt<int>& bench);
void stencilBlitzVersion(BenchmarkExt<int>& bench);
void stencilBlitzExpressionVersion(BenchmarkExt<int>& bench);
void stencilBlitzProductVersion(BenchmarkExt<int>& bench);
void stencilBlitzProductVersion2(BenchmarkExt<int>& bench);
void stencilBlitzProductVersion3(BenchmarkExt<int>& bench);
void stencilBlitzStencilVersion(BenchmarkExt<int>& bench);
void stencilBlitzIndexVersion(BenchmarkExt<int>& bench);

int main()
{
    int numBenchmarks = 10;
#ifndef FORTRAN_90
		numBenchmarks--;   // No fortran 90
#endif

    BenchmarkExt<int> bench("Array stencil", numBenchmarks);

    const int numSizes = 28;

    bench.setNumParameters(numSizes);
    bench.setRateDescription("Mflops/s");

    Vector<int> parameters(numSizes);
    Vector<long> iters(numSizes);
    Vector<double> flops(numSizes);

    for (int i=0; i < numSizes; ++i)
    {
        parameters[i] = (i+1) * 8;
        iters[i] = 32*8*8*8/(i+1)/(i+1)/(i+1)/4;
        if (iters[i] < 2)
            iters[i] = 2;
        int npoints = parameters[i] - 2;
        flops[i] = npoints * npoints * npoints * 7 * 2;
    }

    bench.setParameterVector(parameters);
    bench.setIterations(iters);
    bench.setOpsPerIteration(flops);

    bench.beginBenchmarking();
#ifdef FORTRAN_90
    stencilFortran90Version(bench);
#endif
    stencilBlitzVersion(bench);
    stencilBlitzStencilVersion(bench);
    stencilBlitzExpressionVersion(bench);
    stencilBlitzProductVersion(bench);
    stencilBlitzProductVersion2(bench);
    stencilBlitzProductVersion3(bench);
    stencilBlitzIndexVersion(bench);
    stencilFortran77Version(bench);
    stencilFortran77VersionTiled(bench);
    bench.endBenchmarking();

    bench.saveMatlabGraph("stencil.m","plot");

    return 0;
}

void initializeRandomDouble(double* data, int numElements, int stride = 1)
{
    static Random<Uniform> rnd;

    for (int i=0; i < numElements; ++i)
        data[i*stride] = rnd.random();
}

void stencilBlitzVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Blitz++ Range Expr");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Blitz++: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Array<double,3> A(N,N,N), B(N,N,N);
        initializeRandomDouble(A.data(), N*N*N, A.stride(thirdDim));
        initializeRandomDouble(B.data(), N*N*N, B.stride(thirdDim));
        TinyVector<int,2> size = N-2;
        generateFastTraversalOrder(size);
        double c = 1/7.;
       
        bench.start();
        for (long i=0; i < iters; ++i)
        {
            Range I(1,N-2), J(1,N-2), K(1,N-2);

            A(I,J,K) = c * (B(I,J,K) + B(I+1,J,K) + B(I-1,J,K) + B(I,J+1,K)
                + B(I,J-1,K) + B(I,J,K+1) + B(I,J,K-1));

            B(I,J,K) = c * (A(I,J,K) + A(I+1,J,K) + A(I-1,J,K) + A(I,J+1,K)
                + A(I,J-1,K) + A(I,J,K+1) + A(I,J,K-1));
        }
        bench.stop();
    }

    bench.endImplementation();
}

BZ_DECLARE_STENCIL_OPERATOR1(test1,B)
return (1./7) * ( (*B) + B.shift(1,0) + B.shift(-1,0) + B.shift(1,1)
		  + B.shift(-1,1) + B.shift(1,2) + B.shift(-1,2));
BZ_END_STENCIL_OPERATOR

BZ_ET_STENCIL(test1, double, double,shape(-1,-1,-1),shape(1,1,1))

BZ_DECLARE_STENCIL2(test1stencil,A,B)
  A=test1(B);
BZ_END_STENCIL_WITH_SHAPE(shape(-1,-1,-1),shape(1,1,1))

void stencilBlitzExpressionVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Blitz++ StencilOp");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Blitz++ Stencil Operator: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Array<double,3> A(N,N,N), B(N,N,N);
        initializeRandomDouble(A.data(), N*N*N, A.stride(thirdDim));
        initializeRandomDouble(B.data(), N*N*N, B.stride(thirdDim));
        TinyVector<int,2> size = N-2;
        generateFastTraversalOrder(size);
        double c = 1/7.;
       
        bench.start();
        for (long i=0; i < iters; ++i)
        {
            Range I(1,N-2), J(1,N-2), K(1,N-2);
	    A(I,J,K) = test1(B);

	    B(I,J,K) = test1(A);
        }
        bench.stop();
    }

    bench.endImplementation();
}

void stencilBlitzProductVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Blitz++ product");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Blitz++ Stencil Operator on product: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Array<double,3> A(N,N,N), B(N,N,N),C(N,N,N);
        initializeRandomDouble(A.data(), N*N*N, A.stride(thirdDim));
        initializeRandomDouble(B.data(), N*N*N, B.stride(thirdDim));
        TinyVector<int,2> size = N-2;
        generateFastTraversalOrder(size);
        double c = 1/7.;
       
        bench.start();
        for (long i=0; i < iters; ++i)
        {
            Range I(1,N-2), J(1,N-2), K(1,N-2);
	    C=B*B;
	    A(I,J,K) = test1(C);
	    C=A*A;
	    B(I,J,K) = test1(C);
        }
        bench.stop();
    }

    bench.endImplementation();
}

void stencilBlitzProductVersion2(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Blitz++ product w alloc");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Blitz++ Stencil Operator on product: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Array<double,3> A(N,N,N), B(N,N,N);
        initializeRandomDouble(A.data(), N*N*N, A.stride(thirdDim));
        initializeRandomDouble(B.data(), N*N*N, B.stride(thirdDim));
        TinyVector<int,2> size = N-2;
        generateFastTraversalOrder(size);
        double c = 1/7.;
       
        bench.start();
        for (long i=0; i < iters; ++i)
        {
            Range I(1,N-2), J(1,N-2), K(1,N-2);
	    {
	      Array<double,3>C(B*B);
	      A(I,J,K) = test1(C);
	    }
	    {
	      Array<double,3>C(A*A);
	      B(I,J,K) = test1(C);
	    }
        }
        bench.stop();
    }

    bench.endImplementation();
}

void stencilBlitzProductVersion3(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Blitz++ product expr");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Blitz++ Stencil Operator on product expr: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Array<double,3> A(N,N,N), B(N,N,N);
        initializeRandomDouble(A.data(), N*N*N, A.stride(thirdDim));
        initializeRandomDouble(B.data(), N*N*N, B.stride(thirdDim));
        TinyVector<int,2> size = N-2;
        generateFastTraversalOrder(size);
        double c = 1/7.;
       
        bench.start();
        for (long i=0; i < iters; ++i)
        {
            Range I(1,N-2), J(1,N-2), K(1,N-2);
	    A(I,J,K) = test1(B*B);

	    B(I,J,K) = test1(A*A);
        }
        bench.stop();
    }

    bench.endImplementation();
}

void stencilBlitzIndexVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Blitz++ Indexed StencilOp");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Blitz++ Indexed Stencil Operator: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Array<double,3> A(N,N,N), B(N,N,N);
        initializeRandomDouble(A.data(), N*N*N, A.stride(thirdDim));
        initializeRandomDouble(B.data(), N*N*N, B.stride(thirdDim));
        TinyVector<int,2> size = N-2;
        generateFastTraversalOrder(size);
        double c = 1/7.;
       
        bench.start();
        for (long i=0; i < iters; ++i)
        {
            Range I(1,N-2), J(1,N-2), K(1,N-2);
	    A(I,J,K) = test1(B(tensor::i, tensor::j, tensor::k));

	    B(I,J,K) = test1(A(tensor::i, tensor::j, tensor::k));
        }
        bench.stop();
    }

    bench.endImplementation();
}

void stencilBlitzStencilVersion(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Blitz++ Stencil");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Blitz++ Stencil: N = " << N << endl;
        cout.flush();

        long iters = bench.getIterations();

        Array<double,3> A(N,N,N), B(N,N,N);
        initializeRandomDouble(A.data(), N*N*N, A.stride(thirdDim));
        initializeRandomDouble(B.data(), N*N*N, B.stride(thirdDim));
        TinyVector<int,2> size = N-2;
        generateFastTraversalOrder(size);
        double c = 1/7.;
       
	;        bench.start();
        for (long i=0; i < iters; ++i)
        {
            Range I(1,N-2), J(1,N-2), K(1,N-2);
	    applyStencil(test1stencil(),A,B);
	    applyStencil(test1stencil(),B,A);
        }
        bench.stop();
    }

    bench.endImplementation();
}

void stencilFortran77Version(BenchmarkExt<int>& bench)
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

        bench.start();
        stencilf(A, B, N, iters);
        bench.stop();

        delete [] A;
        delete [] B;
    }

    bench.endImplementation();
}

void stencilFortran77VersionTiled(BenchmarkExt<int>& bench)
{
    bench.beginImplementation("Fortran 77 (tiled)");

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

        bench.start();
        stencilftiled(A, B, N, iters);
        bench.stop();

        delete [] A;
        delete [] B;
    }

    bench.endImplementation();
}

#ifdef FORTRAN_90
void stencilFortran90Version(BenchmarkExt<int>& bench)
{
   bench.beginImplementation("Fortran 90");

    while (!bench.doneImplementationBenchmark())
    {
        int N = bench.getParameter();

        cout << "Fortran 90: N = " << N << endl;
        cout.flush();

        int iters = (int)bench.getIterations();

        size_t arraySize = size_t(N) * size_t(N) * N;

        double* A = new double[arraySize];
        double* B = new double[arraySize];

        initializeRandomDouble(A, arraySize);
        initializeRandomDouble(B, arraySize);

        bench.start();
        stencilf90(A, B, N, iters);
        bench.stop();

        delete [] A;
        delete [] B;
    }

    bench.endImplementation();
}
#endif
