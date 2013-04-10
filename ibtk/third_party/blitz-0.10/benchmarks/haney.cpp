// Haney's induction calculation benchmark.
//
// See: Scott W. Haney, Is C++ Fast Enough for Scientific Computing?
//      Computers in Physics Vol. 8 No. 6 (1994), p. 690
//
//      Arch D. Robison, C++ Gets Faster for Scientific Computing,
//      Computers in Physics Vol. 10 No. 5 (1996), p. 458
//

#include <blitz/vector.h>
#include <blitz/rand-uniform.h>
#include <blitz/benchext.h>
#ifdef BZ_HAVE_STD
#include <valarray>
#else
#include <valarray.h>
#endif

BZ_USING_NAMESPACE(blitz)

#ifndef M_PI
 #define M_PI   3.14159265358979323846
#endif

#ifdef BZ_FORTRAN_SYMBOLS_WITH_TRAILING_UNDERSCORES
#define vecopsf    vecopsf_
#define vecopsfo   vecopsfo_
#endif

extern "C"
{
	void vecopsf(float *li, const float *R, const float *w, const int &N,
	             const int& iters);
	void vecopsfo(float *li, const float *R, const float *w, const int &N,
	              const int& iters);
}

inline float sqr(float x)
{
	return x*x;
}

const float Mu0 = 4.0 * M_PI * 1.0e-7;

void HaneyCVersion(BenchmarkExt<int>& bench);
void HaneyFortranVersion(BenchmarkExt<int>& bench);
void HaneyBlitzVersion(BenchmarkExt<int>& bench);

int main()
{
	BenchmarkExt<int> bench("Haney Inductance Calculation", 3);

	bench.setRateDescription("Operations/s");

	bench.beginBenchmarking();

	HaneyCVersion(bench);
	HaneyFortranVersion(bench);
	HaneyBlitzVersion(bench);

	bench.endBenchmarking();

	bench.saveMatlabGraph("haney.m");

	return 0;
}

void initializeRandom(float* data, int length)
{
	Random<Uniform> unif(1.0, 2.0);
	for (int i=0; i < length; ++i)
		data[i] = unif.random();
}

void HaneyCVersion(BenchmarkExt<int>& bench)
{
	bench.beginImplementation("Inlined C");

	while (!bench.doneImplementationBenchmark()) {
		int length = bench.getParameter();
		long iters = bench.getIterations();

		cout << "length = " << length << " iters = " << iters << endl;

		float* li = new float[length];
		float* R = new float[length];
		float* w = new float[length];

		initializeRandom(li, length);
		initializeRandom(R, length);
		initializeRandom(w, length);

		// Tickle the cache
		for (int i=0; i < length; ++i)
			li[i] = R[i] + log(w[i]);

		bench.start();

		for (long j=0; j < iters; ++j) {
			for (int i=0; i < length; ++i) {
				li[i] = Mu0 * R[i] *
				        (0.5 * (1.0 + (1.0/24.0)
				                * sqr(w[i]/R[i])) * log(32.0 * sqr(R[i]/w[i]))
				         + 0.05 * sqr(w[i]/R[i]) - 0.85);
			}
		}

		bench.stop();

		// Subtract the loop overhead
		bench.startOverhead();

		for (long j=0; j < iters; ++j) {}



		bench.stopOverhead();

		delete [] li;
		delete [] w;
		delete [] R;
	}

	bench.endImplementation();
}

void HaneyFortranVersion(BenchmarkExt<int>& bench)
{
	bench.beginImplementation("Fortran");

	while (!bench.doneImplementationBenchmark()) {
		int length = bench.getParameter();
		int iters = (int)bench.getIterations();

		cout << "length = " << length << " iters = " << iters << endl;

		float* li = new float[length];
		float* R = new float[length];
		float* w = new float[length];

		initializeRandom(li, length);
		initializeRandom(R, length);
		initializeRandom(w, length);

		// Tickle
		int oneIter = 1;
		vecopsf(li, R, w, length, oneIter);

		// Time
		bench.start();
		vecopsf(li, R, w, length, iters);
		bench.stop();

		// Time overhead
		bench.startOverhead();
		vecopsfo(li, R, w, length, iters);
		bench.stopOverhead();

		delete [] li;
		delete [] w;
		delete [] R;
	}

	bench.endImplementation();
}

void HaneyBlitzVersion(BenchmarkExt<int>& bench)
{
	bench.beginImplementation("Blitz++");

	while (!bench.doneImplementationBenchmark()) {
		int length = bench.getParameter();
		int iters = (int)bench.getIterations();

		Vector<float> li(length), R(length), w(length);
		initializeRandom(li.data(), length);
		initializeRandom(R.data(), length);
		initializeRandom(w.data(), length);

		cout << "length = " << length << " iters = " << iters << endl;

		// Tickle
		li = w + log(R);

		// Time
		bench.start();
		for (long i=0; i < iters; ++i) {
#if defined(__GNUC__) && (__GNUC__ < 3)
			li = Mu0 * R * ( (0.5 + (0.5/24.0) * sqr(w/R) ) 
			                 * log(32.0 * sqr(R/w)) + 0.05 * sqr(w/R) - 0.85);
#else
			li = Mu0 * R * (0.5 * (1.0 + (1.0/24.0) * sqr(w/R))
			                * log(32.0 * sqr(R/w)) + 0.05 * sqr(w/R) - 0.85);
#endif
		}
		bench.stop();

		// Time overhead
		bench.startOverhead();
		for (long i=0; i < iters; ++i) {
		}
		bench.stopOverhead();
	}

	bench.endImplementation();
}

