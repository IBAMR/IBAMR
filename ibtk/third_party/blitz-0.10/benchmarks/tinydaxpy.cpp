// TinyVector<T,N> DAXPY benchmark

//#define BZ_DISABLE_KCC_COPY_PROPAGATION_KLUDGE

#include <blitz/array.h>
#include <blitz/timer.h>
#include <random/uniform.h>

BZ_USING_NAMESPACE(blitz)

ranlib::Uniform<double> rnd;

template<class T>
void optimizationSink(T&);

template<int N_rank>
void tinyDAXPYBenchmark(TinyVector<double,N_rank>, int iters, double a)
{
    Timer timer;
   
    TinyVector<double,N_rank> ta, tb, tc, td, te, tf, tg, th, ti, tj;
    for (int i=0; i < N_rank; ++i)
    {
        ta[i] = rnd.random()+1;
        tb[i] = rnd.random()+1;
        tc[i] = rnd.random()+1;
        td[i] = rnd.random()+1;
        te[i] = rnd.random()+1;
        tf[i] = rnd.random()+1;
        tg[i] = rnd.random()+1;
        th[i] = rnd.random()+1;
        ti[i] = rnd.random()+1;
        tj[i] = rnd.random()+1;
    }

    double b = -a;

    double numFlops = 0;

    if (N_rank < 20)
    {
      timer.start();
      for (int i=0; i < iters; ++i)
      {
        ta += a * tb;
        tc += a * td;
        te += a * tf;
        tg += a * th;
        ti += a * tj;
        tb += b * ta;
        td += b * tc;
        tf += b * te;
        th += b * tg;
        tj += b * ti;
        ta += a * tb;
        tc += a * td;
        te += a * tf;
        tg += a * th;
        ti += a * tj;
        tb += b * ta;
        td += b * tc;
        tf += b * te;
        th += b * tg;
        tj += b * ti;
      }
      timer.stop();
      numFlops = 40.0 * N_rank * double(iters);
    }
    else {
      timer.start();
      for (int i=0; i < iters; ++i)
      {
        ta += a * tb;
        tb += b * ta;
      }
      timer.stop();
      numFlops = 4.0 * N_rank * double(iters);
    }

    optimizationSink(ta);
    optimizationSink(tb);
    optimizationSink(tc);
    optimizationSink(td);
    optimizationSink(te);
    optimizationSink(tf);
    optimizationSink(tg);
    optimizationSink(th);
    optimizationSink(ti);
    optimizationSink(tj);

    timer.stop();
    float Gflops = numFlops / (1e9*timer.elapsed());

    if (iters > 1)  
    {
    cout << setw(5) << N_rank << '\t' << Gflops << endl;
    }
}

double a = 0.3429843;

template<class T>
void optimizationSink(T&)
{
}

int main()
{
    cout << "TinyVector<double,N> DAXPY benchmark" << endl
         << setw(5) << "N" << '\t' << "Gflops/" << Timer::indep_var() << endl;
    tinyDAXPYBenchmark(TinyVector<double,1>(), 800000, a);
    tinyDAXPYBenchmark(TinyVector<double,2>(), 800000, a);
    tinyDAXPYBenchmark(TinyVector<double,3>(), 800000, a);
    tinyDAXPYBenchmark(TinyVector<double,4>(), 700000, a);
    tinyDAXPYBenchmark(TinyVector<double,5>(), 600000, a);
    tinyDAXPYBenchmark(TinyVector<double,6>(), 500000, a);
    tinyDAXPYBenchmark(TinyVector<double,7>(), 500000, a);
    tinyDAXPYBenchmark(TinyVector<double,8>(), 500000, a);
    tinyDAXPYBenchmark(TinyVector<double,9>(), 500000, a);
    tinyDAXPYBenchmark(TinyVector<double,10>(), 500000, a);

    return 0;
}

