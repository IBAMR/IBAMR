#include <blitz/tinyvec.h>
#include <blitz/tinymat.h>
#include <blitz/timer.h>

using namespace blitz;

template<class T>
void optimizationSink(T&);

int main()
{
    TinyMatrix<float,3,3> A1, A2;
    TinyVector<float,3> b1, b2, c1, c2;
    Timer timer;

    const int iters = 9000000;

    for (int i=0; i < 3; ++i)
    {
        for (int j=0; j < 3; ++j)
            A1(i,j) = 1.0;
        b1(i) = 1.0;
        b2(i) = 1.0;
    }

    optimizationSink(A1);
    optimizationSink(A2);
    optimizationSink(b1);
    optimizationSink(b2);

    timer.start();
    for (int i=0; i < iters; ++i)
    {
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
        c1 = product(A1,b1);
        c2 = product(A1,b2);
        b1 = product(A1,c1);
        c2 = product(A1,c2);
    }

    timer.stop();

    double ops = 64.0 * iters;
    double flops = ops * 15;
    double seconds = timer.elapsedSeconds();

    double timePerOp = seconds / ops;
    double cycles = timePerOp * (100.0 * 1e+6);

    cout << "ops = " << ops << endl
         << "seconds = " << seconds << endl;

    cout << "timePerOp = " << timePerOp << endl
         << "cycles = " << cycles << endl;

    double Mflops = flops / seconds / 1.0e+6;
    cout << "Mflops = " << Mflops << endl;

    optimizationSink(c1);
    optimizationSink(c2);

    return 0;
}

template<class T>
void optimizationSink(T&)
{
}

