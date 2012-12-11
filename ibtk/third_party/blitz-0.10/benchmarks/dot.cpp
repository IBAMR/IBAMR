#include <blitz/tinyvec.h>
#include <blitz/timer.h>

BZ_USING_NAMESPACE(blitz)

double dot(const double* a, const double* b, int n)
{
    double result = 0.;
    for (int i=0; i < n; ++i)
        result += a[i] * b[i];

    return result;
}

template<class T>
void sink(T&)
{
}

void sink(double,double,double,double,double,double,double,double,double,double)
{
}

int main()
{
    Timer timer;
    const int iterations = 1000000;

    double a1[3] = { 0.3242343, 0.1429833, 0.43988583 };
    double b1[3] = { 0.3989421, 0.9854983, 0.58439328 };
    double a2[3] = { 0.3242343, 0.1429833, 0.43988583 };
    double b2[3] = { 0.3989421, 0.9854983, 0.58439328 };
    double a3[3] = { 0.3242343, 0.1429833, 0.43988583 };
    double b3[3] = { 0.3989421, 0.9854983, 0.58439328 };
    double a4[3] = { 0.3242343, 0.1429833, 0.43988583 };
    double b4[3] = { 0.3989421, 0.9854983, 0.58439328 };
    double a5[3] = { 0.3242343, 0.1429833, 0.43988583 };
    double b5[3] = { 0.3989421, 0.9854983, 0.58439328 };
    double a6[3] = { 0.3242343, 0.1429833, 0.43988583 };
    double b6[3] = { 0.3989421, 0.9854983, 0.58439328 };
    double a7[3] = { 0.3242343, 0.1429833, 0.43988583 };
    double b7[3] = { 0.3989421, 0.9854983, 0.58439328 };
    double a8[3] = { 0.3242343, 0.1429833, 0.43988583 };
    double b8[3] = { 0.3989421, 0.9854983, 0.58439328 };
    double a9[3] = { 0.3242343, 0.1429833, 0.43988583 };
    double b9[3] = { 0.3989421, 0.9854983, 0.58439328 };
    double a10[3] = { 0.3242343, 0.1429833, 0.43988583 };
    double b10[3] = { 0.3989421, 0.9854983, 0.58439328 };

    timer.start();
    for (int i=0; i < iterations; ++i)
    {
        double result1 = dot(a1,b1,3);
        double result2 = dot(a2,b2,3);
        double result3 = dot(a3,b3,3);
        double result4 = dot(a4,b4,3);
        double result5 = dot(a5,b5,3);
        double result6 = dot(a6,b6,3);
        double result7 = dot(a7,b7,3);
        double result8 = dot(a8,b8,3);
        double result9 = dot(a9,b9,3);
        double result10 = dot(a10,b10,3);
        sink(result1,result2,result3,result4,result5,result6,result7,result8,
            result9,result10);
    }
    timer.stop();

    double Mflops = 10.0 * 5.0 * iterations / 1e+6;
    cout << "Nonspecialized algorithm: " << (Mflops/timer.elapsedSeconds())
         << endl;

    TinyVector<double,3> c1, c2, c3, c4, c5, c6, c7, c8, c9, c10,
        d1, d2, d3, d4, d5, d6, d7, d8, d9, d10;
    c1 = 0.3242343, 0.1429833, 0.43988583;
    d1 = 0.3989421, 0.9854983, 0.58439328;
    sink(c1);
    sink(d1);
    c2 = 0.3242343, 0.1429833, 0.43988583;
    d2 = 0.3989421, 0.9854983, 0.58439328;
    sink(c2);
    sink(d2);
    c3 = 0.3242343, 0.1429833, 0.43988583;
    d3 = 0.3989421, 0.9854983, 0.58439328;
    sink(c3);
    sink(d3);
    c4 = 0.3242343, 0.1429833, 0.43988583;
    d4 = 0.3989421, 0.9854983, 0.58439328;
    sink(c4);
    sink(d4);
    c5 = 0.3242343, 0.1429833, 0.43988583;
    d5 = 0.3989421, 0.9854983, 0.58439328;
    sink(c5);
    sink(d5);
    c6 = 0.3242343, 0.1429833, 0.43988583;
    d6 = 0.3989421, 0.9854983, 0.58439328;
    sink(c6);
    sink(d6);
    c7 = 0.3242343, 0.1429833, 0.43988583;
    d7 = 0.3989421, 0.9854983, 0.58439328;
    sink(c7);
    sink(d7);
    c8 = 0.3242343, 0.1429833, 0.43988583;
    d8 = 0.3989421, 0.9854983, 0.58439328;
    sink(c8);
    sink(d8);
    c9 = 0.3242343, 0.1429833, 0.43988583;
    d9 = 0.3989421, 0.9854983, 0.58439328;
    sink(c9);
    sink(d9);
    c10 = 0.3242343, 0.1429833, 0.43988583;
    d10 = 0.3989421, 0.9854983, 0.58439328;
    sink(c10);
    sink(d10);

    timer.start();
    for (int i=0; i < iterations; ++i)
    {
        double result1 = dot(c1, d1);
        double result2 = dot(c2, d2);
        double result3 = dot(c3, d3);
        double result4 = dot(c4, d4);
        double result5 = dot(c5, d5);
        double result6 = dot(c6, d6);
        double result7 = dot(c7, d7);
        double result8 = dot(c8, d8);
        double result9 = dot(c9, d9);
        double result10 = dot(c10, d10);
        sink(result1, result2, result3, result4, result5, result6, result7,
            result8, result9, result10);
    }
    timer.stop();
    cout << "Metaprogram: " << (Mflops/timer.elapsedSeconds())
         << endl;

    return 0;
}

