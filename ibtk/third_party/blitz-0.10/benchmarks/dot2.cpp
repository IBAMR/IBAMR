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

void init(double* x, int n)
{
    // Completely arbitrary
    for (int i=0; i < n; ++i)
        x[i] = 3.4982938192839824982 * i;
}

const int nmax = 40;

int main()
{
    Timer timer;
    const int iterations1 = 5000000;

    double a1[nmax],a2[nmax],a3[nmax],a4[nmax],a5[nmax],a6[nmax],a7[nmax],a8[nmax],a9[nmax],
        a10[nmax],b1[nmax],b2[nmax],b3[nmax],b4[nmax],b5[nmax],b6[nmax],b7[nmax],b8[nmax],
        b9[nmax],b10[nmax];
    init(a1,nmax);
    init(a2,nmax);
    init(a3,nmax);
    init(a4,nmax);
    init(a5,nmax);
    init(a6,nmax);
    init(a7,nmax);
    init(a8,nmax);
    init(a9,nmax);
    init(a10,nmax);
    init(b1,nmax);
    init(b2,nmax);
    init(b3,nmax);
    init(b4,nmax);
    init(b5,nmax);
    init(b6,nmax);
    init(b7,nmax);
    init(b8,nmax);
    init(b9,nmax);
    init(b10,nmax);

    for (int n=1; n < nmax; ++n)
    {
    int iterations = iterations1 / n;

    timer.start();
    for (int i=0; i < iterations; ++i)
    {
        double result1 = dot(a1,b1,n);
        double result2 = dot(a2,b2,n);
        double result3 = dot(a3,b3,n);
        double result4 = dot(a4,b4,n);
        double result5 = dot(a5,b5,n);
        double result6 = dot(a6,b6,n);
        double result7 = dot(a7,b7,n);
        double result8 = dot(a8,b8,n);
        double result9 = dot(a9,b9,n);
        double result10 = dot(a10,b10,n);
        sink(result1,result2,result3,result4,result5,result6,result7,result8,
            result9,result10);
    }
    timer.stop();

    double Mflops = 10.0 * (n + (n-1)) * iterations / 1e+6;
    cout << n << '\t' << (Mflops/timer.elapsedSeconds()) << endl;
    }

    return 0;
}

