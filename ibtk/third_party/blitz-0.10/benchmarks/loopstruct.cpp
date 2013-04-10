// CC -64 -LANG:std -LANG:restrict -Ofast -PHASE:clist -IPA=off -IPA:INLINE=off
//
// When compiled with the above options, these are the results on convex:
// func1: 34.6484
// func2: 24.6603
// func3: 17.2822
//
// func1 is a simple C-style loop.
// func2 has the pointers stuck inside a struct.  Prefetching no longer occurs.
// func3 has the pointers inside "iterator" structs, and a read is done
//       using an inline operator().  This somehow results in the loop
//       being unrolled only twice, instead of 4 times.

#include <iostream.h>
#include <sys/resource.h>

// Paste the Timer class in here so you don't have to have Blitz++

class Timer {

public:
    Timer() 
    { 
        state_ = uninitialized;
    }

    void start()
    { 
        state_ = running;
        t1_ = systemTime();
    }

    void stop()
    {
        t2_ = systemTime();
        state_ = stopped;
    }

/* Compaq cxx compiler in ansi mode cannot print out long double type! */
#if defined(__DECCXX)
    double elapsedSeconds()
#else
    long double elapsedSeconds()
#endif
    {
        return t2_ - t1_;
    }

private:
    Timer(Timer&) { }
    void operator=(Timer&) { }

    long double systemTime()
    {
        getrusage(RUSAGE_SELF, &resourceUsage_);
        double seconds = resourceUsage_.ru_utime.tv_sec 
            + resourceUsage_.ru_stime.tv_sec;
        double micros  = resourceUsage_.ru_utime.tv_usec 
            + resourceUsage_.ru_stime.tv_usec;
        return seconds + micros/1.0e6;
    }

    enum { uninitialized, running, stopped } state_;

    struct rusage resourceUsage_;
    long double t1_, t2_;
};


// func1:  Simple version

void func1(double* restrict x, double* restrict y, double a, int N)
{
    for (int i=0; i < N; ++i)
        y[i] += a*x[i];
}


// func2: With pointers inside a struct

struct A {
    double* restrict x;
    double* restrict y;
    double a;
    int N;
};

void func2(A& z)
{
    for (int i=0; i < z.N; ++i)
        z.y[i] += z.a * z.x[i];
}


// func3: with very simple "iterators" (the B struct).

struct B {
    int q;
    double* restrict data;
    double operator()(int i)
    { return data[i]; }
};

struct C {
    B x;
    B y;
    double a;
    int N;
};

void func3(C& z)
{
    for (int i=0; i < z.N; ++i)
        z.y.data[i] += z.a * z.x(i);
}


// Initialize array

void init(double* x, int N)
{
    for (int i=0; i < N; ++i)
        x[i] = 1.0;
}

int main()
{
    Timer timer;

    int N = 1000000;
    int iters = 20;
    double Mflops = N * iters * 2 / 1000000.;

    double* x = new double[N];
    double* y = new double[N];
    double a = .14989182;
    init(x,N);
    init(y,N);

    timer.start();
    for (int i=0; i < iters; ++i)
        func1(x,y,a,N);
    timer.stop();

    cout << "func1: " << Mflops/timer.elapsedSeconds() << endl;

    timer.start();
    A z;
    z.x = x;
    z.y = y;
    z.a = a;
    z.N = N;
    for (int i=0; i < iters; ++i)
        func2(z);
    timer.stop();

    cout << "func2: " << Mflops/timer.elapsedSeconds() << endl;

    timer.start();
    {
    C z;
    z.x.data = x;
    z.y.data = y;
    z.a = a;
    z.N = N;
    for (int i=0; i < iters; ++i)
        func3(z);
    }
    timer.stop();
    cout << "func3: " << Mflops/timer.elapsedSeconds() << endl;
}

