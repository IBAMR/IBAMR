#include <blitz/timer.h>

BZ_USING_NAMESPACE(blitz)

void initialize(double& c, double& d, double* a, double* b, int& N);

template<class T>
void sink(T&)
{ }

void benchmarkLoops(int, long);

int main()
{
    cout << "This program measures the performance of DAXPY operations" 
         << endl << "using various C loop structures." << endl << endl;

    cout << endl << "In-cache:" << endl;

    benchmarkLoops(400,50000);

    cout << endl << "Out of cache:" << endl;

    benchmarkLoops(1000000,50);

    return 0;
}

void benchmarkLoops(int N, long iterations)
{
    double* _bz_restrict a = new double[N];
    double* _bz_restrict b = new double[N];
    double c, d;
    double t1, t2;

    initialize(c, d, a, b, N);

    double mflops = iterations * 4.0 * N / (1024.0 * 1024.0);

    Timer timer;

    cout << "Mflops/s Description" << endl;

    long iter;
    int i;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        for (i=0; i < N; ++i)
            a[i] += c * b[i];

        for (i=0; i < N; ++i)
            a[i] += d * b[i];
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds()) 
         << "   for, indirection, unit stride" << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        for (i=0; i < N; ++i)
            a[i] = a[i] + c * b[i];

        for (i=0; i < N; ++i)
            a[i] = a[i] + d * b[i];
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "   for, indirection, unit stride, no +=" << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        for (i=N-1; i >= 0; --i)
            a[i] += c * b[i];

        for (i=N-1; i >= 0; --i)
            a[i] += d * b[i];
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "   for, indirection, unit stride, backwards loops" << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        double c2 = c;

        int n1 = N & 3;
        for (i=0; i < n1; ++i)
            a[i] += c2 * b[i];

        for (; i < N; i += 4)
        {
            a[i] += c2 * b[i];
            a[i+1] += c2 * b[i+1];
            a[i+2] += c2 * b[i+2];
            a[i+3] += c2 * b[i+3];
        }

        double d2 = d;
        int n2 = N & 3;
        for (i=0; i < n2; ++i)
            a[i] += d2 * b[i];

        for (; i < N; i += 4)
        {
            a[i] += d2 * b[i];
            a[i+1] += d2 * b[i+1];
            a[i+2] += d2 * b[i+2];
            a[i+3] += d2 * b[i+3];
        } 
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    for, unroll=4, unit stride, constants loaded into temps"
         << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        double c2 = c;

        int n1 = N & 3;
        for (i=0; i < n1; ++i)
            a[i] += c2 * b[i];

        for (; i < N; i += 4)
        {
            double t1 = c2 * b[i];
            double t2 = c2 * b[i+1];
            double t3 = c2 * b[i+2];
            double t4 = c2 * b[i+3];

            a[i] += t1;
            a[i+1] += t2;
            a[i+2] += t3;
            a[i+3] += t4;
        }

        double d2 = d;
        int n2 = N & 3;
        for (i=0; i < n2; ++i)
            a[i] += d2 * b[i];

        for (; i < N; i += 4)
        {
            double t1 = d2 * b[i];
            double t2 = d2 * b[i+1];
            double t3 = d2 * b[i+2];
            double t4 = d2 * b[i+3];

            a[i] += t1;
            a[i+1] += t2;
            a[i+2] += t3;
            a[i+3] += t4;
        }
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    for, unroll=4, unit stride, constants loaded into temps,"
         << endl << "\t\t4 read then 4 write" 
         << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        double c2 = c;

        int n1 = N & 3;
        for (i=0; i < n1; ++i)
            a[i] += c2 * b[i];

        for (; i < N; i += 4)
        {
            a[i] = a[i] + c2 * b[i];
            a[i+1] = a[i+1] + c2 * b[i+1];
            a[i+2] = a[i+2] + c2 * b[i+2];
            a[i+3] = a[i+3] + c2 * b[i+3];
        }

        double d2 = d;
        int n2 = N & 3;
        for (i=0; i < n2; ++i)
            a[i] += d2 * b[i];

        for (; i < N; i += 4)
        {
            a[i] = a[i] + d2 * b[i];
            a[i+1] = a[i+1] + d2 * b[i+1];
            a[i+2] = a[i+2] + d2 * b[i+2];
            a[i+3] = a[i+3] + d2 * b[i+3];
        }
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    for, unroll=4, unit stride, constants loaded into temps,"
         << endl << "            no += "
         << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        double c2 = c;

        int n1 = N & 3;
        for (i=0; i < n1; ++i)
            a[i] += c2 * b[i];

        for (; i < N; i += 4)
        {
            int i1 = i + 1;
            a[i] += c2 * b[i];
            int i2 = i + 2;
            a[i1] += c2 * b[i1];
            int i3 = i + 3;
            a[i2] += c2 * b[i2];
            a[i3] += c2 * b[i3];
        }

        double d2 = d;
        int n2 = N & 3;
        for (i=0; i < n2; ++i)
            a[i] += d2 * b[i];

        for (; i < N; i += 4)
        {
            int i1 = i + 1;
            a[i] += d2 * b[i];
            int i2 = i + 2;
            a[i1] += d2 * b[i1];
            int i3 = i + 3;
            a[i2] += d2 * b[i2];
            a[i3] += d2 * b[i3];
        }
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    for, unroll=4, unit stride, constants loaded into temps,"
         << endl << "        CSE for index offsets"
         << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        double c2 = c;

        int n1 = N & 3;
        for (i=0; i < n1; ++i)
            a[i] += c2 * b[i];

        double* pa = a+n1;
        double* pb = b+n1;
 
        int top = N - n1 - 4;

        for (i=top; i >= 0; i -= 4)
        {
            pa[i] += c2 * pb[i];
            pa[i+1] += c2 * pb[i+1];
            pa[i+2] += c2 * pb[i+2];
            pa[i+3] += c2 * pb[i+3];
        }

        double d2 = d;
        int n2 = N & 3;
        for (i=0; i < n2; ++i)
            a[i] += d2 * b[i];

        pa = a+n2;
        pb = b+n2;

        top = N - n2 - 4;
        for (i=top; i >= 0; i -= 4)
        {
            pa[i] += d2 * pb[i];
            pa[i+1] += d2 * pb[i+1];
            pa[i+2] += d2 * pb[i+2];
            pa[i+3] += d2 * pb[i+3];
        }
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    for, unroll=4, unit stride, constants loaded into temps,"
         << "            backwards"
         << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        double c2 = c;

        int n1 = N & 7;
        for (i=0; i < n1; ++i)
            a[i] += c2 * b[i];

        for (; i < N; i += 8)
        {
            a[i] += c2 * b[i];
            a[i+1] += c2 * b[i+1];
            a[i+2] += c2 * b[i+2];
            a[i+3] += c2 * b[i+3];
            a[i+4] += c2 * b[i+4];
            a[i+5] += c2 * b[i+5];
            a[i+6] += c2 * b[i+6];
            a[i+7] += c2 * b[i+7];
        }

        double d2 = d;
        int n2 = N & 7;
        for (i=0; i < n2; ++i)
            a[i] += d2 * b[i];

        for (; i < N; i += 8)
        {
            a[i] += d2 * b[i];
            a[i+1] += d2 * b[i+1];
            a[i+2] += d2 * b[i+2];
            a[i+3] += d2 * b[i+3];
            a[i+4] += d2 * b[i+4];
            a[i+5] += d2 * b[i+5];
            a[i+6] += d2 * b[i+6];
            a[i+7] += d2 * b[i+7];
        }
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    for, unroll=8, unit stride, constants loaded into temps"
         << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        double c2 = c;
        for (i=0; i < N; ++i)
            a[i] += c2 * b[i];

        double d2 = d;
        for (i=0; i < N; ++i)
            a[i] += d2 * b[i];
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "   for, indirection, unit stride, constants into temps"
         << endl;

    /*********************************************************************/
    
    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        int stride = 1;
        sink(stride);    // Prevent copy propagation

        for (i=0; i < N; i += stride)
            a[i] += c * b[i];

        for (i=0; i < N; i += stride)
            a[i] += d * b[i];
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    for, indirection, non-unit stride" << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        int stride = 1;
        sink(stride);    // Prevent copy propagation

        double c2 = c;
        for (i=0; i < N; i += stride)
            a[i] += c2 * b[i];

        double d2 = d;
        for (i=0; i < N; i += stride)
            a[i] += d2 * b[i];
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    for, indirection, non-unit stride, constants "
            "loaded into temps" << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        double * _bz_restrict pa1 = a,
               * _bz_restrict pb1 = b;
        double * _bz_restrict paend1 = a + N;
        while (pa1 != paend1)
        {
            *pa1 += c * (*pb1);
            ++pa1;
            ++pb1;
        }

        double * _bz_restrict pa2 = a,
               * _bz_restrict pb2 = b;
        double * _bz_restrict paend2 = a + N;
        while (pa2 != paend2)
        {
            *pa2 += d * (*pb2);
            ++pa2;
            ++pb2;
        }
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    while, pointer increment, unit stride" << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        double * _bz_restrict pa1 = a,
               * _bz_restrict pb1 = b;
        double * _bz_restrict paend1 = a + N;
        double c2 = c;
        while (pa1 != paend1)
        {
            *pa1 += c2 * (*pb1);
            ++pa1;
            ++pb1;
        }

        double * _bz_restrict pa2 = a,
               * _bz_restrict pb2 = b;
        double * _bz_restrict paend2 = a + N;
        double d2 = d;
        while (pa2 != paend2)
        {
            *pa2 += d2 * (*pb2);
            ++pa2;
            ++pb2;
        }
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    while, pointer increment, unit stride, " << endl
         << "    constants loaded into temps" 
         << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        int stride = 1;
        sink(stride);

        double * _bz_restrict pa1 = a,
               * _bz_restrict pb1 = b;
        double * _bz_restrict paend1 = a + N * stride;
        while (pa1 != paend1)
        {
            *pa1 += c * (*pb1);
            pa1 += stride;
            pb1 += stride;
        }

        double * _bz_restrict pa2 = a,
               * _bz_restrict pb2 = b;
        double * _bz_restrict paend2 = a + N * stride;
        while (pa2 != paend2)
        {
            *pa2 += d * (*pb2);
            pa2 += stride;
            pb2 += stride;
        }
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    while, pointer increment, non-unit stride" << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        int stride = 1;
        sink(stride);

        double * _bz_restrict pa1 = a,
               * _bz_restrict pb1 = b;
        double * _bz_restrict paend1 = a + N * stride;
        double c2 = c;
        int n1 = N & 3;

        for (i=0; i < n1; ++i)
        {
            *pa1 += c2 * (*pb1);
            pa1 += stride;
            pb1 += stride;
        }

        while (pa1 != paend1)
        {
            pa1[0] += c2 * pb1[0];
            pa1[1] += c2 * pb1[1];
            pa1[2] += c2 * pb1[2];
            pa1[3] += c2 * pb1[3];
            pa1 += 4 * stride;
            pb1 += 4 * stride;
        }

        double * _bz_restrict pa2 = a,
               * _bz_restrict pb2 = b;
        double * _bz_restrict paend2 = a + N * stride;
        double d2 = d;
        int n2 = N & 3;

        for (i=0; i < n2; ++i)
        {
            *pa2 += d2 * (*pb2);
            pa2 += stride;
            pb2 += stride;
        }

        while (pa2 != paend2)
        {
            pa2[0] += d2 * pb2[0];
            pa2[1] += d2 * pb2[1];
            pa2[2] += d2 * pb2[2];
            pa2[3] += d2 * pb2[3];
            pa2 += 4 * stride;
            pb2 += 4 * stride;
        }
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    while, pointer increment, unroll=4, non-unit stride," << endl
         << "     constants loaded into temps" << endl;

    /*********************************************************************/

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        double c2 = c;

        int n1 = N & 3;
        for (i=0; i < n1; ++i)
            a[i] += c2 * b[i];

        for (; i < N; i += 4)
        {
            t1 = a[i+4];
            a[i] += c2 * b[i];
            a[i+1] += c2 * b[i+1];
            t2 = b[i+4];
            a[i+2] += c2 * b[i+2];
            a[i+3] += c2 * b[i+3];
        }

        double d2 = d;
        int n2 = N & 3;
        for (i=0; i < n2; ++i)
            a[i] += d2 * b[i];

        for (; i < N; i += 4)
        {
            t1 = a[i+4];
            a[i] += d2 * b[i];
            a[i+1] += d2 * b[i+1];
            t2 = b[i+4];
            a[i+2] += d2 * b[i+2];
            a[i+3] += d2 * b[i+3];
        }
    }
    timer.stop();

    
    sink(t1);
    sink(t2);

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    for, unroll=4, unit stride, constants loaded into temps,"
         << "            prefetching"
         << endl;

    /********************************************************************/

    struct vectorPair {
        double a;
        double b;
    };
    vectorPair* v = new vectorPair[N];
    int N2 = 2*N;
    initialize(c, d, (double*)v, (double*)v, N2);

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        for (i=0; i < N; ++i)
            v[i].a += c * v[i].b;

        for (i=0; i < N; ++i)
            v[i].a += d * v[i].b;
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "   interlaced, for, indirection, unit stride" << endl;

    /*********************************************************************/

    initialize(c, d, (double*)v, (double*)v, N2);

    timer.start();
    for (iter=0; iter < iterations; ++iter)
    {
        double c2 = c;

        int n1 = N & 3;
        for (i=0; i < n1; ++i)
            v[i].a += c2 * v[i].b;

        for (; i < N; i += 4)
        {
            v[i].a += c2 * v[i].b;
            v[i+1].a += c2 * v[i+1].b;
            v[i+2].a += c2 * v[i+2].b;
            v[i+3].a += c2 * v[i+3].b;
        }

        double d2 = d;
        int n2 = N & 3;
        for (i=0; i < n2; ++i)
            v[i].a += d2 * v[i].b;

        for (; i < N; i += 4)
        {
            v[i].a += d2 * v[i].b;
            v[i+1].a += d2 * v[i+1].b;
            v[i+2].a += d2 * v[i+2].b;
            v[i+3].a += d2 * v[i+3].b;
        }
    }
    timer.stop();

    cout << setw(7) << setprecision(5) << (mflops/timer.elapsedSeconds())
         << "    for, unroll=4, unit stride, interlaced, " << endl
         << "\t\tconstants loaded into temps"
         << endl;

    delete [] v;

    /********************************************************************/

    delete [] a;
    delete [] b;
}

void initialize(double& c, double& d, double* a, double* b, int& N)
{
    for (int i=0; i < N; ++i)
    {
        a[i] = 1/7.; 
        b[i] = 1/3.;
    }
    c = 0.398192839842;
    d = - c;
}

