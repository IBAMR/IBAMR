#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

using namespace blitz;

template<typename T>
void optimizationSink(T&);

void foo()
{
    TinyMatrix<float,3,3> A;
    TinyVector<float,3> b, c;

    optimizationSink(A);
    optimizationSink(b);

    c = product(A,b);

    optimizationSink(c);
}

