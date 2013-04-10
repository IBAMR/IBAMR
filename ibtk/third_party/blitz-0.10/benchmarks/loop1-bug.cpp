// In KAI C++ 3.2, restrict causes problems for copy propagation.
// Temporary fix: disable restrict

#define BZ_DISABLE_RESTRICT

#include <blitz/vector.h>

BZ_USING_NAMESPACE(blitz)

template<class T> void sink(T&) { }

void foo(int N)
{
    Vector<double> x(N);
    Vector<double> y(N);

    x=sqrt(y);

    sink(x);
    sink(y);
}

