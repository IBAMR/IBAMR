#include <blitz/array.h>

using namespace blitz;

template<class T>
void sink(T&) { }

void foo()
{
    Vector<float> A(100), B(100);

    A = B + B;
    sink(A);
    sink(B);
}


