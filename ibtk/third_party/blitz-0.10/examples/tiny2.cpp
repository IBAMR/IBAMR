#define BZ_DISABLE_KCC_COPY_PROPAGATION_KLUDGE

#include <blitz/tinyvec-et.h>

using namespace blitz;

typedef TinyVector<double,3> vec;

inline void reflect(vec& reflection, const vec& incident, 
    const vec& surfaceNormal)
{
    // The surface normal must be unit length to use this equation.

    reflection = incident - 2 * dot(incident,surfaceNormal) * surfaceNormal;
}

template<typename T>
void optimizationSink(T&);

void foo(TinyVector<double,3>& x)
{
    TinyVector<double,3> y, z;

    y =  1.00,  0.40, -1.00;
    z =  0.31,  0.20,  0.93;

    reflect(x, y, z);
}

