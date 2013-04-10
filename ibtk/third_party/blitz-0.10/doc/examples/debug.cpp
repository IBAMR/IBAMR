#include <blitz/array.h>

using namespace blitz;

int main()
{
    Array<complex<float>, 2> Z(4,4);

    Z = complex<float>(0.0, 1.0);

    Z(4,4) = complex<float>(1.0, 0.0);

    return 0;
}

