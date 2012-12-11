/*****************************************************************************
 * simple.cpp        Some simple array operations
 *****************************************************************************/

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<double,1> x(100);
    x = tensor::i;          // x = [ 0, 1, 2, ..., 99 ]

    Array<double,1> z(x + 150);
    Array<double,1> v(z + x * 2);

    cout << v << endl;
}

