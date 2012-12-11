#include <blitz/array.h>

using namespace blitz;

int main()
{
    Array<float,1> x(4), y(4);
    Array<float,2> A(4,4);

    x = 1, 2, 3, 4;
    y = 1, 0, 0, 1;

    firstIndex i;
    secondIndex j;

    A = x(i) * y(j);

    cout << A << endl;

    return 0;
}

