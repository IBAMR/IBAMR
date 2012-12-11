#include <blitz/array.h>
#include <blitz/array/convolve.h>

using namespace blitz;

int main()
{
    Array<float,1> B(Range(-2,+2));
    Array<float,1> C(Range(10,15));

    B = 1, 0, 2, 5, 3;
    C = 10, 2, 4, 1, 7, 2;

    Array<float,1> A = convolve(B,C);

    cout << "A has domain " << A.lbound(0) << "..." << A.ubound(0) << endl
         << A << endl;

    return 0;
}

