#include <blitz/array.h>

using namespace blitz;

int main()
{
    Array<float,2> A(3,3), B(3,3), C(3,3);

    A = 1, 0, 0,
        2, 2, 2,
        1, 0, 0;

    B = 0, 0, 7,
        0, 8, 0,
        9, 9, 9;

    C = A + B;

    cout << "A = " << A << endl
         << "B = " << B << endl
         << "C = " << C << endl;

    return 0;
}

