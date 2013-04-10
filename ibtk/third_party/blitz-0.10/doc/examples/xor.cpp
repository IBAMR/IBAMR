#include <blitz/array.h>

using namespace blitz;

int main()
{
    Array<int,1> A(4), B(4), C(4);
    A = 0, 1, 2, 3;
    B = 0, 1, 0, 3;
    C = A ^ B;
    cout << C << endl;
}

