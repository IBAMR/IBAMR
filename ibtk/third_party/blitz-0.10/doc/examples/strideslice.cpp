#include <blitz/array.h>

using namespace blitz;

int main()
{
    Array<int,2> A(8,8);
    A = 0;

    Array<int,2> B = A(Range(1,7,3), Range(1,5,2));
    B = 1;

    cout << "A = " << A << endl;
    return 0;
}

