#include <blitz/array.h>

using namespace blitz;

int main()
{
    Array<int,2> A(4,5,FortranArray<2>());
    firstIndex i;
    secondIndex j;
    A = 10*i + j;

    cout << "A = " << A << endl;

    Array<float,1> B(20);
    B = exp(-i/100.);
    
    cout << "B = " << endl << B << endl;

    return 0;
}

