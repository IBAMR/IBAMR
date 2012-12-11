#include <blitz/array.h>

using namespace blitz;

int main()
{
    Array<int,2> A(6,6), B(3,3);
  
    // Set the upper left quadrant of A to 5 
    A(Range(0,2), Range(0,2)) = 5; 

    // Set the upper right quadrant of A to an identity matrix
    B = 1, 0, 0,
        0, 1, 0,
        0, 0, 1;
    A(Range(0,2), Range(3,5)) = B;

    // Set the fourth row to 1
    A(3, Range::all()) = 1;

    // Set the last two rows to 0
    A(Range(4, toEnd), Range::all()) = 0;

    // Set the bottom right element to 8
    A(5,5) = 8;

    cout << "A = " << A << endl;

    return 0;
}

