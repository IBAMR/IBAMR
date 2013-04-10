#include <blitz/array.h>

using namespace blitz;

int main()
{
    Array<int,1> A(7);
    A = 0, 1, 2, 3, 4, 5, 6;

    cout << A(Range::all())  << endl          // [ 0 1 2 3 4 5 6 ]
         << A(Range(3,5))    << endl          // [ 3 4 5 ]
         << A(Range(3,toEnd)) << endl         // [ 3 4 5 6 ]
         << A(Range(fromStart,3)) << endl     // [ 0 1 2 3 ]
         << A(Range(1,5,2)) << endl           // [ 1 3 5 ]
         << A(Range(5,1,-2)) << endl          // [ 5 3 1 ]
         << A(Range(fromStart,toEnd,2)) << endl;    // [ 0 2 4 6 ]

    return 0;
}
