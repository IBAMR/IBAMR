/*****************************************************************************
 * slicing.cpp        Blitz++ Array slicing & subarrays example
 *****************************************************************************/

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

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

#ifdef BZ_HAVE_PARTIAL_ORDERING
    A(3, Range::all()) = 1;
#else
    cout << "Warning: your compiler does not support partial ordering of"
         << endl << "member templates; using kludge." << endl;
    A(Range(3,3), Range::all()) = 1;
#endif

    // Set the last two rows to 0
    A(Range(4, toEnd), Range::all()) = 0;

    // Set the bottom right element to 8
    A(5,5) = 8;

    cout << "A = " << A << endl;

    return 0;
}

