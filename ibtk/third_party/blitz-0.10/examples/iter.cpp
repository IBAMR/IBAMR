// Example of using Array<T,N>::iterator 

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    // Create a 4x4 array and fill it with some numbers

    Array<int,2> A(4,4);
    A = tensor::i * 10 + tensor::j; 
    cout << "A = " << A << endl;


    // Use an iterator to list the array elements

    Array<int,2>::iterator iter = A.begin(), end = A.end();

    while (iter != end)
    {
        cout << iter.position() << '\t' << (*iter) << endl;
        ++iter;
    }

    return 0;
}

