/*****************************************************************************
 * where.cpp     Blitz++ 1D Array example, illustrating where(X,Y,Z)
 *               expressions.
 *****************************************************************************/

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
//  Vector<int> x = Range(-3,+3);   // [ -3 -2 -1  0  1  2  3 ]
    Array<int,1> x(Range(-3,+3));
    x = tensor::i;

    // The where(X,Y,Z) function is similar to the X ? Y : Z operator.
    // If X is logical true, then Y is returned; otherwise, Z is
    // returned.

//  Vector<int> y = where(abs(x) > 2, x+10, x-10);
    Array<int,1> y(x.shape());
    y = where(abs(x) > 2, x+10, x-10);

    // The above statement is transformed into something resembling:
    //
    // for (unsigned i=0; i < 7; ++i)
    //     y[i] = (abs(x[i]) > 2) ? (x[i]+10) : (x[i]-10);
    //
    // The first expression (abs(x) > 2) can involve the usual
    // comparison and logical operators: < > <= >= == != && || 

    cout << x << endl
         << y << endl;

    return 0;
}

