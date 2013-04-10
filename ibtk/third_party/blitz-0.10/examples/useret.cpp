#define BZ_NEW_EXPRESSION_TEMPLATES

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

double myFunction(double x)
{ return 1.0 / (1 + x); }

double foobar(double x, double y)
{
    return x*y;
}

BZ_DECLARE_FUNCTION(myFunction)
BZ_DECLARE_FUNCTION2(foobar)

int main()
{
    Array<double,2> A(4,4), B(4,4), C(4,4);

    A = 0,  1,  2,  3,
        4,  5,  6,  7,
        8,  9, 10, 11,
       12, 13, 14, 15;
    C = 3;

    cout << "A = " << A << endl
         << "C = " << C << endl;

    B = myFunction(A);

    cout << "B = myFunction(A) = " << B << endl;

    B = foobar(A,C);
    cout << "B = foobar(A,C) = " << B << endl;

    B = foobar(tensor::i, tensor::j);
    cout << "B = foobar(tensor::i, tensor::j) = " << B << endl;
}

