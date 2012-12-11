/*****************************************************************************
 * cast.cpp        Illustration of the cast() function
 *****************************************************************************/

#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Array<int,1> A(4), B(4);
    Array<float,1> C(4);

    A = 1, 2, 3, 5;
    B = 2, 2, 2, 7;

    C = A / B;
    cout << C << endl;

#ifdef BZ_NEW_EXPRESSION_TEMPLATES
    C = A / cast<float>(B);
#else
    C = A / cast(B, float());
#endif
    cout << C << endl;

    return 0;
}

