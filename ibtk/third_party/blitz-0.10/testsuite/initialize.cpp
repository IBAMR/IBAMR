#include "testsuite.h"
#include <blitz/vector2.h>
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main()
{
    Vector<int> a(6);
    a = 1, 4, 6, 3, 7, 8;

    BZTEST(a(0) == 1);
    BZTEST(a(1) == 4);
    BZTEST(a(2) == 6);
    BZTEST(a(3) == 3);
    BZTEST(a(4) == 7);
    BZTEST(a(5) == 8);
    cout << a << endl;

    Vector<int> b = a(Range(0,4,2));
    cout << b << endl;
//  b = -1, -2, -3;
//  workaround: use list initialization with Vector that has own storage
    Vector<int> bb(3);
    bb = -1, -2, -3;
    b = bb;
    cout << a << endl;
    cout << b << endl;
    BZTEST(a(0) == -1);
    BZTEST(a(1) == 4);
    BZTEST(a(2) == -2);
    BZTEST(a(3) == 3);
    BZTEST(a(4) == -3);
    BZTEST(a(5) == 8);

    Vector<int> c = a.reverse(firstDim);
    c = 0;
    BZTEST(a(0) == 0);
    BZTEST(a(1) == 0);
    BZTEST(a(2) == 0);
    BZTEST(a(3) == 0);
    BZTEST(a(4) == 0);
    BZTEST(a(5) == 0);
//  c = 8, 7, 3, 6, 4, 1;
//  workaround: use list initialization with Vector that has own storage
    Vector<int> cc(6);
    cc = 8, 7, 3, 6, 4, 1;
    c = cc;
    BZTEST(a(0) == 1);
    BZTEST(a(1) == 4);
    BZTEST(a(2) == 6);
    BZTEST(a(3) == 3);
    BZTEST(a(4) == 7);
    BZTEST(a(5) == 8);

    Array<int,1> d(5);
    d=1,-1,2,-2,0;
    BZTEST(d(0) == 1);
    BZTEST(d(1) == -1);
    BZTEST(d(2) == 2);
    BZTEST(d(3) == -2);
    BZTEST(d(4) == 0);

    
    return 0;
}

