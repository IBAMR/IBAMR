#include <blitz/array.h>

using namespace blitz;

int main()
{
    Array<float,2> A(4,4), B(4,4);
    Array<int,2> C(4,4);

    prettyPrintFormat f1;
    prettyPrintFormat f2(true);

    string s1, s2, s3, s4;

    (A+B).prettyPrint(s1,f1);
    cout << s1 << endl;
    (A+B).prettyPrint(s2,f2);
    cout << s2 << endl;
    (A+B+C).prettyPrint(s3,f1);
    cout << s3 << endl;
    (A+5*C).prettyPrint(s4,f1);
    cout << s4 << endl;
    {
    prettyPrintFormat f2(true);
    string s5;
    (A+5*C+cos(B)*A+A+pow2(B)+C).prettyPrint(s5,f2);
    cout << s5 << endl;
    }

    return 0;
}

