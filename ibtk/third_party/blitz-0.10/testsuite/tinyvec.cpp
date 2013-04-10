#include "testsuite.h"
#include <blitz/array.h>

using namespace blitz;

typedef TinyVector<int,3> tv3;

int main()
{
  tv3 A, B;
    A = 1, 2, 3;
    B = 0, 1, 2;

    BZTEST(A[0]==1);
    BZTEST(A[1]==2);
    BZTEST(A[2]==3);
    BZTEST(B[0]==0);
    BZTEST(B[1]==1);
    BZTEST(B[2]==2);

    A += B*B;
    BZTEST(all(A==tv3(1,3,7)));

    //cout << A << endl;

    A *= 10;
    BZTEST(all(A==tv3(10,30,70)));

    //cout << A << endl;

    tv3::iterator it=B.begin(), end=B.end();
    cout << "(" << *it;
    for (it++; it != end; ++it) 
        cout << "," << *it;
    cout << ")" << endl;

// test tinyvector dot(), product(), and sum() functions

    int dotAB = dot(A,B); 
    BZTEST(dotAB==170);
    int dotAB2 = dot(A,B+B);
    BZTEST(dotAB2==2*dotAB);
    int dotA2B = dot(A+A,B);
    BZTEST(dotA2B==2*dotAB);
    int dotA2B2 = dot(A+A,B+B);
    BZTEST(dotA2B2==4*dotAB);
    int prod1 = product(A);
    BZTEST(prod1==21000);
    int prod2 = product(A+B);
    BZTEST(prod2==22320);
    int sum1 = sum(A);
    BZTEST(sum1==110);
    int sum2 = sum(A-B);
    BZTEST(sum2==sum1-sum(B));

    // (cross product is tested in levicivita.cpp)

    // test funcs
    B=-1,-2,-3;
    A=abs(B);
    BZTEST(A[0]==1);
    BZTEST(A[1]==2);
    BZTEST(A[2]==3);

    TinyVector<double,3> C,D;
    D=1,4,9;
    C=sqrt(D);
    BZTEST(C[0]==1);
    BZTEST(C[1]==2);
    BZTEST(C[2]==3);

    // test expr constructor
    tv3 E(A+B*B);
    BZTEST(all(E==tv3(2,6,12)));

    return 0;
}

