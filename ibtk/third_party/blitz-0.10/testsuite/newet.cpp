#define BZ_NEW_EXPRESSION_TEMPLATES

#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

double myFunction(double x)
{
    return 1.0 / (1.0 + x);
}

BZ_DECLARE_FUNCTION(myFunction)

class myFunctor
{
public:
    myFunctor(double p) : period_(p) { }
    double operator()(double x) const
    { return sin(x/period_); }
    BZ_DECLARE_FUNCTOR(myFunctor)
private:
    double period_;
};

class myClass
{
public:
    myClass(double p) : power_(p) { }
    int foo(double x, double y) const
    { return int(x * pow(y,power_)); }
    BZ_DECLARE_MEMBER_FUNCTION2_RET(myClass,foo,int)
private:
    double power_;
};

int max3(int a, int b, int c) { return (a>b?(a>c?a:c):(b>c?b:c)); }
BZ_DECLARE_FUNCTION3(max3)


int main()
{
    Array<int,1> A(4), B(4), C(4), D(4);
    C = 0, 1, 2, 3;
    B = 5, 6, 7, 8;
    A = B + C;
    D = 5, 7, 9, 11;
    BZTEST(count(D==A)==4);

    A = 3 * tensor::i;
    D = 0, 3, 6, 9;
    BZTEST(count(D==A)==4);

    A = 3.0 * tensor::i;
    D = 0, 3, 6, 9;
    BZTEST(count(D==A)==4);

    A = tensor::i * 3.0;
    BZTEST(count(D==A)==4);

    A = 3;
    BZTEST(count(A==3)==4);

    A = tensor::i;
    BZTEST(count(A==C)==4);

    Array<double,1> E(4), F(4);
    E = 2.0 * sqrt(cast<double>(tensor::i * tensor::i));
    F = 2.0 * tensor::i;
    BZTEST(count(E==F)==4);

    E = sqr(tensor::i + 1.0);
    F = pow(tensor::i+1.0,2.0);
    BZTEST(count(E==F)==4);

    E = 4.0 * tensor::i;
    F = myFunction(E);
    E = 1.0 / (1.0 + 4.0 * tensor::i);
    BZTEST(count(E==F)==4);

    myFunctor f1(10.0);
    E = f1(2.5 + tensor::i);
    F = sin((2.5 + tensor::i)/10.0);
    BZTEST(count(E==F)==4);

    myClass f2(3.5);
    A = f2.foo(1.2 - tensor::i, 1.2 + tensor::i);
    B = cast<int>((1.2 - tensor::i) * pow(1.2 + tensor::i, 3.5));
    BZTEST(count(A==B)==4);

    A = 3 * tensor::i;
    B = max3(5,6,7);
    C = 9 - 2 * tensor::i;
    D = 9, 7, 7, 9;
    Array<int,1> DD(4);
    DD = max3(A,B,C);
    BZTEST(count(DD==D)==4);

}

