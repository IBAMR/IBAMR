#include "testsuite.h"
#include "blitz/array.h"
using namespace blitz;

// test functors

class crap {
public:
  double a(const double x) const { return 2*x; }
  double b(const double x, const double y) const { return y*x; }
  double operator()(const double x) const {return a(x); }
  BZ_DECLARE_FUNCTOR_RET(crap, double);
  BZ_DECLARE_MEMBER_FUNCTION_RET(crap, a, double);
  BZ_DECLARE_MEMBER_FUNCTION2_RET(crap, b, double);
};

class crap2 {
public:
  double operator()(const double x, const double y) const { return y*x; }
  BZ_DECLARE_FUNCTOR2_RET(crap2, double);
};

class crap3 {
public:
  double operator()(const double x, const double y, const double z) const
  { return z*y*x; }
  BZ_DECLARE_FUNCTOR3_RET(crap3, double);
};

int main(int, char**)
{
    crap c;
    crap2 c2;
    crap3 c3;

    Array<double,1> lam(10);
    double l = 0.0;
    for (Array<double,1>::iterator i=lam.begin();i!=lam.end();++i)
        *i = l++;

    Array<double,1> a(c.a(lam));
    for (int i=0;i<10;++i)
        BZTEST(a(i)==2.0*i);

    a=c.b(lam,lam);
    for (int i=0;i<10;++i)
        BZTEST(a(i)==i*i);

    a=c(lam);
    for (int i=0;i<10;++i)
        BZTEST(a(i)==2*i);

    a=c2(lam,lam);
    for (int i=0;i<10;++i)
        BZTEST(a(i)==i*i);

    a=c3(lam,lam,lam);
    for (int i=0;i<10;++i)
        BZTEST(a(i)==i*i*i);

    return 0;
}
