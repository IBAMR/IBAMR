#include "testsuite.h"
#include <blitz/array.h>
#include <blitz/array.cc>
#include <blitz/tinyvec2.cc>
#include <blitz/tinymat2.h>
#include <blitz/tinymat2.cc>
// #include <blitz/array/asexpr.cc>
// #include <blitz/array/ops.h>
// #include <blitz/array/ops.cc>
// #include <blitz/array/fastiter.h>
// #include <blitz/array/methods.cc>
// #include <blitz/array/methods.cc>
// #include <blitz/array/io.cc>
// #include <blitz/ops.h>
// #include <blitz/array/funcs.h>

BZ_USING_NAMESPACE(blitz)

// Tests functors on multicomponent arrays.

static const double eps = 0.0001;

#define ISCLOSE(a, b) BZTEST(fabs((a)-(b))<eps)

typedef TinyVector<double,2> tv;
typedef TinyMatrix<double,2,2> tmt;

// types to test
typedef Array<double, 1> a1;
typedef Array<tv, 1> a1v;
typedef Array<tmt, 1> a1m;
typedef TinyVector<tv, 5> tvv;
typedef TinyMatrix<tv, 2,2> tmv;
// and the hardcore multi-level multicomponent types
typedef Array<TinyVector<tv, 2>, 1> a1vv;
typedef Array<TinyMatrix<tv, 2,2>, 1> a1mv;

class ufunc {
public:
  double operator()(double d) const { return 2*(-d)+sqrt(d*d); }
  BZ_DECLARE_FUNCTOR(ufunc);
};

class bifunc {
public:
  double operator()(double a, double b) const { return 2*(-a)+sqrt(b*b); }
  BZ_DECLARE_FUNCTOR2(bifunc);
};

class trifunc {
public:
  double operator()(double a, double b, double c) const
  { return 2*(-a)+sqrt(b*c); }
  BZ_DECLARE_FUNCTOR3(trifunc);
};

class umemfunc {
public:
  double func(double d) const { return 2*(-d)+sqrt(d*d); }
  BZ_DECLARE_MEMBER_FUNCTION(umemfunc, func);
};

class bimemfunc {
public:
  double func(double a, double b) const { return 2*(-a)+sqrt(b*b); }
  BZ_DECLARE_MEMBER_FUNCTION2(bimemfunc, func);
};

class trimemfunc {
public:
  double func(double a, double b, double c) const { return 2*(-a)+sqrt(b*c); }
  BZ_DECLARE_MEMBER_FUNCTION3(trimemfunc, func);
};

const int sz=5;

void test_results(a1& a1A, a1& a1B, a1v& a1vA, a1v& a1vB,
		  a1m& a1mA, a1m& a1mB, tvv& tvvA, tvv& tvvB,
		  tmv& tmvA, tmv& tmvB, a1vv& a1vvA, a1vv& a1vvB,
		  a1mv& a1mvA, a1mv& a1mvB)
{
  for(int i=0; i<sz; ++i) {
    ISCLOSE ( a1B(i), -a1A(i) );

    for(int j=0; j<2; ++j) {
      ISCLOSE(a1vB(i)(j), a1vA(i)(j)>0?-a1vA(i)(j):-3*a1vA(i)(j));
      ISCLOSE(tvvB(i)(j), tvvA(i)(j)>0?-tvvA(i)(j):-3*tvvA(i)(j));
      for(int k=0; k<2; ++k) {
	ISCLOSE(a1mB(i)(j,k),a1mA(i)(j,k)>0?-a1mA(i)(j,k):-3*a1mA(i)(j,k));
	ISCLOSE(a1vvB(i)(j)(k),a1vvA(i)(j)(k)>0?-a1vvA(i)(j)(k):-3*a1vvA(i)(j)(k));
	if(i<2) {
	  ISCLOSE(tmvB(j,k)(i),tmvA(j,k)(i)>0?-tmvA(j,k)(i):-3*tmvA(j,k)(i));
	}
	for(int l=0; l<2; ++l) {
	  ISCLOSE(a1mvB(i)(j,k)(l),a1mvA(i)(j,k)(l)>0?-a1mvA(i)(j,k)(l):-3*a1mvA(i)(j,k)(l));
	}
      }
    }
  }
}


int main()
{

  // create operands
  a1 a1A(sz), a1B(sz), a1C(sz);
  a1v a1vA(sz), a1vB(sz), a1vC(sz);
  a1m a1mA(sz), a1mB(sz), a1mC(sz);
  tvv tvvA, tvvB, tvvC;
  tmv tmvA, tmvB, tmvC;
  a1vv a1vvA(sz), a1vvB(sz), a1vvC(sz);
  a1mv a1mvA(sz), a1mvB(sz), a1mvC(sz);

  // fill them with data
  a1A=1,2,3,4,5;
  a1vA=tv(1,-1),tv(2,-2),tv(3,-3),tv(4),tv(5);
  a1mA=tmt(1),tmt(2),tmt(3),tmt(4),tmt(5);
  a1mA(2)=0,1,2,3;
  tvvA=tv(1,-1),tv(2,-2),tv(3,-3),tv(4),tv(5);
  tmvA=tv(1,-1),tv(2,-2),tv(3,-3),tv(4);
  a1vvA=TinyVector<TinyVector<double,2>,2>(tv(1,-1),tv(2,-2)),
    TinyVector<TinyVector<double,2>,2>(tv(3,-3),tv(3,-3)),
    TinyVector<TinyVector<double,2>,2>(tv(4,-4),tv(5,-5)),
    TinyVector<TinyVector<double,2>,2>(tv(6,6),tv(7,7)),
    TinyVector<TinyVector<double,2>,2>(-8);
  a1mvA=TinyMatrix<TinyVector<double,2>,2,2>(tv(1,-1));
  a1mvA(1)=2;
  a1mvA(2)=tv(1,-1),tv(2,-2),tv(3,-3),tv(4);
  a1mvA(3)(1,1)[1]=-42;

  ufunc uf;
  bifunc bf;
  trifunc tf;
  umemfunc umf;
  bimemfunc bmf;
  trimemfunc tmf;

  // test unary functors
  a1B=uf(a1A);
  a1vB=uf(a1vA);
  a1mB=uf(a1mA);
  tvvB=uf(tvvA);
  tmvB=uf(tmvA);
  a1vvB=uf(a1vvA);
  a1mvB=uf(a1mvA);

  test_results(a1A,a1B,a1vA,a1vB,a1mA,a1mB,tvvA,tvvB,tmvA,
	       tmvB,a1vvA,a1vvB,a1mvA,a1mvB);

  // test binary functors
  a1B=bf(a1A,a1A);
  a1vB=bf(a1vA,a1vA);
  a1mB=bf(a1mA,a1mA);
  tvvB=bf(tvvA,tvvA);
  tmvB=bf(tmvA,tmvA);
  a1vvB=bf(a1vvA,a1vvA);
  a1mvB=bf(a1mvA,a1mvA);

  test_results(a1A,a1B,a1vA,a1vB,a1mA,a1mB,tvvA,tvvB,tmvA,
	       tmvB,a1vvA,a1vvB,a1mvA,a1mvB);

  // test ternary functors
  a1B=tf(a1A,a1A,a1A);
  a1vB=tf(a1vA,a1vA,a1vA);
  a1mB=tf(a1mA,a1mA,a1mA);
  tvvB=tf(tvvA,tvvA,tvvA);
  tmvB=tf(tmvA,tmvA,tmvA);
  a1vvB=tf(a1vvA,a1vvA,a1vvA);
  a1mvB=tf(a1mvA,a1mvA,a1mvA);

  test_results(a1A,a1B,a1vA,a1vB,a1mA,a1mB,tvvA,tvvB,tmvA,
	       tmvB,a1vvA,a1vvB,a1mvA,a1mvB);

  // test unary member functions
  a1B=umf.func(a1A);
  a1vB=umf.func(a1vA);
  a1mB=umf.func(a1mA);
  tvvB=umf.func(tvvA);
  tmvB=umf.func(tmvA);
  a1vvB=umf.func(a1vvA);
  a1mvB=umf.func(a1mvA);

  test_results(a1A,a1B,a1vA,a1vB,a1mA,a1mB,tvvA,tvvB,tmvA,
	       tmvB,a1vvA,a1vvB,a1mvA,a1mvB);

  // test binary member functions
  a1B=bmf.func(a1A,a1A);
  a1vB=bmf.func(a1vA,a1vA);
  a1mB=bmf.func(a1mA,a1mA);
  tvvB=bmf.func(tvvA,tvvA);
  tmvB=bmf.func(tmvA,tmvA);
  a1vvB=bmf.func(a1vvA,a1vvA);
  a1mvB=bmf.func(a1mvA,a1mvA);

  test_results(a1A,a1B,a1vA,a1vB,a1mA,a1mB,tvvA,tvvB,tmvA,
	       tmvB,a1vvA,a1vvB,a1mvA,a1mvB);

  // test ternary member functions
  a1B=tmf.func(a1A,a1A,a1A);
  a1vB=tmf.func(a1vA,a1vA,a1vA);
  a1mB=tmf.func(a1mA,a1mA,a1mA);
  tvvB=tmf.func(tvvA,tvvA,tvvA);
  tmvB=tmf.func(tmvA,tmvA,tmvA);
  a1vvB=tmf.func(a1vvA,a1vvA,a1vvA);
  a1mvB=tmf.func(a1mvA,a1mvA,a1mvA);

  test_results(a1A,a1B,a1vA,a1vB,a1mA,a1mB,tvvA,tvvB,tmvA,
	       tmvB,a1vvA,a1vvB,a1mvA,a1mvB);

  // make sure functors can be applied to expressions too (the
  // functors and member functions use the same machinery so if one
  // works, the others should, too...)
  a1B=uf(1.0*a1A);
  a1vB=uf(1.0*a1vA);
  a1mB=uf(1.0*a1mA);
  tvvB=uf(1.0*tvvA);
  tmvB=uf(1.0*tmvA);
  a1vvB=uf(1.0*a1vvA);
  a1mvB=uf(1.0*a1mvA);

  test_results(a1A,a1B,a1vA,a1vB,a1mA,a1mB,tvvA,tvvB,tmvA,
	       tmvB,a1vvA,a1vvB,a1mvA,a1mvB);

  a1B=bmf.func(1.0*a1A,a1A+0);
  a1vB=bmf.func(1.0*a1vA,a1vA+0);
  a1mB=bmf.func(1.0*a1mA,a1mA+0);
  tvvB=bmf.func(1.0*tvvA,tvvA+0);
  tmvB=bmf.func(1.0*tmvA,tmvA+0);
  a1vvB=bmf.func(1.0*a1vvA,a1vvA+0);
  a1mvB=bmf.func(1.0*a1mvA,a1mvA+0);

  test_results(a1A,a1B,a1vA,a1vB,a1mA,a1mB,tvvA,tvvB,tmvA,
	       tmvB,a1vvA,a1vvB,a1mvA,a1mvB);

  a1B=tmf.func(1.0*a1A,a1A+0,a1A-0);
  a1vB=tmf.func(1.0*a1vA,a1vA+0,a1vA-0);
  a1mB=tmf.func(1.0*a1mA,a1mA+0,a1mA-0);
  tvvB=tmf.func(1.0*tvvA,tvvA+0,tvvA-0);
  tmvB=tmf.func(1.0*tmvA,tmvA+0,tmvA-0);
  a1vvB=tmf.func(1.0*a1vvA,a1vvA+0,a1vvA-0);
  a1mvB=tmf.func(1.0*a1mvA,a1mvA+0,a1mvA-0);

  test_results(a1A,a1B,a1vA,a1vB,a1mA,a1mB,tvvA,tvvB,tmvA,
	       tmvB,a1vvA,a1vvB,a1mvA,a1mvB);

  return 0;
}

