#include "testsuite.h"
#include <blitz/array.h>
#include <blitz/array.cc>
#include <blitz/tinyvec2.cc>
#include <blitz/tinymat2.h>
#include <blitz/tinymat2.cc>

BZ_USING_NAMESPACE(blitz)

// Tests that operations on multicomponent arrays work as
// expected. This is a bit tricky because now they involve two
// recursive ET applications.

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

int main()
{
  const int sz=5;

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

  // test that at least a subset was initialized correctly
  ISCLOSE(a1A(2),3);
  ISCLOSE(a1vA(2)(1),-3);
  ISCLOSE(a1mA(2)(1,0),2);
  ISCLOSE(tvvA(2)(1),-3);
  ISCLOSE(a1vvA(2)(1)(0),5);
  ISCLOSE(a1mvA(2)(1,0)(1),-3);

  cout << a1A << endl;
  cout << a1vA << endl;
  cout << a1mA << endl;
  cout << tvvA << endl;
  cout << tmvA << endl;
  cout << a1vvA << endl;
  cout << a1mvA << endl;

  // evaluate a complicated expression to exercise unary, binary,
  // constant, and funcs

  a1B = 2*(-a1A)+sqrt(a1A*a1A);
  a1vB = 2*(-a1vA)+sqrt(a1vA*a1vA);
  a1mB = 2*(-a1mA)+sqrt(a1mA*a1mA);
  tvvB = 2*(-tvvA)+sqrt(tvvA*tvvA);
  tmvB = 2*(-tmvA)+sqrt(tmvA*tmvA);
  a1vvB = 2*(-a1vvA)+sqrt(a1vvA*a1vvA);
  a1mvB = 2*(-a1mvA)+sqrt(a1mvA*a1mvA);

  cout << "\nTesting element-wise expression evaluation:\n";
  cout << a1B << endl;
  cout << a1vB << endl;
  cout << a1mB << endl;
  cout << tvvB << endl;
  cout << tmvB << endl;
  cout << a1vvB << endl;
  cout << a1mvB << endl;

  // test results (we are not testing reductions here so we loop over
  // elements)
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

  // also test that the update versions of the operators work
  a1B -= 2*(-a1A)+sqrt(a1A*a1A);
  a1vB -= 2*(-a1vA)+sqrt(a1vA*a1vA);
  a1mB -= 2*(-a1mA)+sqrt(a1mA*a1mA);
  tvvB -= 2*(-tvvA)+sqrt(tvvA*tvvA);
  tmvB -= 2*(-tmvA)+sqrt(tmvA*tmvA);
  a1vvB -= 2*(-a1vvA)+sqrt(a1vvA*a1vvA);
  a1mvB -= 2*(-a1mvA)+sqrt(a1mvA*a1mvA);

  for(int i=0; i<sz; ++i) {
    ISCLOSE ( a1B(i), 0 );

    for(int j=0; j<2; ++j) {
      ISCLOSE(a1vB(i)(j), 0);
      ISCLOSE(tvvB(i)(j), 0);
      for(int k=0; k<2; ++k) {
	ISCLOSE(a1mB(i)(j,k), 0);
	ISCLOSE(a1vvB(i)(j)(k), 0);
	if(i<2) {
	  ISCLOSE(tmvB(j,k)(i), 0 );
	}
	for(int l=0; l<2; ++l) {
	  ISCLOSE(a1mvB(i)(j,k)(l), 0);
	}
      }
    }
  }

  cout << "\nTesting scalar wrapper:\n";

  // now we also want to test the "scalar" wrapper.
  a1vB = a1vA*scalar(tv(2,-2));
  tmt sc;
  sc=1,2,3,4;
  a1mB = a1mA*scalar(sc);

  cout << a1vB << endl;
  cout << a1mB << endl;

  ISCLOSE(a1vB(1)(1),4);
  ISCLOSE(a1vB(4)(1),-10);
  ISCLOSE(a1mB(1)(0,1),4);
  ISCLOSE(a1mB(2)(1,1),12);

  return 0;
}

