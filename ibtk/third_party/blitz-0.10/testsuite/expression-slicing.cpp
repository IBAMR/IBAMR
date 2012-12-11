#include "testsuite.h"
#include <blitz/array.h>
#include <blitz/array/stencil-et.h>

BZ_USING_NAMESPACE(blitz)

typedef blitz::Array<double,1> array_1;
typedef blitz::Array<double,2> array_2;
typedef blitz::Array<double,3> array_3;

/* Test slicing functionality for expressions. */

#define test_expr(d1,d2) BZTEST(all(abs(d1-d2)<1e-5));

// test with functors
class doubler {
public:
  double operator()(double x) const {return 2.0*x;}
  BZ_DECLARE_FUNCTOR(doubler);
};

class multiplier {
public:
  double operator()(double a, double b) const {return a*b;}
  BZ_DECLARE_FUNCTOR2(multiplier);
};


int main()
{
  array_2 a(5);
  a=tensor::i+10.*tensor::j;

  test_expr((2.*a)(1,Range::all()), array_2(2.*a)(1,Range::all()));
  test_expr(sin(a)(1,Range::all()), array_2(sin(a))(1,Range::all()));
  test_expr(where(a<5.,1./a,a)(1,Range::all()), 
	    array_2(where(a<5.,1./a,a))(1,Range::all()));
  // stencils can't be sliced down in rank, so integers are converted to unit ranges
  test_expr(Laplacian2D(sin(a))(1,Range::all()),
	    array_2(Laplacian2D(sin(a)))(Range(1,1),Range::all()));
  doubler d;
  multiplier m;
  test_expr(d(a)(1,Range::all()), array_2(d(a))(1,Range::all()));
  test_expr(m(a,a)(1,Range::all()), array_2(m(a,a))(1,Range::all()));

  //complete slicing to scalar
  BZTEST((2.*a)(1,1)==array_2(2.*a)(1,1));
  BZTEST((2.*a(1,Range::all())(1))==(array_1(2.*a(1,Range::all()))(1)));

  //more dims
  Array<double,5> big(2,3,4,5,6);
  big=tensor::i+10.*tensor::j+100*tensor::k+1000.*tensor::l+10000.*tensor::m;
  BZTEST(((2.*big)(1,2,3,4,5))==2.*big(1,2,3,4,5));
  
  test_expr((2.*big)(0,Range(0,1),1,Range::all(), Range(2,4,2)),
	    (Array<double,3>(2.*big(0,Range(0,1),1,Range::all(), Range(2,4,2)))));

  // slicing reductions, index remappings and index placeholder
  // expressions does not work.

  return 0;
}


