/*
 * Bug found by Derrick Bass <derrick@tapir.caltech.edu>.  In
 * ../blitz/arrayeval.cc, evaluateWithStackTraversal1(), there
 * was an error in unrolling loops with common strides which caused
 * the wrong loop bounds to be chosen.  As a result, more elements
 * were overwritten than were supposed to be.  The appropriate
 * fix was changing this line:
 *
 *             int n1 = (ubound & 3) * commonStride;
 *
 * to
 *
 *             int n1 = (length(firstRank) & 3) * commonStride;
 *
 * -- TV 980402
 */

#include "testsuite.h"

#include <blitz/tuning.h>

#define BZ_ARRAY_STACK_TRAVERSAL_UNROLL
#define BZ_ARRAY_EXPR_USE_COMMON_STRIDE

#include <blitz/array.h>

BZ_USING_NAMESPACE(std)
BZ_USING_NAMESPACE(blitz)

int assignCount = 0;

class ddouble {
  double d;

public:
  ddouble() : d(0) { }
  ddouble(double _d) : d(_d) { }
  
  ddouble& operator=(double _d) { 
    d = _d; 
    return *this; }
  ddouble& operator=(ddouble _d) { 
    d = _d.d; 
    ++assignCount;
    return *this; }
  
  operator double() { return d; }
};

int main() 
{
  Range fr(0,11);
  Range bw2(0,2);
  
  Array<ddouble, 2> D(bw2,fr);

  D(bw2, 11) = 0;

  BZTEST(assignCount == 3);
  return 0;
}

