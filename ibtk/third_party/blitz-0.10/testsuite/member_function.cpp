#include "testsuite.h"
#include <blitz/array.h>

//#undef BZ_FAST_COMPILE

BZ_USING_NAMESPACE(blitz)

// test that member function functors work. (There was a bug where
// they would use an out-of-scope temporary copy of the class itself
// if BZ_ETPARM was pass-by-value.)

class testclass {
private:
  testclass(const testclass&) {
    BZTEST(0);
  };

public:
  int d_;

  testclass(int d) : d_(d) {};

  int func(int d) const {cout << "In func, this is " << this << endl; return d_*d; };

  BZ_DECLARE_MEMBER_FUNCTION(testclass, func);
};


int main()
{
  Array<int,1> a(3), b(3);
  a=1,2,3;

  testclass t(2);
  cout << "&t is " << &t << endl;
  // this expression should not require copying t.
  b=t.func(2*a);

  return 0;
}

