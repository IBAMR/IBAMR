#include "testsuite.h"
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

int main(int argc, char** argv)
{
  // test 64-bit compatibility if enabled
  if(sizeof(sizeType)>4) {

    // wrap in try so we can catch allocation failure if test platform
    // has insufficient memory. To minimize memory reqs we use arrays
    // of char, but we still need at least 4GB. (We want it to be a
    // contiguous array so the loops are coalesced.)
    Array<unsigned char,4> test(contiguousArray);
    try {
      test.resize(shape(256,256,256,256));
    }
    catch (std::bad_alloc&) {
      cerr << "Allocation failure - not enough memory to test 64-bit compatibility\n";
      return 0;
    }

    // we MUST do a full assignment, because if we use a subarray the
    // offset is wrt the subarray. But we avoid using some complicated
    // index expression since that makes it very slow.
    test = 42;
    unsigned char crap= test(255,255,255,255);
    cerr << int(crap) << endl;
    BZTEST(test(255,255,255,255)==42);
    test = test + test;
    BZTEST(test(255,255,255,255)==84);
    test.free();

#ifdef BZ_FULLY64BIT
    // test arrays with long dimensions
    Array<char, 1> test2(1<<32);
    test2 = 1;
    BZTEST(sum(test2)==1<<32);
    test2.free();
#endif
  }
  else {
    cerr << "Not a 64-bit machine, can't test 64-bit compatibility.\n";
  }

  return 0;
}

