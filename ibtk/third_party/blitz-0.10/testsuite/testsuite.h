#ifndef TESTSUITE_H
#define TESTSUITE_H

#ifdef BZ_DEBUG_H
 #error <testsuite/testsuite.h> must be included before any other header files.
#endif

#define BZ_TESTSUITE

#include <blitz/blitz.h>

#define BZTEST(x)    { if(!(x)) { BZ_STD_SCOPE(cout) << "Test failed: " << __FILE__ << " line " << __LINE__ << BZ_STD_SCOPE(endl); BZ_STD_SCOPE(exit)(1); } }

#endif // TESTSUITE_H
