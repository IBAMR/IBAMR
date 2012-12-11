#define BZ_GENERATE_GLOBAL_INSTANCES

/*
 * This module contains instances of global variables and objects.
 * The preprocessor symbol _bz_global normally expands to "extern",
 * but the declaration of BZ_GENERATE_GLOBAL_INSTANCES causes it
 * to expand to "" for this source file.  This ensures that globals
 * will have one unique instance, in libblitz.a.
 */

#include <blitz/blitz.h>          // Basic Blitz++ includes
#include <blitz/array.h>          // Dummy array object (_dummyArray)
#include <blitz/indexexpr.h>      // Tensor index objects
#include <blitz/bzdebug.h>        // Test suite globals
#include <blitz/timer.h>          // Timer instances
#include <random/mtparam.cc>      // MT19937 parameter sets
#include <string>

BZ_NAMESPACE(blitz);

// used to avoid loop elimination in the benchmarks
void sink() {};

char _blitz_id[] = BZ_PACKAGE_STRING;

#ifdef BZ_HAVE_LIBPAPI
#include <papiStdEventDefs.h>
const int Timer::Events[] =  {PAPI_TOT_CYC, PAPI_TOT_INS, PAPI_FP_INS};
const string Timer::ivar_="c";
#else
const string Timer::ivar_="s";
#endif

BZ_NAMESPACE_END
