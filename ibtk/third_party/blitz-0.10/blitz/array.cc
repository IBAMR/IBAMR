
#ifndef BZ_ARRAY_CC
#define BZ_ARRAY_CC

// not-quite-declarations-of-array files

#include <blitz/array/iter.h>       // Array iterators
#include <blitz/array/fastiter.h>   // Fast Array iterators (for et)
#include <blitz/array/asexpr.cc>       // ET type functions
#include <blitz/array/expr.cc>       // Array expression objects
#include <blitz/range.cc>       // range objects
#include <blitz/array/methods.cc>   // Member functions
#include <blitz/array/ops.cc>       // Assignment operators
#include <blitz/array/io.cc>        // Output formatting
#include <blitz/array/interlace.cc> // Allocation of interlaced arrays
#include <blitz/array/resize.cc>    // Array resize, resizeAndPreserve
#include <blitz/array/slicing.cc>   // Slicing and subarrays
#include <blitz/array/cycle.cc>     // Cycling arrays
#include <blitz/array/complex.cc>   // Special support for complex arrays
#include <blitz/array/et.h>         // Expression templates
#include <blitz/array/reduce.h>     // Array reduction expression templates
#include <blitz/array/zip.h>        // Zipping multicomponent types
#include <blitz/array/where.h>      // where(X,Y,Z)
#include <blitz/array/indirect.h>   // Indirection
#include <blitz/array/stencils.h>   // Stencil objects
#include <blitz/tinymat2.cc>        // TinyMatrix implementations
#endif
