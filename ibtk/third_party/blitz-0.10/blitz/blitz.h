// -*- C++ -*-
/***************************************************************************
 * blitz/blitz.h      Includes all the important header files
 *
 * $Id$
 *
 * Copyright (C) 1997-2011 Todd Veldhuizen <tveldhui@acm.org>
 *
 * This file is a part of Blitz.
 *
 * Blitz is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation, either version 3
 * of the License, or (at your option) any later version.
 *
 * Blitz is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Blitz.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Suggestions:          blitz-devel@lists.sourceforge.net
 * Bugs:                 blitz-support@lists.sourceforge.net
 *
 * For more information, please see the Blitz++ Home Page:
 *    https://sourceforge.net/projects/blitz/
 *
 ***************************************************************************/

#ifndef BZ_BLITZ_H
#define BZ_BLITZ_H

/*
 * These symbols allow use of the IEEE and System V math libraries
 * (libm.a and libmsaa.a) on some platforms.
 */

#ifdef BZ_ENABLE_XOPEN_SOURCE
 #ifndef _ALL_SOURCE
  #define _ALL_SOURCE
 #endif
 #ifndef _XOPEN_SOURCE
  #define _XOPEN_SOURCE
 #endif
 #ifndef _XOPEN_SOURCE_EXTENDED
  #define _XOPEN_SOURCE_EXTENDED 1
 #endif
#endif

#include <blitz/bzconfig.h>
#include <blitz/compiler.h>          // Compiler-specific directives
#include <blitz/tuning.h>            // Performance tuning
#include <blitz/tau.h>               // Profiling

#ifdef BZ_HAVE_STL
  #include <string>
#endif

#ifdef BZ_HAVE_STD
  #include <iostream>
  #include <iomanip>
  #include <cstdio>                    // sprintf, etc.
  #include <cmath>
#else
  #include <iostream.h>
  #include <iomanip.h>
  #include <stdio.h>                   // sprintf, etc.
  #include <math.h>
#endif

#ifdef BZ_HAVE_COMPLEX
  #include <complex>
#endif

#define BZ_THROW                     // Needed in <blitz/numinquire.h>

// This macro is needed to pass template types to macros, since macros
// don't recognize <> as parentheses.
#define bzCC(...) __VA_ARGS__

BZ_NAMESPACE(blitz)

#ifdef BZ_HAVE_STD
 BZ_USING_NAMESPACE(std)
#endif

#ifdef BZ_GENERATE_GLOBAL_INSTANCES
 #define _bz_global
 #define BZ_GLOBAL_INIT(X)   =X
#else
 #define _bz_global extern
 #define BZ_GLOBAL_INIT(X)
#endif

/* Define types for indexing, depending on whether 64- or 32-bit
   indices are desired. There are separate typedefs for sizeType and
   indexType, because it might be useful to have possibility of arrays
   with 64-bit numbers of elements without paying the size overhead of
   making all dimensional indexes 64-bit.
 */
// Used for dimensional indexes (not implemented yet).
#ifdef BZ_FULLY64BIT
#warning 64-bit array dimensions not yet implemented
typedef ptrdiff_t indexType;
#else
typedef int indexType;
#endif
typedef size_t sizeType; // Used for memory indexing
typedef ptrdiff_t diffType; // Used for memory index differences, ie strides

// set the default padding policy
#ifdef BZ_PAD_ARRAYS
#define BZ_PADDING_DEFAULT paddedData
#else
#define BZ_PADDING_DEFAULT contiguousData
#endif

BZ_NAMESPACE_END

/*
 * Thread safety issues.  Compiling with -pthread under gcc, or -mt
 * under solaris, should automatically define _REENTRANT. Also have
 * support for OpenMP (which defines _OPENMP) or Windows thread
 * implementation.  The --enable-threadsafe configure option now
 * defines BZ_THREADSAFE. If this is defined but no thread support is
 * detected when compiling, we call #error below.
 */

/*
 * Which mutex implementation should be used for synchronizing
 * reference counts.  Options are Thread Building Block Atomics (which
 * is preferred over the others), pthreads, OpenMP, or Windows
 * threads. If we use TBB, the mutex macros are empty since it
 * requires no locking.
 */
#ifdef BZ_THREADSAFE
 #ifdef BZ_THREADSAFE_USE_TBB
  #include "tbb/atomic.h"
  #define BZ_THREADSAFE_USE_ATOMIC
  #define BZ_REFCOUNT_DECLARE(name) tbb::atomic<int> name;
 #else
  #define BZ_REFCOUNT_DECLARE(name) volatile int name;
  #if defined(_REENTRANT)
   #define BZ_THREADSAFE_USE_PTHREADS
  #elif defined (_OPENMP)
   #define BZ_THREADSAFE_USE_OPENMP
  #elif defined(_WIN32)
   #define BZ_THREADSAFE_USE_WINDOWS
  #else
   #error Blitz is configured with --enable-threadsafe, but no compiler thread support is found. Did you forget, e.g., "--pthread"?
  #endif
 #endif
#else
 #define BZ_REFCOUNT_DECLARE(name) int name;
#endif


#ifdef BZ_THREADSAFE_USE_PTHREADS
 #include <pthread.h>

 #define BZ_MUTEX_DECLARE(name)   mutable pthread_mutex_t name;
 #define BZ_MUTEX_INIT(name)      pthread_mutex_init(&name,NULL); mutexLocking_ = true;

 #define BZ_MUTEX_LOCK(name)      if (mutexLocking_) pthread_mutex_lock(&name);
 #define BZ_MUTEX_UNLOCK(name)    if (mutexLocking_) pthread_mutex_unlock(&name);
 #define BZ_MUTEX_DESTROY(name)   pthread_mutex_destroy(&name);
#elif defined (BZ_THREADSAFE_USE_WINDOWS)
 // Include Windows.h header in case user has not already done so.
 // Disable Windows min/max macro definitions
 #define NOMINMAX
 #include <Windows.h>

 #define BZ_MUTEX_DECLARE(name)   mutable CRITICAL_SECTION name;
 #define BZ_MUTEX_INIT(name)      ::InitializeCriticalSection(&name);  mutexLocking_ = true;
 #define BZ_MUTEX_LOCK(name)      if (mutexLocking_) ::EnterCriticalSection(&name);
 #define BZ_MUTEX_UNLOCK(name)    if (mutexLocking_) ::LeaveCriticalSection(&name);
 #define BZ_MUTEX_DESTROY(name)   ::DeleteCriticalSection(&name);
#elif defined (BZ_THREADSAFE_USE_OPENMP)
 #include <omp.h>

 #define BZ_MUTEX_DECLARE(name)   mutable omp_lock_t name;
 #define BZ_MUTEX_INIT(name)      omp_init_lock(&name); mutexLocking_ = true;
 #define BZ_MUTEX_LOCK(name)      if (mutexLocking_) omp_set_lock(&name);
 #define BZ_MUTEX_UNLOCK(name)    if (mutexLocking_) omp_unset_lock(&name);
 #define BZ_MUTEX_DESTROY(name)   omp_destroy_lock(&name);
#else
 #define BZ_MUTEX_DECLARE(name)
 #define BZ_MUTEX_INIT(name)
 #define BZ_MUTEX_LOCK(name)
 #define BZ_MUTEX_UNLOCK(name)
 #define BZ_MUTEX_DESTROY(name)
#endif

#include <blitz/bzdebug.h>           // Debugging macros

#endif // BZ_BLITZ_H
