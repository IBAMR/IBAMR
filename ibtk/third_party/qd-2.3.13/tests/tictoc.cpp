/*
 * tests/tictoc.cpp
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2006
 *
 * Contains function used for timing.
 */

#include "tictoc.h"


#ifndef _WIN32

#ifdef HAVE_CLOCK_GETTIME

#ifdef CLOCK_HIGHRES
#define SAMPLED_CLOCK CLOCK_HIGHRES
#else
#define SAMPLED_CLOCK CLOCK_REALTIME
#endif

void tic(tictoc *tv) {
  if (clock_gettime(SAMPLED_CLOCK, tv)) 
    tv->tv_sec = tv->tv_nsec = -1;
}

double toc(tictoc *tv) {
  struct timespec tv2;

  if (clock_gettime(SAMPLED_CLOCK, &tv2)) 
    tv2.tv_sec = tv2.tv_nsec = -1;

  double  sec = static_cast<double>(tv2.tv_sec - tv->tv_sec);
  double nsec = static_cast<double>(tv2.tv_nsec - tv->tv_nsec);

  return (sec + 1.0e-9 * nsec);
}
#else

#ifdef HAVE_GETTIMEOFDAY

void tic(tictoc *tv) {
  gettimeofday(tv, 0L);
}

double toc(tictoc *tv) {
  tictoc tv2;

  gettimeofday(&tv2, 0L);
  double  sec = static_cast<double>(tv2.tv_sec - tv->tv_sec);
  double usec = static_cast<double>(tv2.tv_usec - tv->tv_usec);

  return (sec + 1.0e-6 * usec);
}

#else
// Fall back to C/C++ low resolution time function.

void tic(tictoc *tv) {
  time(tv);
}

double toc(tictoc *tv) {
  tictoc tv2;
  time(&tv2);
  return difftime(tv2, *tv);
}

#endif

#endif

#else

// Windows.

void tic(tictoc *tv) {
  *tv = GetTickCount();
}

double toc(tictoc *tv) {
  tictoc tv2;
  tv2 = GetTickCount();
  return 1.0e-3 * (tv2 - *tv);
}

#endif
