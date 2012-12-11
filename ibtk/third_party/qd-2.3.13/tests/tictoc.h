/*
 * tests/timer.h
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains function used for timing.
 */

#ifndef TICTOC_H__
#define TICTOC_H__

#include "config.h"

#ifdef _WIN32

#include <windows.h>
typedef DWORD tictoc;
#else

#ifdef HAVE_CLOCK_GETTIME
#include <time.h>
typedef struct timespec tictoc;
#else

#ifdef HAVE_GETTIMEOFDAY
#include <sys/time.h>
typedef struct timeval tictoc;
#else
#include <ctime>
typedef time_t tictoc;
#endif

#endif

#endif

void   tic(tictoc *tv);   /* start timing. */
double toc(tictoc *tv);   /* stop  timing. */

#endif
