// -*- C++ -*-
/***************************************************************************
 * blitz/Timer.h        Timer class, for use in benchmarking
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

// This class is not portable to non System V platforms.
// It will need to be rewritten for Windows, NT, Mac.
// NEEDS_WORK

#ifndef BZ_TIMER_H
#define BZ_TIMER_H

#ifndef BZ_BLITZ_H
 #include <blitz/blitz.h>
#endif

#ifdef BZ_HAVE_RUSAGE
 #include <sys/resource.h>
#else
 #include <time.h>
#endif

BZ_NAMESPACE(blitz)

#ifndef BZ_HAVE_LIBPAPI

class Timer {

public:
    Timer()
    {
        state_ = uninitialized;
    }

    void start()
    {
        state_ = running;
        t1_ = systemTime();
    }

    void stop()
    {
        t2_ = systemTime();
        BZPRECONDITION(state_ == running);
        state_ = stopped;
    }

/* Compaq cxx compiler in ansi mode cannot print out long double type! */
#if defined(__DECCXX)
    double elapsed() const
#else
    /** Return elapsed time in seconds. */
    double elapsed() const
#endif
    {
        BZPRECONDITION(state_ == stopped);
        return (t2_ - t1_)/1e6;
    }

  long long instr() const { return 0; };
  long long flops() const { return 0; };

  static const string& indep_var() { return ivar_; };

private:
    Timer(Timer&) { }
    void operator=(Timer&) { }

    long int systemTime()
    {
#ifdef BZ_HAVE_RUSAGE
        getrusage(RUSAGE_SELF, &resourceUsage_);
        long int sec = resourceUsage_.ru_utime.tv_sec
            + resourceUsage_.ru_stime.tv_sec;
        long int usec  = resourceUsage_.ru_utime.tv_usec
            + resourceUsage_.ru_stime.tv_usec;
        return sec*1000000+usec;
#else
        return 1000000*clock() / (long int) CLOCKS_PER_SEC;
#endif
    }

    enum { uninitialized, running, stopped } state_;

  static const string ivar_;

#ifdef BZ_HAVE_RUSAGE
    struct rusage resourceUsage_;
#endif

    long int t1_, t2_;
};

#else

// implementation using PAPI performance counters

#include <papi.h>

class Timer {

public:
    Timer()
    {
      // maybe overhead is less from just reading counters
      if(PAPI_start_counters((int*)Events, nevents)!=PAPI_OK) {
	cerr << "Error starting counters\n";
      }
      state_ = uninitialized;
    }
  ~Timer()
    {
      PAPI_stop_counters(counters_.data(), nevents);
    }

    void start()
    {
        state_ = running;
	// this resets the counters
        PAPI_read_counters(counters_.data(), nevents);
    }

    void stop()
    {
        PAPI_read_counters(counters_.data(), nevents);
	BZPRECONDITION(state_ == running);
	state_ = stopped;
    }

  /** since we don't know the clock frequency of the processor, we
      instead output "flops/clock cycle" which seems like a better
      measure of code performance and not machine performance. */
    long long elapsed() const
    {
        BZPRECONDITION(state_ == stopped);
        return counters_[0];
    }

  long long instr() const { return counters_[1]; };
  long long flops() const { return counters_[2]; };

  static const string& indep_var() { return ivar_; };

private:
    Timer(Timer&) { }
    void operator=(Timer&) { }

    enum { uninitialized, running, stopped } state_;

  static const int nevents=3;
  static const int Events[nevents];
  static const string ivar_;

  TinyVector<long long, nevents> counters_;
};

#endif




BZ_NAMESPACE_END

#endif // BZ_TIMER_H

