// -*- C++ -*-
/***************************************************************************
 * blitz/bench.h      Benchmark classes
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

#ifndef BZ_BENCH_H
#define BZ_BENCH_H

#include <blitz/matrix.h>
#include <blitz/timer.h>

#if defined(BZ_HAVE_STD)
#include <cmath>
#else
#include <math.h>
#endif

BZ_NAMESPACE(blitz)

// Forward declaration
template<typename P_parameter = unsigned>
class BenchmarkImplementation;


// Declaration of class Benchmark<T>
// The template parameter T is the parameter type which is varied in
// the benchmark.  Typically T will be an unsigned, and will represent
// the length of a vector, size of an array, etc.

template<typename P_parameter = unsigned>
class Benchmark {

public:
    typedef P_parameter T_parameter;

    Benchmark(unsigned numImplementations);

    ~Benchmark();

    void addImplementation(BenchmarkImplementation<T_parameter>*
        implementation);

    void run(ostream& log = cout);

    double getMflops(unsigned implementation, unsigned setting) const;

    double getRate(unsigned implementation, unsigned setting) const;

    void saveMatlabGraph(const char* filename) const;
    void savePylabGraph(const char* filename) const;

public:
    // Virtual functions

    virtual const char* description() const
    { return ""; }

    virtual const char* parameterDescription() const
    { return "Vector length"; }

    virtual unsigned numParameterSettings() const
    { return 19; }

    virtual T_parameter getParameterSetting(unsigned i) const
    { return BZ_MATHFN_SCOPE(pow)(10.0, (i+1)/4.0); }

    virtual long getIterationSetting(unsigned i) const
    { return 1000000L / getParameterSetting(i); }

private:
    Benchmark(const Benchmark<P_parameter>&) { }
    void operator=(const Benchmark<P_parameter>&) { }

    enum { uninitialized, initialized, running, done } state_;

    unsigned numImplementations_;
    unsigned numStoredImplementations_;

    BenchmarkImplementation<T_parameter>** implementations_;

    Matrix<double,RowMajor> rates_;       // Iterations per second array
    Matrix<double,RowMajor> Mflops_;
};

template<typename P_parameter>
class BenchmarkImplementation {

public:
    typedef P_parameter T_parameter;

    virtual void initialize(P_parameter parameter) { }

    virtual void done() { }

    virtual const char* implementationName() const
    { return ""; }

    virtual void run(long iterations) = 0;

    virtual void runOverhead(long iterations)
    {
        for (long i=0; i < iterations; ++i)
        {
        }
    };

    virtual void tickle() { }

    virtual long flopsPerIteration() const
    { return 0; }
};

BZ_NAMESPACE_END

#include <blitz/bench.cc>

#endif // BZ_BENCH_H
