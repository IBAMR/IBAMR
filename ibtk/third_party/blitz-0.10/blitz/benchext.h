// -*- C++ -*-
/***************************************************************************
 * blitz/benchext.h      BenchmarkExt classes (Benchmarks with external
 *                       control)
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

#ifndef BZ_BENCHEXT_H
#define BZ_BENCHEXT_H

//#ifndef BZ_MATRIX_H
#include <blitz/array.h>
//#endif
 #include <vector>

#ifndef BZ_TIMER_H
 #include <blitz/timer.h>
#endif

#include <math.h>
#include <string>

// NEEDS_WORK: replace use of const char* with <string>, once standard
// library is widely supported.

BZ_NAMESPACE(blitz)

// Declaration of class BenchmarkExt<T>
// The template parameter T is the parameter type which is varied in
// the benchmark.  Typically T will be an unsigned, and will represent
// the length of a vector, size of an array, etc.

template<typename P_parameter = unsigned>
class BenchmarkExt {

public:
    typedef P_parameter T_parameter;

    BenchmarkExt(const char* description, int numImplementations);

    ~BenchmarkExt();

    void setNumParameters(int numParameters);
  void setParameterVector(Array<T_parameter,1> parms);
    void setParameterDescription(const char* string);
  void setIterations(Array<long,1> iters);
  void setOpsPerIteration(Array<double,1> flopsPerIteration);
  void setDependentVariable(const char* string);

    void beginBenchmarking();

    void beginImplementation(const char* description);
    bool doneImplementationBenchmark() const;
  const string& currentImplementation() const;

    T_parameter getParameter() const;
    long        getIterations() const;

    inline void start();
    inline void stop();

    void startOverhead();
    void stopOverhead();
    void skip();


    void endImplementation();

    void endBenchmarking();
    double getMflops(unsigned implementation, unsigned parameterNum) const;
    double getinstrperc(int implementation, int parameterNum) const;
    double getflopsperc(int implementation, int parameterNum) const;

    void saveMatlabGraph(const char* filename, const char* graphType="semilogx") const;
    void savePylabGraph(const char* filename, const char* graphType="semilogx") const;

protected:
    BenchmarkExt(const BenchmarkExt<P_parameter>&) { }
    void operator=(const BenchmarkExt<P_parameter>&) { }

    enum { initializing, benchmarking, benchmarkingImplementation,
       running, runningOverhead, done } state_;

    unsigned numImplementations_;
    unsigned implementationNumber_;

  std::string description_;
  std::vector<std::string> implementationDescriptions_;

  Array<double,2> times_;       // Elapsed time
  Array<long long,2> instr_;       // instructions according to timer
  Array<long long,2> flops_;       // flops according to timer

  Array<T_parameter,1> parameters_;
  Array<long,1> iterations_;
  Array<double,1> flopsPerIteration_;

    Timer timer_;

  std::string parameterDescription_;
  std::string depvar_;
  double timerconversion_;

    unsigned numParameters_;
    unsigned parameterNumber_;
};

BZ_NAMESPACE_END

#include <blitz/benchext.cc>

#endif // BZ_BENCHEXT_H
